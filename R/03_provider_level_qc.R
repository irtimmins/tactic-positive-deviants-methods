# 03  provider-level_qc (run once, then curate)
# -----------------------------------------------------------------------------
# Reconciles three sources into harmonised 5-digit diagnosing and treating sites:
#   cohort tabulations  -> the codes that must be placed (with volume)
#   ODS API             -> identity: trust vs site, active vs closed, parent trust
#   curated Excel       -> eligibility: Bowel_ca_surgery, merger colour
# Diagnosing and treating are both 5-digit sites. A 3-char trust code is placed
# only when its trust has a single hospital site; otherwise the patients are
# unattributable and excluded. Treating sites are kept only if they do bowel
# cancer surgery. Outputs feed 04.
#
# This step needs the ODS API and the curated Excel, so run it in the secure
# environment. Review the two printed tables and the code_aliases before trusting
# the crosswalks; ODS lookups are cached so reruns are cheap.

library(tidyverse)
library(readxl)
library(httr2)

source("R/01_config.R")

# Excel column names (adjust to match your file)
xl_site   <- "Hospital_site_code"
xl_trust  <- "trust_nacs"
xl_trustn <- "Trust_Name"
xl_siten  <- "Hospital_Name"
xl_bowel  <- "Bowel_ca_surgery"
xl_colour <- "trust_nacs_colour"

vol_year_guide <- 10   # per-year guidance only; the window cut is applied in 04
# site_min_vol and dominance_share come from R/01_config.R (shared with the
# calibration helper); calibrate them there.

# code unification: fold a code that is really the same hospital under another
# code into its canonical 5-digit site - elective surgical hubs, old or merged
# codes, or a trust code you want to force to a specific site. Applied BEFORE
# colon-site detection (so the merged codes count as one site) and written into
# the crosswalk (so 04 remaps the raw data). Seed it from nonstandard_codes.csv
# and the calibration helper, then rerun.
code_aliases <- tribble(
  ~alt_code, ~canonical_code, ~note
  # "E0A3H", "RXH01", "Royal Sussex elective hub -> Royal Sussex main site"
  # "RV820", "R1K01", "old Northwick Park code -> current site"
)
alias_lu    <- setNames(code_aliases$canonical_code, code_aliases$alt_code)
canon_alias <- function(x) { y <- unname(alias_lu[x]); ifelse(is.na(y), x, y) }

# treating sites confirmed (or rejected) as bowel-cancer surgery by clinical
# review, overriding the Excel flag. Fixes code-mismatch false negatives, where a
# real colorectal unit failed only because its cohort code differs from the Excel.
bowel_confirm <- c(
  "RRK99", "RX1CC", "RWP01", "RM317", "RH5A8", "R0D02", "RNN62", "Z1R8K", "R0D01",
  "RYJ02", "R0B01", "RBN51", "E0Z3F", "RH5O4", "RH880", "T9F5R", "E0A3H", "Q8U9Z",
  "R0B0Q", "R8T1Q", "I3W1A", "M7L2T", "E5E1O", "C2P9J", "I1Z8O", "F3E3F", "D8W9O", "B0D8Z")
bowel_reject <- c("RN543", "V8O2H", "RTFFS")   # reviewed: not bowel-surgery sites

# 1. observed codes with per-year volume -------------------------------------
coh <- readRDS(in_rds)
per_year <- function(d, code_col, only_resected = FALSE) {
  if (only_resected) d <- filter(d, had_surgery)
  d %>% transmute(code = trimws(as.character(.data[[code_col]])), ydiag) %>%
    filter(!is.na(code), code != "") %>%
    count(code, ydiag, name = "n") %>%
    group_by(code) %>%
    summarise(n_total = sum(n), n_years = n_distinct(ydiag),
              mean_per_year = round(sum(n) / n_distinct(ydiag), 1), .groups = "drop")
}
diag_obs  <- per_year(coh, "diag_hosp")
treat_obs <- per_year(coh, "SITETRET", only_resected = TRUE)

# merged view for colon-site detection: fold aliased codes into their canonical
# site so two codes for one hospital count as one site, not two competing ones
diag_grp <- diag_obs %>% mutate(code = canon_alias(code)) %>%
  group_by(code) %>% summarise(n_total = sum(n_total), .groups = "drop")

# 2. curated Excel site list -------------------------------------------------
site <- read_excel(site_xlsx) %>%
  transmute(site_code  = trimws(.data[[xl_site]]),
            trust_code = trimws(.data[[xl_trust]]),
            trust_name = .data[[xl_trustn]],
            site_name  = .data[[xl_siten]],
            bowel_ca_surgery = suppressWarnings(as.integer(.data[[xl_bowel]])),
            colour     = tolower(trimws(.data[[xl_colour]]))) %>%
  mutate(site_code = if_else(is.na(site_code) | site_code == "",
                             str_match(site_name, "\\(([A-Z0-9]{3,5})\\)\\s*$")[, 2], site_code),
         colour    = if_else(colour %in% c("", "blank", NA), "unchanged", colour),
         dissolved = colour == "red") %>%
  filter(!is.na(site_code))

# some hospital names carry an alternate code in brackets, e.g. "Northwick Park
# Hospital (RV820)" - register it as an alias row so a cohort code that uses the
# bracketed form still finds the trust / bowel-surgery flag.
aliases <- site %>%
  mutate(bracket = str_match(site_name, "\\(([A-Z0-9]{5})\\)\\s*$")[, 2]) %>%
  filter(!is.na(bracket), bracket != site_code) %>%
  mutate(site_code = bracket) %>% select(-bracket)
site <- bind_rows(site, aliases)

# one row per site code (your Excel can have several, e.g. recoloured duplicates);
# prefer the active, bowel-surgery row so the kept flag is the favourable one
site_lu <- site %>%
  arrange(dissolved, desc(coalesce(bowel_ca_surgery, 0L))) %>%
  distinct(site_code, .keep_all = TRUE)

# single-colon-site detection is done after the ODS lookup (it needs ODS parents
# to group observed sites under their trust); see the trust_colon block below.

# 3. ODS lookup (cached) -----------------------------------------------------
ods_base <- "https://directory.spineservices.nhs.uk/ORD/2-0-0/organisations/"
role_labels <- c(RO198 = "NHS trust", RO197 = "NHS trust site",
                 RO157 = "non-NHS organisation",
                 RO172 = "independent sector site", RO182 = "independent sector")
as_list_of <- function(x) if (is.null(x)) list() else if (!is.null(names(x))) list(x) else x
lookup_ods <- function(code) {
  empty <- tibble(code = code, ods_name = NA_character_, status = NA_character_,
                  record_class = NA_character_, primary_role = NA_character_,
                  parent_code = NA_character_, postcode = NA_character_,
                  last_change = NA_character_)
  resp <- request(paste0(ods_base, code)) |>
    req_user_agent("tactic-provider-qc") |> req_retry(max_tries = 3) |>
    req_error(is_error = \(r) FALSE) |> req_perform()
  if (resp_status(resp) != 200) return(empty)
  org <- resp_body_json(resp)$Organisation
  if (is.null(org)) return(empty)
  roles   <- as_list_of(org$Roles$Role)
  prim    <- keep(roles, \(r) isTRUE(r$primaryRole))
  prim_id <- if (length(prim)) prim[[1]]$id else NA_character_
  rels    <- as_list_of(org$Rels$Rel)
  parent  <- rels |>
    keep(\(r) identical(r$Status, "Active") && !is.null(r$Target$OrgId$extension) &&
           r$Target$OrgId$extension != code) |>
    map_chr(\(r) r$Target$OrgId$extension) |> head(1)
  `%||%` <- rlang::`%||%`
  tibble(code = code, ods_name = org$Name %||% NA_character_,
         status = org$Status %||% NA_character_,
         record_class = recode(org$orgRecordClass %||% NA_character_,
                               RC1 = "trust", RC2 = "site", .default = NA_character_),
         primary_role = recode(prim_id, !!!role_labels, .default = prim_id),
         parent_code  = if (length(parent)) parent else NA_character_,
         postcode     = org$GeoLoc$Location$PostCode %||% NA_character_,
         last_change  = org$LastChangeDate %||% NA_character_)
}

all_codes <- unique(c(diag_obs$code, treat_obs$code, site$site_code, site$trust_code))
all_codes <- all_codes[!is.na(all_codes) & all_codes != ""]
ods_cache_path <- file.path(qc_dir, "ods_cache.csv")
if (file.exists(ods_cache_path)) {
  ods <- read_csv(ods_cache_path, col_types = cols(.default = "c"))
  todo <- setdiff(all_codes, ods$code)
} else { ods <- tibble(); todo <- all_codes }
if (length(todo) > 0) {
  ods <- bind_rows(ods, todo %>% map(possibly(lookup_ods, otherwise = NULL),
                                     .progress = TRUE) %>% list_rbind()) %>%
    distinct(code, .keep_all = TRUE)
  write_csv(ods, ods_cache_path)
}

is_site_code  <- function(code, rc) (!is.na(rc) & rc == "site")  | (is.na(rc) & nchar(code) == 5)
is_trust_code <- function(code, rc) (!is.na(rc) & rc == "trust") | (is.na(rc) & nchar(code) == 3)

# one colon site per trust -----------------------------------------------------
# group the observed 5-digit diagnosing sites under their ODS parent trust, and
# look at how the trust's colon volume is spread across its sites. A trust code
# is mapped to a single site only when one site clearly dominates; trusts with
# colon activity genuinely spread across sites stay unattributable.
obs_sites <- diag_grp %>% filter(nchar(code) == 5) %>%
  left_join(ods %>% select(code, parent_code), by = "code") %>%
  filter(!is.na(parent_code)) %>%
  rename(trust = parent_code)

trust_colon <- obs_sites %>%
  arrange(trust, desc(n_total)) %>%
  group_by(trust) %>%
  summarise(n_sites          = n(),
            n_sites_real     = sum(n_total >= site_min_vol),
            trust_site_vol   = sum(n_total),
            top_site         = first(code),
            top_vol          = first(n_total),
            second_site      = nth(code, 2),
            second_vol       = coalesce(nth(n_total, 2), 0L),
            top_share        = first(n_total) / sum(n_total),
            .groups = "drop") %>%
  mutate(classification = case_when(
    top_vol < site_min_vol       ~ "no real colon site",
    top_share >= dominance_share ~ "single/dominant colon site",
    TRUE                         ~ "ambiguous multi-site"),
    colon_site = if_else(classification == "single/dominant colon site",
                         top_site, NA_character_))

# the trust's own (3-char) cohort volume, i.e. what mapping recovers / excluding loses
trust_code_vol <- diag_obs %>% filter(nchar(code) == 3) %>%
  transmute(trust = code, trust_code_vol = n_total)

# readable names for the evidence tables
name_tbl <- bind_rows(site_lu %>% transmute(code = site_code, name = site_name),
                      ods %>% transmute(code, name = ods_name)) %>%
  filter(!is.na(name), name != "") %>% distinct(code, .keep_all = TRUE)

trust_colon <- trust_colon %>%
  left_join(trust_code_vol, by = "trust") %>%
  left_join(name_tbl %>% rename(top_site = code, top_name = name), by = "top_site") %>%
  left_join(name_tbl %>% rename(second_site = code, second_name = name), by = "second_site")
write_csv(trust_colon, file.path(qc_dir, "trust_colon_site.csv"))

# per-trust per-site volumes: the raw material for calibrating site_min_vol and
# dominance_share, and for spotting two codes that are really one hospital
trust_site_volumes <- obs_sites %>% left_join(name_tbl, by = "code") %>%
  arrange(trust, desc(n_total)) %>% transmute(trust, site = code, site_name = name, n_total)
write_csv(trust_site_volumes, file.path(qc_dir, "trust_site_volumes.csv"))

# non-standard codes (not a 3-char trust or 5-char site code) - elective surgical
# hubs and recoded sites - the candidates for code_aliases
all_obs <- bind_rows(diag_obs %>% transmute(code, n = n_total, side = "diag"),
                     treat_obs %>% transmute(code, n = n_total, side = "treat"))
nonstandard <- all_obs %>%
  filter(!grepl("^R[A-Z0-9]{4}$", code), !grepl("^R[A-Z0-9]{2}$", code)) %>%
  left_join(name_tbl, by = "code") %>%
  group_by(code, name) %>%
  summarise(diag_vol = sum(n[side == "diag"]), treat_vol = sum(n[side == "treat"]), .groups = "drop") %>%
  arrange(desc(diag_vol + treat_vol))
write_csv(nonstandard, file.path(qc_dir, "nonstandard_codes.csv"))

trust_colon_lu <- setNames(trust_colon$colon_site, trust_colon$trust)

# resolve an observed diagnosing code to a 5-digit site, or NA if unattributable
resolve_site <- function(code) {
  if (code %in% names(alias_lu)) return(unname(alias_lu[code]))
  o  <- ods %>% filter(code == !!code)
  rc <- if (nrow(o)) o$record_class[1] else NA_character_
  if (is_site_code(code, rc))  return(code)
  if (is_trust_code(code, rc)) {
    cs <- trust_colon_lu[[code]]
    return(if (is.null(cs) || is.na(cs)) NA_character_ else cs)
  }
  NA_character_
}

# 4. diagnosing crosswalk (5-digit site) -------------------------------------
diag <- diag_obs %>% left_join(ods, by = "code") %>%
  mutate(resolved_site = map_chr(code, resolve_site),
         code_type = case_when(is_site_code(code, record_class) ~ "site",
                               is_trust_code(code, record_class) ~ "trust", TRUE ~ "unknown")) %>%
  mutate(canonical_code = resolved_site) %>%
  left_join(site_lu %>% select(canonical_code = site_code, site_name, trust_code,
                               colour, dissolved_excel = dissolved), by = "canonical_code") %>%
  mutate(canonical_name = coalesce(site_name, ods_name),
         is_dissolved = coalesce(dissolved_excel, status == "Inactive"),
         exclude_reason = case_when(!is.na(canonical_code) ~ NA_character_,
                                    code_type == "trust" ~ "multi-site trust (unattributable)",
                                    TRUE ~ "unknown / unresolved code"),
         needs_review = is.na(canonical_code) | is.na(canonical_name)) %>%
  select(raw_code = code, canonical_code, canonical_name, code_type, ods_status = status,
         is_dissolved, exclude_reason, n_total, mean_per_year, trust_code, needs_review)

diag_canon <- diag %>% filter(!is.na(canonical_code)) %>%
  group_by(canonical_code, canonical_name, trust_code) %>%
  summarise(accepted_codes = paste(sort(unique(raw_code)), collapse = "; "),
            n_codes = n_distinct(raw_code), n_total = sum(n_total),
            mean_per_year = round(sum(n_total) / max(mean_per_year, na.rm = TRUE), 1),
            any_dissolved = any(is_dissolved, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_total))

write_csv(diag, file.path(qc_dir, "diagnosing_providers.csv"))
write_csv(diag_canon, file.path(qc_dir, "diagnosing_providers_canonical.csv"))
write_csv(diag %>% filter(!is.na(canonical_code)) %>%
            transmute(raw_code, canonical_code, canonical_name), diag_xwalk_csv)
# diagnosing include: resolved sites that are not dissolved
write_csv(diag_canon %>% filter(!any_dissolved) %>% transmute(canonical_code), diag_include_csv)

# 5. treating crosswalk (5-digit site, bowel surgery) ------------------------
treat <- treat_obs %>% left_join(ods, by = "code") %>%
  left_join(site_lu %>% select(code = site_code, trust_code, site_name, bowel_ca_surgery,
                               colour, dissolved_excel = dissolved), by = "code") %>%
  mutate(canonical_code = canon_alias(code),
         canonical_name = coalesce(site_name, ods_name),
         trust_code = coalesce(trust_code, parent_code),
         is_dissolved = coalesce(dissolved_excel, status == "Inactive"),
         bowel_ca_surgery = case_when(
           code %in% bowel_confirm ~ 1L,
           code %in% bowel_reject  ~ 0L,
           TRUE ~ coalesce(bowel_ca_surgery, 0L)),
         needs_review = is.na(canonical_name) | is.na(bowel_ca_surgery)) %>%
  select(raw_code = code, canonical_code, canonical_name, trust_code, bowel_ca_surgery,
         ods_status = status, is_dissolved, n_total, mean_per_year, needs_review)

treat_canon <- treat %>%
  group_by(canonical_code, canonical_name, trust_code) %>%
  summarise(accepted_codes = paste(sort(unique(raw_code)), collapse = "; "),
            n_codes = n_distinct(raw_code), n_total = sum(n_total),
            mean_per_year = round(sum(n_total) / max(mean_per_year, na.rm = TRUE), 1),
            bowel_ca_surgery = max(bowel_ca_surgery),
            any_dissolved = any(is_dissolved, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_total))

write_csv(treat, file.path(qc_dir, "treating_providers.csv"))
write_csv(treat_canon, file.path(qc_dir, "treating_providers_canonical.csv"))
write_csv(treat %>% transmute(raw_code, canonical_code, canonical_name), treat_xwalk_csv)
# treating include: bowel-surgery sites that are not dissolved
write_csv(treat_canon %>% filter(bowel_ca_surgery == 1, !any_dissolved) %>%
            transmute(canonical_code), treat_include_csv)

# 6. review summary ----------------------------------------------------------
cat(sprintf("diagnosing: %d codes -> %d sites (%d unresolved/excluded)\n",
            nrow(diag), nrow(diag_canon), sum(is.na(diag$canonical_code))))
excl_n <- sum(diag$n_total[is.na(diag$canonical_code)])
cat(sprintf("diagnosing volume unattributable to a 5-digit site: %d patients (%.1f%% of diagnosing volume)\n",
            excl_n, 100 * excl_n / sum(diag$n_total)))
cat("\nTrust-code resolution (calibrate site_min_vol / dominance_share on trust_colon_site.csv):\n")
print(trust_colon %>% count(classification, wt = coalesce(trust_code_vol, 0L), name = "trust_coded_patients") %>%
        left_join(count(trust_colon, classification, name = "n_trusts"), by = "classification"))
cat(sprintf("treating:   %d codes -> %d sites (%d not bowel-surgery)\n",
            nrow(treat), nrow(treat_canon), sum(treat_canon$bowel_ca_surgery == 0)))
cat("\nDiagnosing codes not placed at a 5-digit site (check for single-site trusts):\n")
print(diag %>% filter(is.na(canonical_code)) %>% arrange(desc(n_total)) %>%
        select(raw_code, code_type, n_total, exclude_reason), n = 40)
cat("\nTreating codes with volume but not flagged bowel-surgery:\n")
print(treat %>% filter(bowel_ca_surgery == 0, n_total >= 20) %>% arrange(desc(n_total)) %>%
        select(raw_code, canonical_name, n_total), n = 30)