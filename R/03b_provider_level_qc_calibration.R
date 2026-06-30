# qc  calibrate the single-colon-site detection (run after 03, not in the pipeline)
# -----------------------------------------------------------------------------
# Reads what 03 writes and helps you set two numbers and find code aliases:
#   site_min_vol     - colon volume that makes a trust's site "real" not noise
#   dominance_share  - share one site must hold for the trust code to map to it
# It sweeps both, shows the recovered-vs-excluded trade-off, lists the trusts
# sitting near the dominance line (where the call could go either way), and flags
# trusts whose top two sites look like the same hospital under two codes - the
# candidates for code_aliases back in 03.

library(tidyverse)

source("R/01_config.R")

tsv  <- read_csv(file.path(qc_dir, "trust_site_volumes.csv"), show_col_types = FALSE)
nons <- read_csv(file.path(qc_dir, "nonstandard_codes.csv"),  show_col_types = FALSE)

# re-derive the classification for any pair of thresholds, from the raw per-site
# volumes - so the sweep matches exactly what 03 would do
classify <- function(d, min_vol, dom_share) {
  d %>% arrange(trust, desc(n_total)) %>%
    group_by(trust) %>%
    summarise(n_sites_real = sum(n_total >= min_vol),
              trust_site_vol = sum(n_total),
              top_site = first(site), top_vol = first(n_total),
              top_share = first(n_total) / sum(n_total), .groups = "drop") %>%
    mutate(resolved = top_vol >= min_vol & top_share >= dom_share)
}

# 1. sweep ---------------------------------------------------------------------
grid <- expand_grid(site_min_vol    = c(10, 15, 20, 30),
                    dominance_share = c(0.60, 0.70, 0.80, 0.90))
sweep <- grid %>%
  mutate(res = map2(site_min_vol, dominance_share, ~ classify(tsv, .x, .y))) %>%
  mutate(trusts_resolved = map_int(res, ~ sum(.x$resolved)),
         trusts_total    = map_int(res, ~ nrow(.x))) %>%
  select(-res)
cat("Trusts resolved (mapped to one colon site) by threshold:\n")
print(sweep %>% pivot_wider(names_from = dominance_share, values_from = trusts_resolved,
                            id_cols = site_min_vol, names_prefix = "dom_"))

# 2. borderline trusts at the current settings --------------------------------
# trusts whose dominant site sits just below the line - lowering dominance_share
# would recover them, so check they really are one colon site
cur <- classify(tsv, site_min_vol, dominance_share)
borderline <- cur %>%
  filter(!resolved, top_share >= dominance_share - 0.20) %>%
  arrange(desc(top_vol)) %>%
  left_join(distinct(tsv, site, site_name), by = c("top_site" = "site"))
cat(sprintf("\nBorderline trusts (top share within 0.20 below dominance = %.2f):\n", dominance_share))
print(borderline %>% transmute(trust, top_site, top_name = site_name,
                               top_share = round(top_share, 2), top_vol, trust_site_vol), n = 30)

# 3. same-hospital alias candidates -------------------------------------------
# within a trust, two sites whose names look alike are probably one hospital under
# two codes; fold the smaller into the larger via code_aliases in 03
norm <- function(x) str_squish(str_to_lower(str_replace_all(coalesce(x, ""), "\\(.*\\)|hospital|nhs|trust|general|elective surgical hub", "")))
pairs <- tsv %>% group_by(trust) %>% filter(n() >= 2) %>%
  arrange(desc(n_total), .by_group = TRUE) %>%
  summarise(a_site = site[1], a_name = site_name[1], a_n = n_total[1],
            b_site = site[2], b_name = site_name[2], b_n = n_total[2], .groups = "drop") %>%
  mutate(name_overlap = norm(a_name) == norm(b_name) & norm(a_name) != "") %>%
  filter(name_overlap)
cat("\nPossible same-hospital pairs (consider aliasing b_site -> a_site):\n")
print(pairs %>% transmute(trust, a_site, b_site, a_name, b_name, a_n, b_n), n = 30)

# 4. non-standard codes still unaliased ---------------------------------------
cat("\nNon-standard codes (hubs / recoded sites) to map in code_aliases:\n")
print(nons %>% arrange(desc(diag_vol + treat_vol)), n = 40)

# 5. postcode-based reassignment ----------------------------------------------
# the strongest objective evidence that a hub / recoded code is really an existing
# hospital is a shared ODS postcode (same building). Match each non-standard code
# to the standard 5-char site at the same postcode: equal postcode plus agreeing
# name is a safe alias. No postcode match -> cannot reassign reliably -> exclude.
ods <- read_csv(file.path(qc_dir, "ods_cache.csv"), col_types = cols(.default = "c"))
std <- ods %>% filter(grepl("^R[A-Z0-9]{4}$", code), !is.na(postcode)) %>%
  transmute(std_code = code, postcode, std_name = ods_name)
reassign <- nons %>%
  left_join(ods %>% select(code, postcode, ns_name = ods_name), by = "code") %>%
  left_join(std, by = "postcode", relationship = "many-to-many") %>%
  filter(!is.na(std_code), std_code != code) %>%
  arrange(desc(diag_vol + treat_vol)) %>%
  transmute(alt_code = code, alt_name = ns_name, postcode, std_code, std_name, diag_vol, treat_vol)
cat("\nNon-standard codes sharing a postcode with a standard site (evidenced alias candidates):\n")
print(reassign, n = 40)
unmatched <- setdiff(nons$code, reassign$alt_code)
cat(sprintf("\nNon-standard codes with no postcode match (cannot reassign, stay excluded): %s\n",
            if (length(unmatched)) paste(unmatched, collapse = ", ") else "none"))

# 6. any two codes at one postcode (same building) ----------------------------
# generalises section 5 to standard codes too: within a trust, two ordinary site
# codes at the same postcode are one hospital fragmented across codes, which can
# make a single-site trust look multi-site. These are merge candidates that may
# let a borderline trust resolve cleanly rather than be excluded on assumption.
allsites <- bind_rows(tsv %>% transmute(code = site, vol = n_total)) %>%
  group_by(code) %>% summarise(vol = sum(vol), .groups = "drop") %>%
  left_join(ods %>% select(code, postcode, parent_code, name = ods_name), by = "code") %>%
  filter(!is.na(postcode))
samepc <- allsites %>% group_by(parent_code, postcode) %>%
  filter(n_distinct(code) >= 2) %>%
  arrange(parent_code, postcode, desc(vol)) %>% ungroup()
cat("\nMultiple codes at one postcode within a trust (same hospital, merge candidates):\n")
print(samepc %>% transmute(trust = parent_code, postcode, code, name, vol), n = 50)