# 02  build the analysis dataset
# -----------------------------------------------------------------------------
# Draws the analysis sample and a flowchart that starts from all colon (C18)
# adults diagnosed in the window, then excludes stage 4 / non-elective / DCO /
# diagnosed-at-death / later primaries to reach the eligible cohort, then joins
# the built cohort for the surgical and decision-to-treat steps and the volume
# threshold. The top of the funnel comes from the registry file (one row per
# tumour, with inclusion flags); the bottom comes from the built cohort.

library(dplyr)
library(lubridate)

source("R/01_config.R")

reg    <- readRDS(registry_rds)   # all C18 adults 2015+, inclusion flags
cohort <- readRDS(in_rds)         # eligible cohort + surgery + CWT + covariates

flow <- tibble(step = character(), n = integer(), n_hosp = integer())
add  <- function(label, n, n_hosp = NA_integer_)
  flow <<- bind_rows(flow, tibble(step = label, n = n, n_hosp = n_hosp))

# analysis window ------------------------------------------------------------
end_date   <- if (is.na(window_end)) max(reg$diagmdy, na.rm = TRUE) else as.Date(window_end)
start_date <- end_date %m-% months(window_months) %m+% days(1)
cat(sprintf("window: %s to %s\n", start_date, end_date))

w  <- reg %>% filter(diagmdy >= start_date, diagmdy <= end_date)
np <- function(mask) n_distinct(w$pseudo_patientid[mask])   # distinct patients

# registry funnel: all C18 adults in window -> eligible cohort ---------------
add(sprintf("All colon (C18) adults diagnosed in window (%s to %s)", start_date, end_date),
    n_distinct(w$pseudo_patientid))
m <- w$incl_stage13;            add("Stage 1-3 (excl stage 4 / X / U)", np(m))
m <- m & w$incl_elective_route; add("Elective / known route to diagnosis", np(m))
m <- m & w$incl_not_dco;        add("Excl death-certificate-only", np(m))
m <- m & w$incl_not_diag_death; add("Excl diagnosed at death", np(m))
m <- m & w$is_first_primary;    add("First primary (eligible cohort)", np(m))

elig_ids <- unique(w$pseudo_patientid[m])

# join the eligible in-window patients to the built cohort for surgery, CWT and
# covariates (one row per patient). Sizes should match the eligible count above.
df <- cohort %>% filter(pseudo_patientid %in% elig_ids)
if (nrow(df) != length(elig_ids))
  cat(sprintf("note: %d eligible patients, %d matched in the built cohort\n",
              length(elig_ids), nrow(df)))

# treatment requirement ------------------------------------------------------
df <- df %>% filter(had_surgery)
add("Underwent a colon resection", nrow(df), n_distinct(df$diag_hosp))
df <- df %>% filter(!emergency)
add("Elective resection (non-emergency admission)", nrow(df), n_distinct(df$diag_hosp))
df <- df %>% filter(dtt_valid)
add("Valid decision-to-treat (clock node)", nrow(df), n_distinct(df$diag_hosp))

# outcome plausibility -------------------------------------------------------
# valid DTT already guarantees a non-negative wait; this bounds the upper tail
# and (optionally) removes exact zero-day waits, which can be genuine same-day
# fast-track activity - left as the drop_zero_wait toggle.
df$wait <- df[[outcome_var]]
df <- df %>% filter(wait <= max_wait)
add(sprintf("Wait <= %d days", max_wait), nrow(df), n_distinct(df$diag_hosp))
if (drop_zero_wait) {
  df <- df %>% filter(wait > 0)
  add("Wait > 0 days", nrow(df), n_distinct(df$diag_hosp))
}

# hospital identity ----------------------------------------------------------
# resolve recorded codes to canonical hospitals (identity unless a crosswalk is
# supplied), then optionally keep only a curated list of hospitals.
xwalk <- load_provider_xwalk()
df <- df %>% mutate(diag_hosp = trimws(as.character(diag_hosp)),
                    diag_hosp_canon = canonicalise_hosp(diag_hosp, xwalk))
if (file.exists(provider_include_csv)) {
  inc <- read.csv(provider_include_csv, colClasses = "character")$canonical_code
  df  <- df %>% filter(diag_hosp_canon %in% trimws(inc))
  add("Included hospitals (curated list)", nrow(df), n_distinct(df$diag_hosp_canon))
}

# covariates -----------------------------------------------------------------
# sex is the raw NCRAS coding (1/2); assumed 1 = male, 2 = female - confirm.
df <- df %>%
  mutate(
    age   = agediag,
    male  = as.integer(sex == 1),
    sexf  = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    
    cci_n_conditions = as.numeric(cci_n_conditions),
    cci_strata = factor(ifelse(cci_n_conditions >= 2, "2+", "0-1"),
                        levels = c("0-1", "2+")),
    
    stage = factor(stage, levels = c("1", "2", "3")),
    stage_2 = as.integer(stage == "2"),
    stage_3 = as.integer(stage == "3"),
    
    eth = case_when(
      ethnicity_group_broad == "White" ~ "White",
      ethnicity_group_broad == "Asian" ~ "Asian",
      ethnicity_group_broad == "Black" ~ "Black",
      ethnicity_group_broad == "Mixed" ~ "Mixed",
      ethnicity_group_broad %in% c("Chinese", "Other") ~ "Other",
      TRUE ~ "Unknown"
    ),
    eth = factor(eth, levels = c("White","Asian","Black","Mixed","Other","Unknown")),
    eth_asian   = as.integer(eth == "Asian"),
    eth_black   = as.integer(eth == "Black"),
    eth_mixed   = as.integer(eth == "Mixed"),
    eth_other   = as.integer(eth == "Other"),
    eth_unknown = as.integer(eth == "Unknown"),
    
    imd = factor(NHSE_reversed_imd_quintile_lsoas),
    imd_q = case_when(
      grepl("^1", NHSE_reversed_imd_quintile_lsoas) ~ 1L,
      grepl("^2", NHSE_reversed_imd_quintile_lsoas) ~ 2L,
      grepl("^3", NHSE_reversed_imd_quintile_lsoas) ~ 3L,
      grepl("^4", NHSE_reversed_imd_quintile_lsoas) ~ 4L,
      grepl("^5", NHSE_reversed_imd_quintile_lsoas) ~ 5L,
      TRUE ~ NA_integer_
    ),
    imd_2 = as.integer(imd_q == 2),
    imd_3 = as.integer(imd_q == 3),
    imd_4 = as.integer(imd_q == 4),
    imd_5 = as.integer(imd_q == 5),
    
    route = factor(route_combined),
    
    mid_date = start_date %m+% months(window_months / 2),
    period   = factor(ifelse(diagmdy < mid_date, "first", "second"),
                      levels = c("first", "second")),
    yr_late  = as.integer(period == "second"),
    qtr      = quarter(diagmdy),
    q2 = as.integer(qtr == 2),
    q3 = as.integer(qtr == 3),
    q4 = as.integer(qtr == 4),
    
    canalliance = canalliance_2024_name
  )

df <- df %>% mutate(across(c(imd_2, imd_3, imd_4, imd_5), ~ ifelse(is.na(.), 0L, .)))

# volume threshold on the analysis unit (canonical hospital) -----------------
vol <- df %>% count(diag_hosp_canon, name = "volume")
df  <- df %>% left_join(vol, by = "diag_hosp_canon") %>%
  filter(volume > min_volume)
add(sprintf("Hospital volume > %d in window", min_volume),
    nrow(df), n_distinct(df$diag_hosp_canon))

# numeric hospital id for the estimation step (contiguous after the filters)
df <- df %>% mutate(hospid = as.integer(factor(diag_hosp_canon)), hosp = hospid)

hosp_lookup <- df %>%
  distinct(hospid, diag_hosp_canon, diag_hosp, diag_hosp_name, canalliance, volume) %>%
  arrange(hospid)

cat("\nflowchart:\n"); print(as.data.frame(flow))
cat(sprintf("\nanalysis sample: %d patients, %d hospitals\n",
            nrow(df), n_distinct(df$diag_hosp_canon)))
cat(sprintf("median wait %.0f days (IQR %.0f-%.0f)\n",
            median(df$wait), quantile(df$wait, .25), quantile(df$wait, .75)))

saveRDS(df,          file.path(out_dir, "analysis_data.rds"))
saveRDS(hosp_lookup, file.path(out_dir, "hosp_lookup.rds"))
write.csv(as.data.frame(flow), file.path(out_dir, "analysis_flowchart.csv"), row.names = FALSE)