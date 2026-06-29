# Build the analysis dataset: restrict to the latest two years of diagnoses,
# pick the diagnosing hospital as the unit, define the outcome and the balance
# covariates, drop low-volume hospitals, and split the window into two halves
# for the improvement dimension.

library(dplyr)
library(lubridate)

source("00_config_funcs.R")

df <- readRDS(in_rds)

# attrition tracking
att <- tibble(step = character(), n = integer(), n_hosp = integer())
log_step <- function(d, label, hosp = "diag_hosp") {
  att <<- bind_rows(att, tibble(step = label, n = nrow(d),
                                n_hosp = n_distinct(d[[hosp]])))
  d
}

df <- df %>% log_step("start")

# valid diagnosing hospital (5-char ODS code), as in the prep scripts
df <- df %>%
  mutate(diag_hosp = trimws(as.character(diag_hosp))) %>%
  filter(grepl("^R[A-Z0-9]{4}$", diag_hosp)) %>%
  log_step("valid hospital code")

# analysis window
end_date   <- if (is.na(window_end)) max(df$diagmdy, na.rm = TRUE) else as.Date(window_end)
start_date <- end_date %m-% months(window_months) %m+% days(1)
cat(sprintf("window: %s to %s\n", start_date, end_date))

df <- df %>%
  filter(diagmdy >= start_date, diagmdy <= end_date) %>%
  log_step("within window")

# outcome and implausible-wait exclusions.
# the dataset already enforces diag <= dtt (no negatives) and the max observed
# wait is under 180, so in practice only the zero-day waits are affected here.
# zero-day diagnosis-to-DTT can be genuine fast-track activity, so dropping it
# removes some of the fastest pathways - left as a toggle (drop_zero_wait).
df$wait <- df[[outcome_var]]
df <- df %>% filter(wait <= max_wait) %>% log_step("wait <= max")
if (drop_zero_wait) df <- df %>% filter(wait > 0) %>% log_step("wait > 0")

# covariates ----------------------------------------------------------------
# sex is the raw NCRAS coding (1/2). assumed 1 = male, 2 = female - confirm.
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
      ethnicity_group_broad == "White"     ~ "White",
      ethnicity_group_broad == "Asian"     ~ "Asian",
      ethnicity_group_broad == "Black"     ~ "Black",
      ethnicity_group_broad == "Mixed"     ~ "Mixed",
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
    
    # calendar terms: late-half-of-window indicator and season (quarter)
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

# any imd missing -> set dummies to 0 and flag (keeps rows in the primary model,
# which does not use imd; the full model can be re-run on complete imd if wanted)
df <- df %>% mutate(across(c(imd_2, imd_3, imd_4, imd_5),
                           ~ ifelse(is.na(.), 0L, .)))

# hospital index for balancer (numeric group id), keep the code as a label
df <- df %>%
  mutate(hospid = as.integer(factor(diag_hosp)),
         hosp   = hospid)

# drop low-volume hospitals
vol <- df %>% count(diag_hosp, name = "volume")
df  <- df %>% left_join(vol, by = "diag_hosp") %>%
  filter(volume > min_volume) %>%
  log_step("volume > min")

# re-index hospitals after the volume filter so ids are contiguous
df <- df %>% mutate(hospid = as.integer(factor(diag_hosp)), hosp = hospid)

hosp_lookup <- df %>%
  distinct(hospid, diag_hosp, diag_hosp_name, canalliance, volume) %>%
  arrange(hospid)

cat("\nattrition:\n"); print(as.data.frame(att))
cat(sprintf("\nfinal: %d patients, %d hospitals\n",
            nrow(df), n_distinct(df$diag_hosp)))
cat(sprintf("median wait %.0f days (IQR %.0f-%.0f)\n",
            median(df$wait), quantile(df$wait, .25), quantile(df$wait, .75)))

saveRDS(df,          file.path(out_dir, "analysis_data.rds"))
saveRDS(hosp_lookup, file.path(out_dir, "hosp_lookup.rds"))
write.csv(as.data.frame(att), file.path(out_dir, "attrition.csv"), row.names = FALSE)