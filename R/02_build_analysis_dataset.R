# 02  build the analysis dataset (clinical eligibility)
# -----------------------------------------------------------------------------
# Draws the analysis sample through the clinical funnel: all colon (C18) adults
# diagnosed in the window, then stage 4 / non-elective / DCO / diagnosed-at-death
# / later primaries removed to the eligible cohort, then the surgical and
# decision-to-treat requirements. Hospital codes are left raw here; provider
# identity and volume are handled in 04, after the crosswalks are built in 03.
# Writes the clinical analysis set and the top of the flowchart.

library(dplyr)
library(lubridate)

source("R/01_config.R")

reg    <- readRDS(registry_rds)   # all C18 adults 2015+, inclusion flags
cohort <- readRDS(in_rds)         # eligible cohort + surgery + CWT + covariates

flow <- tibble(box = integer(), kind = character(), label = character(),
               n = integer(), n_hosp = integer())
.prev <- NA_integer_
fc_box <- function(box, label, n, n_hosp = NA_integer_) {
  flow  <<- bind_rows(flow, tibble(box = box, kind = "box", label = label,
                                   n = as.integer(n), n_hosp = as.integer(n_hosp)))
  .prev <<- as.integer(n)
}
fc_excl <- function(box, reason, n_after) {
  flow  <<- bind_rows(flow, tibble(box = box, kind = "exclusion", label = reason,
                                   n = .prev - as.integer(n_after), n_hosp = NA_integer_))
  .prev <<- as.integer(n_after)
}

# analysis window ------------------------------------------------------------
end_date   <- if (is.na(window_end)) max(reg$diagmdy, na.rm = TRUE) else as.Date(window_end)
start_date <- end_date %m-% months(window_months) %m+% days(1)
cat(sprintf("window: %s to %s\n", start_date, end_date))

w  <- reg %>% filter(diagmdy >= start_date, diagmdy <= end_date)
np <- function(mask) n_distinct(w$pseudo_patientid[mask])

# box 1: all colon adults in the window -------------------------------------
fc_box(1, sprintf("All adult patients with colon cancer (C18) diagnosed between %s and %s",
                  format(start_date, "%b %Y"), format(end_date, "%b %Y")),
       n_distinct(w$pseudo_patientid))

# box 2: stage 1-3 -----------------------------------------------------------
m <- w$incl_stage13
fc_excl(2, "Stage 4 / unknown / unstaged", np(m))
fc_box(2, "Patients with Stage 1-3 colon cancer", np(m))

# box 3: eligible incident colon cancer --------------------------------------
m <- m & w$incl_elective_route; fc_excl(3, "Non-elective / unknown route to diagnosis", np(m))
m <- m & w$incl_not_dco;        fc_excl(3, "Death-certificate-only diagnosis", np(m))
m <- m & w$incl_not_diag_death; fc_excl(3, "Diagnosed at death", np(m))
m <- m & w$is_first_primary;    fc_excl(3, "Not the first primary tumour", np(m))
fc_box(3, "Patients with non-emergency routes to diagnosis", np(m))

elig_ids <- unique(w$pseudo_patientid[m])

# join eligible in-window patients to the built cohort for surgery / CWT / covars
df <- cohort %>% filter(pseudo_patientid %in% elig_ids)
if (nrow(df) != length(elig_ids))
  cat(sprintf("note: %d eligible patients, %d matched in the built cohort\n",
              length(elig_ids), nrow(df)))

# box 4: underwent an elective colon resection -------------------------------
df <- df %>% filter(had_surgery, !emergency)
fc_excl(4, "No elective colon resection", nrow(df))
fc_box(4, "Underwent elective colon resection", nrow(df), n_distinct(df$diag_hosp))

# box 5: valid, plausible waiting time ---------------------------------------
df <- df %>% filter(dtt_valid)
fc_excl(5, "No valid decision-to-treat recorded", nrow(df))
df$wait <- df[[outcome_var]]
df <- df %>% filter(wait <= max_wait)
fc_excl(5, sprintf("Wait > %d days", max_wait), nrow(df))
if (drop_zero_wait) {
  df <- df %>% filter(wait > 0)
  fc_excl(5, "Wait <= 0 days", nrow(df))
}
fc_box(5, "Valid waiting time (DTT recorded, plausible wait)", nrow(df), n_distinct(df$diag_hosp))

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

# NHS England region for representativeness checks. Registry geography fields
# follow the <geo>_2024_name pattern (as canalliance does); take the first region
# field present, else leave region blank (Table 1 then omits the region block).
region_src <- intersect(c("nhser_2024_name", "nhs_england_region_2024_name",
                          "region_2024_name", "nhser_name"), names(df))
df$region <- if (length(region_src)) as.character(df[[region_src[1]]]) else NA_character_
if (!length(region_src))
  cat("no NHS region field found in the registry - set region_src in 02 if needed\n")

cat("\nflowchart (boxes and exclusions):\n"); print(as.data.frame(flow))
cat(sprintf("\nclinically eligible: %d patients, %d raw diagnosing codes\n",
            nrow(df), n_distinct(df$diag_hosp)))

saveRDS(df, clinical_rds)
write.csv(as.data.frame(flow), flow_clinical_csv, row.names = FALSE)