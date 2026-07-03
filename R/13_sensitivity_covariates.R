# 13  sensitivity of the ranking to the adjustment set (quick, unshrunk)
# -----------------------------------------------------------------------------
# A fast robustness check on the sustained ranking. The primary model standardises
# for age, comorbidity, season and calendar year. Here we re-estimate the sustained
# hospital means two other ways and check the ranking is stable:
#   - age + comorbidity only : the primary set with season and calendar year removed
#   - unstandardised         : raw hospital mean wait, no case-mix adjustment
# We report rank agreement with the primary and whether the fast-tail candidate
# set (fastest quintile) is stable, since that set is what the positive-deviance
# investigation acts on. These are plain weighted means, no shrinkage; the shrunk,
# per-hospital version of the calendar comparison is the table in 10, and the fit
# and effective-sample-size justification of the set is in 15.

library(balancer)
library(dplyr)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

std_site <- function(cont, bin, d = df) {
  d  <- d %>% arrange(hosp)
  X  <- make_std_matrix(d, cont, bin)
  so <- standardize(X, rep(0, ncol(X)), d$hosp, lambda = lambda_main, exact_global = FALSE)
  w  <- extract_weights(so)
  tibble(hosp = d$hosp, wait = d$wait, w = w) %>%
    group_by(hosp) %>% summarise(est = weighted.mean(wait, w), .groups = "drop")
}

cv <- code_covariates(df)
df <- cv$data                     # augmented with the cci dummies std_site needs

# primary set uses time_bin (season + calendar year); the reduced set drops it
est <- std_site(cv$cont, c(cv$bin, time_bin)) %>% rename(primary = est) %>%
  inner_join(std_site(cv$cont, cv$bin) %>% rename(age_cci = est), by = "hosp") %>%
  inner_join(df %>% group_by(hosp) %>% summarise(unstandardised = mean(wait), .groups = "drop"),
             by = "hosp")

# rank agreement (1 = fastest) -----------------------------------------------
ranks <- est %>% mutate(across(c(primary, age_cci, unstandardised),
                               ~ rank(.x, ties.method = "average"), .names = "r_{.col}"))
rho <- cor(ranks %>% select(starts_with("r_")), method = "spearman")
cat("Spearman rank correlation with the primary model:\n")
print(round(rho, 3))

# stability of the fast-tail candidate set (fastest quintile) -----------------
q <- 0.20
tag <- function(x) x <= quantile(x, q)
cand <- est %>% transmute(hosp,
                          primary = tag(primary), age_cci = tag(age_cci),
                          unstandardised = tag(unstandardised))
overlap <- function(a, b) sum(a & b) / sum(a)
cat(sprintf("\nFastest-quintile candidates: primary %d hospitals.\n", sum(cand$primary)))
cat(sprintf("Retained under age + comorbidity only: %.0f%%\n", 100 * overlap(cand$primary, cand$age_cci)))
cat(sprintf("Retained unstandardised:               %.0f%%\n", 100 * overlap(cand$primary, cand$unstandardised)))

write.csv(est,   file.path(out_dir, "sensitivity_estimates.csv"), row.names = FALSE)
write.csv(round(rho, 3), file.path(out_dir, "sensitivity_rank_correlation.csv"))
cat("\nsensitivity outputs written.\n")