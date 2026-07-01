# 12  sensitivity analyses on the sustained estimand
# -----------------------------------------------------------------------------
# The main model standardises for age and comorbidity only. Here we re-estimate
# the sustained hospital means under two pre-specified alternatives and check the
# hospital ranking is robust:
#   - time-adjusted : main case-mix plus season (quarter) and calendar year
#   - unstandardised: raw hospital mean wait, no case-mix adjustment
# We report rank agreement with the main model and whether the fast-tail
# candidate set (fastest quintile) is stable, since that set is what the
# positive-deviance investigation acts on.

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

est <- std_site(cv$cont, cv$bin) %>% rename(main = est) %>%
  inner_join(std_site(cv$cont, c(cv$bin, time_bin)) %>% rename(time_adj = est), by = "hosp") %>%
  inner_join(df %>% group_by(hosp) %>% summarise(unstandardised = mean(wait), .groups = "drop"),
             by = "hosp")

# rank agreement (1 = fastest) -----------------------------------------------
ranks <- est %>% mutate(across(c(main, time_adj, unstandardised),
                               ~ rank(.x, ties.method = "average"), .names = "r_{.col}"))
rho <- cor(ranks %>% select(starts_with("r_")), method = "spearman")
cat("Spearman rank correlation with the main model:\n")
print(round(rho, 3))

# stability of the fast-tail candidate set (fastest quintile) -----------------
q <- 0.20
tag <- function(x) x <= quantile(x, q)
cand <- est %>% transmute(hosp,
                          main = tag(main), time_adj = tag(time_adj),
                          unstandardised = tag(unstandardised))
overlap <- function(a, b) sum(a & b) / sum(a)
cat(sprintf("\nFastest-quintile candidates: main %d hospitals.\n", sum(cand$main)))
cat(sprintf("Retained under time-adjusted:  %.0f%%\n", 100 * overlap(cand$main, cand$time_adj)))
cat(sprintf("Retained under unstandardised: %.0f%%\n", 100 * overlap(cand$main, cand$unstandardised)))

write.csv(est,   file.path(out_dir, "sensitivity_estimates.csv"), row.names = FALSE)
write.csv(round(rho, 3), file.path(out_dir, "sensitivity_rank_correlation.csv"))
cat("\nsensitivity outputs written.\n")