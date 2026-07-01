# explore_covariate_coding  choose how to code age and cci (run once, not in the
# pipeline). Two questions: (1) is wait non-linear in age / cci, so that mean-
# balancing a continuous term leaves residual imbalance, and (2) do hospital ranks
# actually move between codings. If ranks are stable, prefer the simpler, higher
# effective-sample-size coding and report robustness.

library(balancer)
library(dplyr)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# 1. functional form: mean wait across the natural categories -----------------
# if the step changes between adjacent levels are uneven, a linear (mean-balanced)
# term will not capture the relationship and a categorical coding is safer.
cat("Mean wait by age band:\n")
df %>% mutate(band = cut(agediag, c(-Inf, 50, 60, 70, 80, Inf),
                         labels = c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  group_by(band) %>% summarise(n = n(), mean_wait = mean(wait), .groups = "drop") %>%
  mutate(step = mean_wait - lag(mean_wait)) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nMean wait by number of Charlson conditions:\n")
df %>% mutate(cci = pmin(cci_n_conditions, 5)) %>%
  group_by(cci) %>% summarise(n = n(), mean_wait = mean(wait), .groups = "drop") %>%
  mutate(step = mean_wait - lag(mean_wait)) %>% as.data.frame() %>% print(row.names = FALSE)

# a linearity check: does adding squared / categorical terms improve fit over
# the linear term (lower AIC = better)?
cat("\nLinearity check (AIC; lower is better):\n")
aic_tbl <- tibble(
  model = c("age linear", "age + age^2", "age banded",
            "cci linear", "cci + cci^2", "cci categorical (0/1/2/3+)"),
  AIC = c(
    AIC(lm(wait ~ agediag, df)),
    AIC(lm(wait ~ poly(agediag, 2), df)),
    AIC(lm(wait ~ cut(agediag, c(-Inf,50,60,70,80,Inf)), df)),
    AIC(lm(wait ~ cci_n_conditions, df)),
    AIC(lm(wait ~ poly(cci_n_conditions, 2), df)),
    AIC(lm(wait ~ factor(pmin(cci_n_conditions, 3)), df))))
print(as.data.frame(aic_tbl), row.names = FALSE)

# 2. does the coding move the hospital estimates / ranks? ---------------------
# standardise under each coding and compare the standardised means (rank Spearman
# against the continuous/continuous reference) and the mean effective sample size.
site_est <- function(age, cci) {
  cv <- code_covariates(df, age = age, cci = cci)
  d  <- cv$data %>% arrange(hosp)
  X  <- make_std_matrix(d, cv$cont, cv$bin)
  so <- standardize(X, rep(0, ncol(X)), d$hosp, lambda = lambda_main, exact_global = FALSE)
  w  <- extract_weights(so)
  tibble(hosp = d$hosp, wait = d$wait, w = w) %>%
    group_by(hosp) %>%
    summarise(stand = weighted.mean(wait, w),
              eff_n = sum(w)^2 / sum(w^2), .groups = "drop")
}

specs <- expand.grid(age = c("cont", "band"),
                     cci = c("cont", "0_1_2plus", "0_1_2_3plus"),
                     stringsAsFactors = FALSE)
ref <- site_est("cont", "cont")

cat("\nCoding comparison (rank agreement vs continuous/continuous, and ESS):\n")
res <- lapply(seq_len(nrow(specs)), function(i) {
  s <- site_est(specs$age[i], specs$cci[i])
  m <- inner_join(ref, s, by = "hosp", suffix = c("_ref", ""))
  tibble(age = specs$age[i], cci = specs$cci[i],
         rank_rho = cor(rank(m$stand_ref), rank(m$stand), method = "spearman"),
         mean_eff_n = mean(s$eff_n),
         median_shift = median(abs(m$stand - m$stand_ref)))
}) %>% bind_rows()
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

cat("\nRead: pick the simplest coding whose rank_rho vs the alternatives is high\n",
    "(ranks robust) and whose mean_eff_n is not badly eroded. A low rho or a large\n",
    "median_shift means the coding choice matters and should follow the AIC / step\n",
    "evidence above.\n")