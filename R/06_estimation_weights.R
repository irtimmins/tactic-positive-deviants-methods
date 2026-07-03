# Directly standardised hospital waiting times via balancing weights.
# Mirrors the Estimation_Step / Post_Estimation steps of the template:
#  - build the standardised balance matrix and reweight each hospital to the
#    population covariate means using balancer::standardize()
#  - a pooled prognostic model gives residuals for the augmented (residual
#    balancing) estimate and the population-mean prediction
#  - per-hospital weighted and augmented means, effective n and pooled SEs
# Produces inputs for the shrinkage step for: the sustained estimand (whole
# window), the improvement estimand (second half minus first), and the two
# comorbidity strata.

library(balancer)
library(dplyr)
library(ggplot2)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# balance vs effective-sample-size trade-off across lambda --------------------
# prognostic weights for the imbalance metric come from the pooled model on the
# standardised covariates, so imbalance is summarised on the outcome scale.
balance_tradeoff <- function(d, cont, bin, grid = lambda_grid) {
  d  <- d %>% arrange(hosp)
  X  <- make_std_matrix(d, cont, bin)
  Z  <- d$hosp
  pm <- lm(as.formula(paste("y_std ~", paste(colnames(X), collapse = " + "))),
           data = data.frame(y_std = d$wait, X))
  beta <- coef(pm)[-1]
  beta[is.na(beta)] <- 0
  
  hosp_means <- rowsum(X, Z) / as.numeric(table(Z))    # hospital x covariate
  unw <- as.numeric(abs((hosp_means %*% beta)))
  
  res <- t(sapply(grid, function(l) {
    so <- standardize(X, rep(0, ncol(X)), Z, lambda = l, exact_global = FALSE)
    w  <- extract_weights(so)
    wm <- (t(so$weights) %*% X)                        # hospital x covariate
    wt <- as.numeric(abs(wm %*% beta))
    ne <- tapply(w, Z, function(x) sum(x)^2 / sum(x^2))
    c(lambda = l,
      bias_removed = 1 - mean(wt) / mean(unw),
      mean_eff_n   = mean(ne),
      mean_deff    = mean(ne / as.numeric(table(Z))))
  }))
  as.data.frame(res)
}

# main analysis: primary standardisation (age + cci + season + calendar year) -
# season and calendar year are part of the primary adjustment set; their
# inclusion, and the choice to keep age linear, are justified against outcome
# fit and effective sample size in script 15. Other patient factors (sex,
# ethnicity, deprivation, stage) remain excluded.
cv          <- code_covariates(df)
primary_bin <- c(cv$bin, season_terms, year_term)
trade       <- balance_tradeoff(cv$data, cv$cont, primary_bin)
cat("balance vs effective n by lambda:\n"); print(round(trade, 3))
write.csv(trade, file.path(out_dir, "lambda_tradeoff.csv"), row.names = FALSE)

# bias-variance trade-off curve (after Keele et al). Each point is one lambda:
# how much case-mix imbalance the weights remove (percent bias reduced) against
# the average effective sample size they keep. Small lambda sits top-left (most
# bias removed, fewest effective patients); larger lambda moves down and right.
# The elbow is the practical sweet spot; the working value lambda_main is marked.
trade_curve <- trade %>%
  mutate(pct_bias_removed = 100 * bias_removed,
         is_main = abs(lambda - lambda_main) < 1e-9)

p_trade <- ggplot(trade_curve, aes(mean_eff_n, pct_bias_removed)) +
  geom_path(colour = "grey60") +
  geom_point(aes(colour = is_main), size = 2) +
  geom_text(aes(label = lambda), vjust = -0.8, size = 2.8) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick"),
                      guide = "none") +
  labs(x = "Average effective sample size per hospital",
       y = "Case-mix bias removed (%)",
       title = "Bias reduction against effective sample size across lambda",
       subtitle = "each point is a lambda value; the working value is highlighted") +
  theme_bw()

ggsave(file.path(out_dir, "lambda_tradeoff.pdf"), p_trade, width = 7, height = 5)

fit_main <- run_standardise(patient_data          = cv$data,
                            continuous_covariates = cv$cont,
                            binary_covariates     = primary_bin,
                            lambda                = lambda_main)
site_sustained <- fit_main$site
saveRDS(fit_main,       file.path(out_dir, "fit_primary.rds"))
saveRDS(site_sustained, file.path(out_dir, "site_sustained.rds"))

# improvement estimand -------------------------------------------------------
# baseline period = first half of the window; later period = second half. Both
# are standardised to the SAME reference: the case-mix of the baseline period.
# The later-year indicator is constant within a half and is dropped by
# run_standardise, so the improvement estimand adjusts age + cci + season.
site_improve <- standardise_change(patient_data          = cv$data,
                                   continuous_covariates = cv$cont,
                                   binary_covariates     = primary_bin)
saveRDS(site_improve, file.path(out_dir, "site_improve.rds"))

# comorbidity strata ---------------------------------------------------------
# within a stratum comorbidity is near-constant, so adjust for age only and
# standardise to the stratum population. Season and calendar year are NOT added
# here, to preserve effective sample size in the smaller stratum samples; add
# c(season_terms, year_term) below if you prefer the strata to mirror the
# primary adjustment set exactly.
for (st in levels(df$cci_strata)) {
  cvs <- code_covariates(filter(df, cci_strata == st), cci = "none")
  fit_st <- run_standardise(patient_data          = cvs$data,
                            continuous_covariates = cvs$cont,
                            binary_covariates     = cvs$bin,
                            lambda                = lambda_main)
  saveRDS(fit_st$site,
          file.path(out_dir, sprintf("site_strata_%s.rds", gsub("[^0-9a-zA-Z]", "", st))))
}