# Fit the normal-normal shrinkage model to the hospital-level summaries, as in
# Shrinkage_Fits. The augmented weighted estimate and its pooled SE are the
# inputs, matching the template (stand.comp.adj / SE.adj.hat.pool). One fit for
# the sustained estimand, one for improvement, one per comorbidity stratum. The
# shrinkage routine itself (fit_shrink) is shared, in 01_config.R.

library(rstan)
library(dplyr)

source("R/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# sustained ------------------------------------------------------------------
site <- readRDS(file.path(out_dir, "site_sustained.rds")) %>% arrange(hosp)
fit_sustained <- fit_shrink(site$stand_adj, site$se_adj_pool)
print(fit_sustained, pars = c("mu_true", "sigma_true"))
saveRDS(list(site = site, fit = fit_sustained),
        file.path(out_dir, "stan_sustained.rds"))

# improvement ----------------------------------------------------------------
imp <- readRDS(file.path(out_dir, "site_improve.rds")) %>% arrange(hosp)
fit_improve <- fit_shrink(imp$delta, imp$se_delta, mu_mean = 0)
print(fit_improve, pars = c("mu_true", "sigma_true"))
saveRDS(list(site = imp, fit = fit_improve),
        file.path(out_dir, "stan_improve.rds"))

# strata ---------------------------------------------------------------------
for (tag in c("01", "2")) {
  f <- file.path(out_dir, sprintf("site_strata_%s.rds", tag))
  if (!file.exists(f)) next
  s <- readRDS(f) %>% arrange(hosp)
  fit_s <- fit_shrink(s$stand_adj, s$se_adj_pool)
  saveRDS(list(site = s, fit = fit_s),
          file.path(out_dir, sprintf("stan_strata_%s.rds", tag)))
}

cat("shrinkage fits saved.\n")