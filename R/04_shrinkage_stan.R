# Fit the normal-normal shrinkage model to the hospital-level summaries, as in
# Shrinkage_Fits. The augmented weighted estimate and its pooled SE are the
# inputs, matching the template (stand.comp.adj / SE.adj.hat.pool). One fit for
# the sustained estimand, one for improvement, one per comorbidity stratum.

library(rstan)
library(dplyr)

source("R/00_config_funcs.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_file <- file.path(stan_dir, "dp_normal_cont.stan")

fit_shrink <- function(y, se, mu_mean = mean(y),
                       mu_sd = prior_mu_sd, tau_scale = prior_tau_scale,
                       seed = 8675309, refresh = 2000, cores = 1,
                       adapt_delta = 0.95) {
  dat <- list(J = length(y), y_site_obs = y, sigma_site_obs = se,
              prior_mu_mean = mu_mean, prior_mu_sd = mu_sd,
              prior_tau_scale = tau_scale)
  # cores = 1 keeps sampling in the main process so the per-iteration progress
  # prints to the console; on Windows parallel workers hide it. these models are
  # tiny so sequential chains cost nothing. raise cores for speed if you do not
  # need the live progress. adapt_delta above the 0.8 default takes smaller steps
  # to suppress divergences in low-signal fits (e.g. the high-comorbidity stratum).
  stan(stan_file, data = dat, seed = seed,
       chains = 4, iter = 4000, warmup = 2000, refresh = refresh, cores = cores,
       control = list(adapt_delta = adapt_delta, max_treedepth = 12))
}

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