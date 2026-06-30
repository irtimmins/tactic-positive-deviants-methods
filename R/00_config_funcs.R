# Shared settings and small helpers for the colon waiting-time
# positive-deviance analysis. Source this at the top of every other script.
#
# The analysis follows the balancer standardize() workflow (Ben-Michael / Keele):
# directly standardised hospital means via balancing weights, an augmented
# residual-balancing estimate, pooled SEs, then a normal-normal shrinkage model
# fitted to the hospital-level summaries, and posterior ranks from that fit.

# paths ---------------------------------------------------------------------
# edit these two for your machine
base_dir <- "D:/Projects/#2045_ICON_TACTIC/Project1_interim_bowel/tactic-bowel-quantifying-variation/Data/ICON"
in_rds   <- file.path(base_dir, "colon_cohort_cWT_2015_2022.rds")
out_dir  <- file.path("Output")
stan_dir <- file.path("R")                     # where the .stan files live
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# analysis window -----------------------------------------------------------
# restrict to the most recent two years of diagnoses. window_end defaults to the
# latest diagnosis date in the data; set window_months to change the length.
window_months <- 24
window_end    <- NA          # NA = use max(diagmdy); or set e.g. as.Date("2022-12-31")

# minimum diagnosing-hospital volume over the window to be profiled
min_volume <- 30

# outcome -------------------------------------------------------------------
outcome_var <- "wt_dx_to_dtt"   # days from diagnosis to decision-to-treat
max_wait    <- 180              # drop waits above this (implausible)
drop_zero_wait <- TRUE          # drop exact zero-day waits (see notes in script 01)

# balance sets --------------------------------------------------------------
# continuous covariates are z-scored; binary covariates are scaled by
# 1/sqrt(p(1-p)) with a floor at p = 0.05, as in the template.
cont_vars <- c("agediag", "cci_n_conditions")

# binary covariates are built as 0/1 dummies in script 01. the "primary" set is
# age + comorbidity + calendar year + season; "full" adds the patient mix.
bin_primary <- c("yr_late", "q2", "q3", "q4")
bin_full    <- c(bin_primary,
                 "male",
                 "eth_asian", "eth_black", "eth_mixed", "eth_other", "eth_unknown",
                 "imd_2", "imd_3", "imd_4", "imd_5",
                 "stage_2", "stage_3")

# shrinkage priors (day scale) ----------------------------------------------
# dp_normal.stan in the template is written for a proportion in [0,1]; for a
# waiting time in days we drop the [0,1] bounds and set the prior scale to the
# day scale. these feed dp_normal_cont.stan.
prior_mu_sd      <- 50    # sd of the normal prior on the grand mean
prior_tau_scale  <- 10    # half-cauchy scale for the between-hospital sd

# lambda grid for the balance / effective-sample-size trade-off
lambda_grid <- c(0, .01, .05, .1, .25, .5, 1, 1.5, 2, 2.5, 3)
lambda_main <- 0.05       # working value used for the headline results

# helpers -------------------------------------------------------------------

# z-score a continuous covariate
z_std <- function(x) (x - mean(x)) / sd(x)

# binary scaling with a floor on rare categories
bin_std <- function(x) {
  p <- mean(x)
  if (p < 0.05) p <- 0.05
  (x - mean(x)) / sqrt(p * (1 - p))
}

# build the standardised balance matrix from a data frame
make_std_matrix <- function(df, cont, bin) {
  Xc <- sapply(df[cont], z_std)
  Xb <- sapply(df[bin],  bin_std)
  cbind(Xc, Xb)
}

# reference moments (means, sds, proportions) of a target population, used when
# standardising one group to a different population than its own.
ref_moments <- function(df, cont, bin) {
  list(cont_mean = sapply(df[cont], mean),
       cont_sd   = sapply(df[cont], sd),
       bin_p     = sapply(df[bin],  mean))
}

# standardised balance matrix using externally supplied reference moments, so a
# zero target balances each group to that reference population rather than to
# its own means. continuous: (x - ref_mean) / ref_sd; binary: (x - ref_p) scaled
# by the floored reference proportion.
make_std_matrix_ref <- function(df, cont, bin, ref) {
  Xc <- mapply(function(v, m, s) (df[[v]] - m) / s,
               cont, ref$cont_mean[cont], ref$cont_sd[cont])
  ctr <- ref$bin_p[bin]
  scl <- sqrt(pmax(ctr, 0.05) * (1 - pmax(ctr, 0.05)))
  Xb <- mapply(function(v, c0, s) (df[[v]] - c0) / s, bin, ctr, scl)
  out <- cbind(Xc, Xb)
  colnames(out) <- c(cont, bin)
  out
}

# pull the per-patient weight out of a balancer standardize() result.
# weights come back as a patients-by-hospitals matrix with the weight sitting in
# the patient's own hospital column, so the row max recovers it.
extract_weights <- function(std_out) apply(std_out$weights, 1, max)

# weighted quantile (used for standardised medians). robust to small or
# degenerate cells: ignores zero-weight rows, returns the single value when a
# group has no spread, and collapses tied cumulative weights before interpolating.
w_quantile <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  if (length(unique(x)) == 1) return(x[1])
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  dup <- duplicated(cw, fromLast = TRUE)
  cw <- cw[!dup]; x <- x[!dup]
  if (length(x) < 2) return(x[length(x)])
  approx(cw, x, xout = p, rule = 2, ties = "ordered")$y
}

# per-hospital standardised summaries for a continuous outcome.
# df must contain: hosp (id), y (outcome), w (weight), resid (outcome-model
# residual), canonical (population-mean model prediction, constant for a pooled
# model). returns one row per hospital with the weighted and augmented estimates,
# effective n and pooled SEs, matching the template's hosp.data construction.
site_summary <- function(df) {
  library(dplyr)
  s <- df %>%
    group_by(hosp) %>%
    summarise(
      n          = n(),
      n_eff      = sum(w)^2 / sum(w^2),
      raw_mean   = mean(y),
      raw_median = median(y),
      stand      = weighted.mean(y, w),
      stand_med  = w_quantile(y, w, 0.5),
      stand_adj  = weighted.mean(resid, w) + mean(canonical),
      sd_w       = sqrt(sum(w^2 * (y - weighted.mean(y, w))^2) / sum(w^2)),
      sd_adj     = sqrt(sum(w^2 * resid^2) / sum(w^2)),
      .groups = "drop"
    )
  sd_pool     <- sqrt(weighted.mean(s$sd_w^2,   s$n_eff))
  sd_pool_adj <- sqrt(weighted.mean(s$sd_adj^2, s$n_eff))
  s %>% mutate(
    se          = sd_w   / sqrt(n_eff),
    se_pool     = sd_pool / sqrt(n_eff),
    se_adj      = sd_adj / sqrt(n_eff),
    se_adj_pool = sd_pool_adj / sqrt(n_eff),
    sd_pool     = sd_pool,
    sd_pool_adj = sd_pool_adj
  )
}

# posterior ranking metrics from a draws-by-hospital matrix of latent means.
# lower outcome is better (shorter wait), so rank 1 is the best performer.
rank_metrics <- function(draws, tops = c(.05, .10, .20, .50)) {
  J <- ncol(draws)
  R <- t(apply(draws, 1, rank, ties.method = "average"))   # draws x J
  out <- data.frame(
    exp_rank = colMeans(R),
    rank_lo  = apply(R, 2, quantile, 0.025),
    rank_hi  = apply(R, 2, quantile, 0.975)
  )
  for (p in tops) {
    thr <- ceiling(p * J)
    out[[paste0("p_top", p * 100)]] <- colMeans(R <= thr)
  }
  out
}