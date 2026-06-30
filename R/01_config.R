# 01  configuration and shared helpers
# -----------------------------------------------------------------------------
# Sourced at the top of every analysis script. Part A holds the settings you may
# want to change; Part B holds helper functions you should not need to edit.
#
# What the analysis does, in plain terms:
#   For each diagnosing hospital we estimate the average time from diagnosis to
#   the decision to treat, adjusted for differences in patient mix (case-mix
#   standardisation), then stabilise small-hospital estimates with a Bayesian
#   shrinkage step before ranking hospitals and flagging consistently fast
#   ("positive deviant") providers. Case-mix adjustment is done two ways for
#   comparison: by reweighting each hospital's patients to the overall patient
#   mix (balancing weights), and by regression standardisation.

# ============================================================================
# PART A  -  settings
# ============================================================================

# paths ----------------------------------------------------------------------
# edit base_dir for your machine; it points at the folder holding the built
# cohort from the data-build repo.
base_dir <- "D:/Projects/#2045_ICON_TACTIC/Project1_interim_bowel/tactic-bowel-quantifying-variation/Data/ICON"
in_rds   <- file.path(base_dir, "colon_cohort_cwt_2015_2022.rds")  # eligible cohort + surgery + CWT
registry_rds <- file.path(base_dir, "colon_registry_2015_2022.rds")  # all C18 adults + inclusion flags
out_dir  <- "Output"      # all analysis outputs are written here
clinical_rds <- file.path(out_dir, "analysis_clinical.rds")  # 02 -> 04 hand-off
flow_clinical_csv <- file.path(out_dir, "flow_clinical.csv")
stan_dir <- "R"           # folder holding the Stan model file
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# analysis window ------------------------------------------------------------
# the analysis uses the most recent stretch of diagnoses. window_end defaults to
# the latest diagnosis date in the cohort; set it explicitly to fix the window.
window_months <- 24
window_end    <- NA       # NA = latest diagnosis date; or e.g. as.Date("2022-12-31")

# inclusion volume: a site must have at least min_per_year patients in EVERY
# calendar year of the window. Applied in 04 to both the diagnosing unit and the
# treating site, and to each half-window in the improvement estimand (06).
min_per_year <- 10

# single-colon-site detection in 03 (calibrate with qc_calibrate_colon_site.R).
# A trust code is mapped to one 5-digit site only when colon activity clearly
# concentrates there: a site needs site_min_vol patients to count as real, and one
# site must hold dominance_share of the trust's site-coded colon volume.
site_min_vol    <- 20
dominance_share <- 0.95

# outcome --------------------------------------------------------------------
outcome_var    <- "wt_dx_to_dtt"   # days from diagnosis to decision-to-treat
max_wait       <- 180              # exclude waits longer than this as implausible
drop_zero_wait <- TRUE             # exclude exact zero-day waits (see note in 02)

# provider QC (built by 03, applied by 04) ----------------------------------
# 03 reconciles the cohort codes, the curated site Excel and the ODS API into
# canonical 5-digit diagnosing and treating sites. 04 reads these to canonicalise
# and to keep only valid providers. Curate the crosswalks once, then rerun 04.
qc_dir          <- file.path(out_dir, "provider_qc")
site_xlsx       <- file.path("Data/Site_level", "NHSHospitals_services_5.3.26_with_colours.xlsx")
diag_xwalk_csv  <- file.path(qc_dir, "diag_crosswalk.csv")    # raw_code, canonical_code, canonical_name
treat_xwalk_csv <- file.path(qc_dir, "treat_crosswalk.csv")
diag_include_csv  <- file.path(qc_dir, "diag_include.csv")    # canonical_code (valid diagnosing sites)
treat_include_csv <- file.path(qc_dir, "treat_include.csv")  # canonical_code (bowel-surgery sites)
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

# covariates to adjust for ---------------------------------------------------
# continuous covariates (entered as standardised values in the weighting step)
cont_vars <- c("agediag", "cci_n_conditions")

# binary covariates, built as 0/1 indicators in 02. The "primary" set adjusts for
# age band, comorbidity, calendar period and season; the "full" set also adjusts
# for sex, ethnicity, deprivation and stage.
bin_primary <- c("yr_late", "q2", "q3", "q4")
bin_full    <- c(bin_primary,
                 "male",
                 "eth_asian", "eth_black", "eth_mixed", "eth_other", "eth_unknown",
                 "imd_2", "imd_3", "imd_4", "imd_5",
                 "stage_2", "stage_3")

# case-mix balance strictness (balancing-weights method) ---------------------
# lambda controls how hard the weights push for an exact match to the overall
# patient mix. 0 = match as closely as possible (can reduce the effective sample
# size sharply); larger = a gentler match that keeps more effective sample size.
# lambda_main is the working value; lambda_grid is scanned to show the trade-off.
lambda_grid <- c(0, .01, .05, .1, .25, .5, 1, 1.5, 2, 2.5, 3)
lambda_main <- 0.05

# Bayesian shrinkage priors (in days) ----------------------------------------
# prior_mu_sd: how far the overall average wait could plausibly sit from the data
#   mean. prior_tau_scale: the typical size of genuine between-hospital
#   differences. Both are weakly informative on the day scale.
prior_mu_sd     <- 50
prior_tau_scale <- 10

# ============================================================================
# PART B  -  helper functions (no need to edit)
# ============================================================================

# standardise a continuous covariate to mean 0, sd 1
z_std <- function(x) (x - mean(x)) / sd(x)

# centre and scale a 0/1 covariate, with a floor so very rare categories do not
# dominate the balance
bin_std <- function(x) {
  p <- mean(x)
  if (p < 0.05) p <- 0.05
  (x - mean(x)) / sqrt(p * (1 - p))
}

# build the standardised covariate matrix from a data frame
make_std_matrix <- function(df, cont, bin) {
  Xc <- sapply(df[cont], z_std)
  Xb <- sapply(df[bin],  bin_std)
  cbind(Xc, Xb)
}

# means / sds / proportions of a target population, used when standardising one
# group to a different population than its own (e.g. the improvement analysis)
ref_moments <- function(df, cont, bin) {
  list(cont_mean = sapply(df[cont], mean),
       cont_sd   = sapply(df[cont], sd),
       bin_p     = sapply(df[bin],  mean))
}

# standardised covariate matrix using externally supplied reference moments, so
# a zero target balances each group to that reference population rather than to
# its own means
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

# pull the per-patient weight out of a balancer standardize() result (the weight
# sits in the patient's own hospital column, so the row max recovers it)
extract_weights <- function(std_out) apply(std_out$weights, 1, max)

# weighted quantile (used for standardised medians); robust to small/degenerate
# cells: ignores zero-weight rows, returns the single value when there is no
# spread, and collapses tied cumulative weights before interpolating
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
# df needs: hosp (id), y (outcome), w (weight), resid (outcome-model residual),
# canonical (population-mean model prediction). Returns one row per hospital with
# the weighted and regression-adjusted estimates, effective sample size and
# pooled standard errors.
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
# shorter wait is better, so rank 1 is the best performer.
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

# resolve recorded hospital codes to canonical hospitals using the optional
# crosswalk; identity (trimmed code) where there is no entry or no file
load_provider_xwalk <- function(path = provider_xwalk_csv) {
  if (!file.exists(path)) {
    message("hospital crosswalk not found at ", path,
            "\n  - using codes as recorded. Curate it (raw_code, canonical_code, ",
            "canonical_name) to merge hub / superseded codes.")
    return(data.frame(raw_code = character(), canonical_code = character(),
                      canonical_name = character()))
  }
  x <- read.csv(path, colClasses = "character", stringsAsFactors = FALSE)
  data.frame(lapply(x, trimws), stringsAsFactors = FALSE)
}

canonicalise_hosp <- function(codes, xwalk) {
  raw <- trimws(toupper(as.character(codes)))
  if (is.null(xwalk) || nrow(xwalk) == 0) return(raw)
  m   <- setNames(xwalk$canonical_code, toupper(xwalk$raw_code))
  out <- unname(m[raw])
  ifelse(is.na(out), raw, out)
}