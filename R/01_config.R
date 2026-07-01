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
stan_file <- file.path(stan_dir, "dp_normal_cont.stan")   # the normal-normal shrinkage model
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

# case-mix adjustment --------------------------------------------------------
# The main analysis adjusts for two clinical factors only: age at diagnosis and
# comorbidity (Charlson). No other patient factors (sex, ethnicity, deprivation,
# stage) are adjusted for. Season and calendar year are held OUT of the main
# model and enter only in the time sensitivity analysis (script 12).
#
# Coding of age and cci is set here and explored in explore_covariate_coding.R.
#   age_coding : "cont" linear age, or "band" (<50 / 50-59 / 60-69 / 70-79 / 80+)
#   cci_coding : "cont" linear Charlson count, "0_1_2plus", "0_1_2_3plus", or "none"
#
# cci_coding is set to "0_1_2plus": it matches how Charlson comorbidity burden is
# naturally described (0, 1, 2+ conditions), and explore_covariate_coding.R shows
# hospital rankings are essentially unchanged across every coding tested (rank
# Spearman >= 0.99 vs continuous), so the natural coding is preferred over the
# marginally higher effective sample size of the continuous term.
age_coding <- "cont"
cci_coding <- "0_1_2plus"

# season (quarter) and calendar-period indicators, added ONLY in the time
# sensitivity analysis, never in the main model
time_bin <- c("q2", "q3", "q4", "yr_late")

# build the case-mix covariates for a given coding: returns the (augmented) data
# and the continuous / binary covariate names to hand to the weighting step.
# reference categories are 60-69 for age bands and 0 conditions for cci.
#
# To add a patient factor to the MAIN model, append its column name(s) to cont
# (continuous, mean-balanced) or bin (0/1 dummies, proportion-balanced) below.
# The dummy columns already exist in the data from 02, so it is a one-liner, e.g.
#   bin <- c(bin, "male")                    # add sex
#   bin <- c(bin, "stage_2", "stage_3")      # add stage (reference = stage 1)
# Stage and sex are deliberately excluded here; add them only if you decide they
# are confounders you want removed rather than part of the pathway.
code_covariates <- function(d, age = age_coding, cci = cci_coding) {
  cont <- character(0); bin <- character(0)
  if (age == "cont") {
    cont <- c(cont, "agediag")
  } else if (age == "band") {
    b <- cut(d$agediag, c(-Inf, 50, 60, 70, 80, Inf),
             labels = c("u50", "50_59", "60_69", "70_79", "80p"))
    for (lv in c("u50", "50_59", "70_79", "80p")) {
      col <- paste0("age_", lv); d[[col]] <- as.integer(b == lv); bin <- c(bin, col)
    }
  }
  if (cci == "cont") {
    cont <- c(cont, "cci_n_conditions")
  } else if (cci == "0_1_2plus") {
    d$cci_1  <- as.integer(d$cci_n_conditions == 1)
    d$cci_2p <- as.integer(d$cci_n_conditions >= 2)
    bin <- c(bin, "cci_1", "cci_2p")
  } else if (cci == "0_1_2_3plus") {
    d$cci_1  <- as.integer(d$cci_n_conditions == 1)
    d$cci_2  <- as.integer(d$cci_n_conditions == 2)
    d$cci_3p <- as.integer(d$cci_n_conditions >= 3)
    bin <- c(bin, "cci_1", "cci_2", "cci_3p")
  }
  list(data = d, cont = cont, bin = bin)
}

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
# Shared building blocks used across the analysis scripts. Each is written to be
# read straight through and checked by hand.

# --- standardising covariates -----------------------------------------------
# The weighting step needs every covariate on a comparable scale. A continuous
# covariate is centred on its mean and scaled by its sd; a 0/1 covariate is
# centred on its proportion and scaled by its binomial sd.

# a continuous covariate, shifted to mean 0 and scaled to sd 1
z_std <- function(x) (x - mean(x)) / sd(x)

# a 0/1 covariate, centred on its proportion p and scaled by sqrt(p(1-p)). The
# 0.05 floor applies to the scale only, so a very rare category is not blown up
# enough to swamp the balance; the centre stays the true proportion.
bin_std <- function(x) {
  p     <- mean(x)
  scale <- sqrt(max(p, 0.05) * (1 - max(p, 0.05)))
  (x - p) / scale
}

# build the standardised covariate matrix: take the covariate columns, z-score
# each continuous one, centre/scale each binary one, and return a numeric matrix
# with one column per covariate. code_covariates() always includes at least age,
# so c(cont, bin) is never empty.
make_std_matrix <- function(df, cont, bin) {
  std <- df[c(cont, bin)]                       # the covariate columns, in order
  for (v in cont) std[[v]] <- z_std(df[[v]])
  for (v in bin)  std[[v]] <- bin_std(df[[v]])
  as.matrix(std)
}

# the mean and sd of each continuous covariate, and the proportion of each binary
# one, in a chosen reference population. Used when a group is standardised to a
# different population than its own (e.g. baseline vs later period in 06).
ref_moments <- function(df, cont, bin) {
  cont_mean <- numeric(0); cont_sd <- numeric(0); bin_p <- numeric(0)
  for (v in cont) {
    cont_mean[v] <- mean(df[[v]])
    cont_sd[v]   <- sd(df[[v]])
  }
  for (v in bin) bin_p[v] <- mean(df[[v]])
  list(cont_mean = cont_mean, cont_sd = cont_sd, bin_p = bin_p)
}

# the same standardised matrix as make_std_matrix, but centred and scaled to
# externally supplied reference moments rather than the group's own values. A
# zero balance target then pulls each group towards that reference population.
make_std_matrix_ref <- function(df, cont, bin, ref) {
  std <- df[c(cont, bin)]
  for (v in cont) {
    std[[v]] <- (df[[v]] - ref$cont_mean[v]) / ref$cont_sd[v]
  }
  for (v in bin) {
    p        <- ref$bin_p[v]
    scale    <- sqrt(max(p, 0.05) * (1 - max(p, 0.05)))
    std[[v]] <- (df[[v]] - p) / scale
  }
  as.matrix(std)
}

# pull the per-patient weight out of a balancer standardize() result. Each
# patient's weight sits in the column for their own hospital, so the row maximum
# recovers it.
extract_weights <- function(std_out) apply(std_out$weights, 1, max)

# weighted quantile, used for standardised medians. Written to be robust in small
# or degenerate cells: drop zero-weight rows, return the single value when there
# is no spread, and collapse tied cumulative weights before interpolating.
w_quantile <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x  <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  if (length(unique(x)) == 1) return(x[1])
  o  <- order(x); x <- x[o]; w <- w[o]
  cum_w <- cumsum(w) / sum(w)
  keep  <- !duplicated(cum_w, fromLast = TRUE)
  cum_w <- cum_w[keep]; x <- x[keep]
  if (length(x) < 2) return(x[length(x)])
  approx(cum_w, x, xout = p, rule = 2, ties = "ordered")$y
}

# per-hospital standardised summaries for a continuous outcome.
# df needs: hosp (id), y (outcome), w (weight), resid (outcome-model residual),
# canonical (population-mean model prediction). Returns one row per hospital with
# the weighted and regression-adjusted estimates, effective sample size and both
# per-hospital and pooled standard errors.
site_summary <- function(df) {
  library(dplyr)
  s <- df %>%
    group_by(hosp) %>%
    summarise(
      n          = n(),
      n_eff      = sum(w)^2 / sum(w^2),                 # effective sample size
      raw_mean   = mean(y),
      raw_median = median(y),
      stand      = weighted.mean(y, w),                 # standardised mean
      stand_med  = w_quantile(y, w, 0.5),               # standardised median
      stand_adj  = weighted.mean(resid, w) + mean(canonical),   # augmented mean
      sd_w       = sqrt(sum(w^2 * (y - weighted.mean(y, w))^2) / sum(w^2)),
      sd_adj     = sqrt(sum(w^2 * resid^2) / sum(w^2)),
      .groups = "drop"
    )
  # pool the within-hospital spread across hospitals (weighted by effective n),
  # giving a more stable SE for small hospitals than each one's own sd
  sd_pool     <- sqrt(weighted.mean(s$sd_w^2,   s$n_eff))
  sd_pool_adj <- sqrt(weighted.mean(s$sd_adj^2, s$n_eff))
  s %>% mutate(
    se          = sd_w        / sqrt(n_eff),
    se_pool     = sd_pool     / sqrt(n_eff),
    se_adj      = sd_adj      / sqrt(n_eff),
    se_adj_pool = sd_pool_adj / sqrt(n_eff),
    sd_pool     = sd_pool,
    sd_pool_adj = sd_pool_adj
  )
}

# posterior ranking metrics from a draws-by-hospital matrix of latent means.
# shorter wait is better, so rank 1 is the best performer. For each posterior
# draw we rank the hospitals, then summarise each hospital's rank across draws.
rank_metrics <- function(draws, tops = c(.05, .10, .20, .25, .50)) {
  J <- ncol(draws)
  ranks <- t(apply(draws, 1, rank, ties.method = "average"))   # draws x hospitals
  out <- data.frame(
    exp_rank = colMeans(ranks),                    # posterior mean rank
    rank_lo  = apply(ranks, 2, quantile, 0.025),   # 95% credible interval
    rank_hi  = apply(ranks, 2, quantile, 0.975)
  )
  # probability each hospital sits in the fastest p% of the ranking
  for (p in tops) {
    threshold <- ceiling(p * J)
    out[[paste0("p_top", p * 100)]] <- colMeans(ranks <= threshold)
  }
  out
}

# --- Bayesian shrinkage (the one routine used everywhere) --------------------
# The normal-normal shrinkage model. Feed it per-hospital point estimates y and
# their standard errors se, and it pulls each hospital towards the overall mean
# by an amount set by its precision. This is the single shrinkage step in the
# analysis: whatever the point estimates are (weighted-standardised means, raw
# means, a change score), they are shrunk the same way, by this function. rstan
# must be loaded by the calling script.
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
  rstan::stan(stan_file, data = dat, seed = seed,
              chains = 4, iter = 4000, warmup = 2000, refresh = refresh, cores = cores,
              control = list(adapt_delta = adapt_delta, max_treedepth = 12))
}

# shrink a set of point estimates and return the shrunk posterior mean and the
# posterior mean rank (1 = fastest). A thin wrapper over fit_shrink for callers
# that only need the shrunk estimate and its rank, in the input hospital order.
stan_shrink_rank <- function(y, se, mu_mean = mean(y)) {
  fit   <- fit_shrink(y, se, mu_mean = mu_mean)
  draws <- rstan::extract(fit, pars = "y_site_true")$y_site_true   # draws x hospitals
  data.frame(post_mean = colMeans(draws), exp_rank = rank_metrics(draws)$exp_rank)
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
  for (v in names(x)) x[[v]] <- trimws(x[[v]])       # tidy stray whitespace
  x
}

# apply the crosswalk: look each recorded code up in the table and return its
# canonical code, keeping the original code where there is no match.
canonicalise_hosp <- function(codes, xwalk) {
  raw <- trimws(toupper(as.character(codes)))
  if (is.null(xwalk) || nrow(xwalk) == 0) return(raw)
  lookup <- setNames(xwalk$canonical_code, toupper(xwalk$raw_code))
  mapped <- lookup[raw]                              # canonical code, or NA
  ifelse(is.na(mapped), raw, unname(mapped))
}

# --- hospital display names -------------------------------------------------
# Shared name tidying, used wherever hospital names are shown (08, 09). Title
# Case a name with small joining words left lower-case (except the first word)
# and NHS kept upper-case; apply manual corrections; and drop any trailing
# "(code)" since the site code has its own column.
name_small_words <- c("and", "of", "the", "for", "in", "on", "at", "to", "by", "an", "a", "or")
name_acronyms    <- c("nhs")
title_case <- function(x) {
  vapply(x, function(one) {
    words <- strsplit(tolower(one), " ")[[1]]
    for (i in seq_along(words)) {
      if (words[i] %in% name_acronyms) {
        words[i] <- toupper(words[i])
      } else if (!(words[i] %in% name_small_words) || i == 1) {
        substr(words[i], 1, 1) <- toupper(substr(words[i], 1, 1))
      }
    }
    paste(words, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

# manual name corrections: replace a whole name (e.g. a crosswalk entry that
# carries the trust name rather than the site) and fix a spelling error that
# comes through from the provider master. Both are best corrected at source.
name_fixes <- c(
  "Tameside and Glossop Integrated Care NHS Foundation Trust" = "Tameside General Hospital"
)
fix_names <- function(x) {
  for (bad in names(name_fixes)) x[x == bad] <- name_fixes[[bad]]
  gsub("Westminister", "Westminster", x)
}

strip_code_suffix <- function(x) trimws(sub("\\s*\\([^)]*\\)\\s*$", "", x))

# named vector: canonical site code -> cleaned display name, from the diagnosing
# crosswalk; NULL if the crosswalk is not present.
hospital_names <- function(path = diag_xwalk_csv) {
  if (!file.exists(path)) return(NULL)
  nm <- read.csv(path, colClasses = "character")
  nm <- nm[!duplicated(nm$canonical_code), ]
  setNames(fix_names(title_case(strip_code_suffix(nm$canonical_name))), nm$canonical_code)
}