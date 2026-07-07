# 09  ranking-method comparison and funnel plots
# -----------------------------------------------------------------------------
# Compares hospital rankings on the sustained and improvement estimands (the
# basis for Table 3), draws funnel plots, and builds rank-movement tables.
#
# The comparison is deliberately built from two clean steps, so every method is
# on the same footing:
#   step 1 - generate a per-hospital mean waiting time, either case-mix adjusted
#            (balancing-weights direct standardisation, the primary model) or not
#            adjusted at all (raw mean);
#   step 2 - optionally shrink those means with the one Bayesian normal-normal
#            routine (fit_shrink, in 01_config.R) - the same routine that feeds
#            the headline figures.
# So each generator appears twice, unshrunk and shrunk, and the shrinkage is
# identical across them. Indirect standardisation (observed minus expected) is
# also shown, as the basis for the funnel plots.
#
# Funnel columns. The movement table carries three funnel-plot columns, all on
# the precision-weighted funnel z-score (deviation from target over its standard
# error): a standard funnel, and two that correct for over-dispersion with an
# additive random-effects model after Spiegelhalter (2005) - one with the
# between-hospital variance tau^2 by method of moments, one with tau^2 estimated
# from Winsorised residuals so a few extreme hospitals do not inflate it. The
# random-effects correction adds tau^2 to each hospital's own variance; because
# that variance is larger for low-volume hospitals, adding a constant tau^2
# shrinks the high-volume z-scores more, so the ranking genuinely changes -
# unlike a multiplicative factor, which rescales every z equally and cannot
# reorder. The two tau estimates differ, so the two random-effects columns give
# genuinely different orderings.
# Set base_col to "wt_mean" to reference the un-shrunk balancing weights instead
# of the headline (weighted + shrinkage) column.

library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)
library(flextable)
library(officer)

source("R/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

df    <- readRDS(file.path(out_dir, "analysis_data.rds"))
fit_p <- readRDS(file.path(out_dir, "fit_primary.rds"))
sus   <- read.csv(file.path(out_dir, "ranks_sustained.csv"))   # headline weighted + shrinkage

cv   <- code_covariates(df)          # age + comorbidity
df   <- cv$data
cont <- cv$cont; bin <- cv$bin
# the primary adjustment set (age + cci + season + calendar year) is what the
# headline weighted column uses, so the regression and indirect methods below
# standardise on the same set - the comparison is then of method, not of set.
# For the change score the later-year indicator coincides with the period split,
# so only season is added there (year would alias the hospital-half terms).
bin_primary <- c(bin, season_terms, year_term)
bin_change  <- c(bin, season_terms)

# funnel statistics for a continuous outcome (Spiegelhalter 2005, adapted) -----
# Given each hospital's indicator value, the common target, and the standard
# error of the indicator, returns the standard funnel z-score and a random-
# effects adjusted z-score, plus the significance flag under each. A hospital is
# flagged significantly FAST (below the lower funnel limit) with ** at the 99.8%
# limit and * at the 95% limit.
#
# The between-hospital variance tau^2 is estimated by the DerSimonian-Laird
# method of moments, taken about the precision-weighted mean so it captures
# genuine heterogeneity and not any offset between that mean and the target (this
# matters for the change score, where the target is zero but the mean change is
# not). If the heterogeneity statistic falls below its null expectation, tau^2 is
# held at zero. A robust tau (Winsorising the standardised residuals before
# forming the statistic) is also returned: the paper recommends it so that a few
# extreme hospitals do not inflate tau and let the widened funnel accommodate the
# very hospitals one is trying to detect. Two random-effects z-scores are
# returned, one using the method-of-moments tau and one the robust tau.
funnel_stats <- function(value, target, se, winsor_q = 0.1) {
  n_units <- length(value)
  w       <- 1 / se^2
  w_mean  <- sum(w * value) / sum(w)               # precision-weighted centre
  z       <- (value - target) / se                 # funnel z-score about the target
  phi     <- sum(z^2) / n_units                     # multiplicative factor, for reference
  
  denom   <- sum(w) - sum(w^2) / sum(w)             # usual method-of-moments denominator
  het     <- sum(w * (value - w_mean)^2)            # heterogeneity statistic (Cochran Q)
  tau2    <- max((het - (n_units - 1)) / denom, 0)
  
  resid_std <- (value - w_mean) * sqrt(w)           # standardised residual about the centre
  r_lo <- as.numeric(quantile(resid_std, winsor_q))
  r_hi <- as.numeric(quantile(resid_std, 1 - winsor_q))
  resid_win <- pmin(pmax(resid_std, r_lo), r_hi)    # Winsorised for the robust tau
  tau2_rob  <- max((sum(resid_win^2) - (n_units - 1)) / denom, 0)
  
  z_re     <- (value - target) / sqrt(se^2 + tau2)      # random-effects z, method of moments
  z_re_rob <- (value - target) / sqrt(se^2 + tau2_rob)  # random-effects z, robust (Winsorised) tau
  
  fast_flag <- function(mod_z) ifelse(mod_z < -3.09, "**", ifelse(mod_z < -1.96, "*", ""))
  data.frame(hosp_z      = z,
             z_re        = z_re,
             z_re_rob    = z_re_rob,
             flag_std    = fast_flag(z),
             flag_re     = fast_flag(z_re),
             flag_re_rob = fast_flag(z_re_rob),
             phi = phi, tau = sqrt(tau2), tau_robust = sqrt(tau2_rob),
             stringsAsFactors = FALSE)
}

# model-based direct standardisation (g-computation) for the sustained estimand.
# Fit one linear model with a separate term for every hospital (no intercept, so
# each hospital gets its own coefficient) plus the same primary-set covariates
# the weighting uses. A hospital's directly-standardised mean is the average
# predicted wait if the whole sample were treated at that hospital. Because the
# model is linear and the covariates enter additively, that average is simply:
# (this hospital's own coefficient) + (the population-mean covariate vector)
# times (the covariate coefficients). Writing that as a linear combination L of
# the coefficients gives the mean as L %*% coef and its standard error as
# sqrt(L %*% vcov %*% L') - the exact model-based se, with no bootstrap.
model_based_standardise <- function(patient_data, continuous_covariates, binary_covariates) {
  data       <- patient_data %>% mutate(hospital = factor(hosp))
  covariates <- c(continuous_covariates, binary_covariates)
  
  model <- lm(as.formula(paste("wait ~ 0 + hospital +", paste(covariates, collapse = " + "))),
              data = data)
  coefficients    <- coef(model)
  covariance      <- vcov(model)
  hospital_levels <- levels(data$hospital)
  covariate_means <- colMeans(data[covariates])       # the population case-mix
  
  # one row of L per hospital: 1 on its own hospital coefficient, and the
  # population covariate means on the covariate coefficients
  L <- matrix(0, nrow = length(hospital_levels), ncol = length(coefficients),
              dimnames = list(hospital_levels, names(coefficients)))
  for (i in seq_along(hospital_levels)) {
    L[i, paste0("hospital", hospital_levels[i])] <- 1
    L[i, names(covariate_means)]                 <- covariate_means
  }
  
  data.frame(hosp     = as.integer(hospital_levels),
             reg_mean = as.numeric(L %*% coefficients),
             reg_se   = sqrt(diag(L %*% covariance %*% t(L))))
}

# The improvement analogue: the change in a hospital's model-based standardised
# mean between the two halves. Fit one coefficient per hospital-and-half (again
# no intercept) plus the covariates; the change for a hospital is its later-half
# coefficient minus its earlier-half coefficient. The additive covariate part is
# identical in both halves and cancels in the difference, so no covariate means
# are needed here. The se again comes from that contrast's model covariance.
model_based_change <- function(patient_data, continuous_covariates, binary_covariates) {
  data       <- patient_data %>% mutate(hospital_half = factor(paste(hosp, period, sep = "_")))
  covariates <- c(continuous_covariates, binary_covariates)
  
  model <- lm(as.formula(paste("wait ~ 0 + hospital_half +", paste(covariates, collapse = " + "))),
              data = data)
  coefficients <- coef(model)
  covariance   <- vcov(model)
  hospitals    <- sort(unique(data$hosp))
  
  # one row of L per hospital: +1 on its later-half term, -1 on its earlier-half term
  L <- matrix(0, nrow = length(hospitals), ncol = length(coefficients),
              dimnames = list(hospitals, names(coefficients)))
  for (i in seq_along(hospitals)) {
    L[i, paste0("hospital_half", hospitals[i], "_second")] <-  1
    L[i, paste0("hospital_half", hospitals[i], "_first")]  <- -1
  }
  
  data.frame(hosp     = hospitals,
             reg_mean = as.numeric(L %*% coefficients),
             reg_se   = sqrt(diag(L %*% covariance %*% t(L))))
}

# indirect standardisation: expected wait from a pooled case-mix model, then
# observed - expected + grand mean. Also the basis for the sustained funnel plot.
pm <- lm(as.formula(paste("wait ~", paste(c(cont, bin_primary), collapse = " + "))),
         data = df)
df$pred <- predict(pm)
grand   <- mean(df$wait)
sigma_r <- sd(resid(pm))             # residual sd, used for the funnel standard error
ind <- df %>% group_by(hosp) %>%
  summarise(n = n(), obs = mean(wait), exp = mean(pred), .groups = "drop") %>%
  mutate(indirect = obs - exp + grand)

# funnel z-scores (standard and random-effects) for the sustained indicator; the
# se is the residual sd over root-n, matching the sustained funnel plot below.
fs_sus      <- funnel_stats(ind$indirect, grand, sigma_r / sqrt(ind$n))
tau_sus     <- fs_sus$tau[1]
tau_sus_rob <- fs_sus$tau_robust[1]
ind$funnel_z        <- fs_sus$hosp_z
ind$funnel_z_re     <- fs_sus$z_re
ind$funnel_z_re_rob <- fs_sus$z_re_rob
ind$flag_std        <- fs_sus$flag_std
ind$flag_re         <- fs_sus$flag_re
ind$flag_re_rob     <- fs_sus$flag_re_rob

# raw standard error per hospital, for the no-adjustment shrinkage. The within-
# hospital sd is pooled across hospitals (weighted by n) so small hospitals get a
# stable se, matching how the weighted estimate's pooled se is built in 06.
raw_sd <- df %>% group_by(hosp) %>% summarise(nn = n(), sdw = sd(wait), .groups = "drop")
sd_pool_raw <- sqrt(weighted.mean(raw_sd$sdw^2, raw_sd$nn))
raw_sd <- raw_sd %>% mutate(se_raw = sd_pool_raw / sqrt(nn))

# sustained point estimates per hospital: the augmented weighted mean (the main
# model estimate) and the raw mean, both in hospital order.
site <- fit_p$site %>% arrange(hosp) %>%
  left_join(raw_sd %>% select(hosp, se_raw), by = "hosp")

# step 2 for the raw means: shrink with the same normal-normal routine as the
# headline. Keep the shrunk posterior mean (for the method scatter) and its rank.
raw_shr <- stan_shrink_rank(site$raw_mean, site$se_raw)
site$raw_shrunk_rank <- raw_shr$exp_rank
site$raw_shrunk_mean <- raw_shr$post_mean

# model-based direct standardisation, shrunk with the same routine again
reg_sus <- model_based_standardise(patient_data = df, continuous_covariates = cont, binary_covariates = bin_primary)
reg_shr_sus <- stan_shrink_rank(reg_sus$reg_mean, reg_sus$reg_se)
reg_sus$reg_shrunk_rank <- reg_shr_sus$exp_rank
reg_sus$reg_shrunk_mean <- reg_shr_sus$post_mean

comp <- site %>%
  transmute(hosp, diag_hosp,
            wt_mean = stand_adj,       # weighted direct-standardised mean (unshrunk)
            raw_mean,                  # raw mean (unshrunk)
            raw_shrunk_rank) %>%       # raw mean, shrunk
  left_join(reg_sus %>% select(hosp, reg = reg_mean, reg_shrunk_rank), by = "hosp") %>%
  left_join(ind %>% select(hosp, funnel_z, funnel_z_re, funnel_z_re_rob,
                           flag_std, flag_re, flag_re_rob), by = "hosp") %>%
  left_join(sus %>% select(hosp, shrunk_rank = exp_rank), by = "hosp")   # weighted, shrunk (headline)

# rank agreement between methods (1 = fastest); the shrunk columns are already
# posterior ranks, so rank() just re-expresses them on the same scale. Both
# funnel z-scores enter, so the correlation table shows how much the random-
# effects correction moves the ordering relative to the plain funnel.
methods <- c("raw_mean", "wt_mean", "reg", "funnel_z", "funnel_z_re", "funnel_z_re_rob")
ranks <- comp %>%
  mutate(across(all_of(methods), ~ rank(.x, ties.method = "average"), .names = "rank_{.col}"),
         rank_wt_shrunk    = rank(shrunk_rank,     ties.method = "average"),
         rank_model_shrunk = rank(reg_shrunk_rank, ties.method = "average"),
         rank_raw_shrunk   = rank(raw_shrunk_rank, ties.method = "average"))
rank_cols <- c(paste0("rank_", methods), "rank_wt_shrunk", "rank_model_shrunk", "rank_raw_shrunk")
rho <- cor(ranks[rank_cols], method = "spearman")
cat("Spearman rank correlation between methods:\n"); print(round(rho, 2))

write.csv(comp, file.path(out_dir, "method_comparison_estimates.csv"), row.names = FALSE)
write.csv(ranks %>% select(hosp, diag_hosp, all_of(rank_cols)),
          file.path(out_dir, "method_comparison_ranks.csv"), row.names = FALSE)
write.csv(round(rho, 3), file.path(out_dir, "method_comparison_rho.csv"))

# rank-movement tables -------------------------------------------------------
# For each estimand, take the hospitals the PRIMARY model ranks best, then show
# where every other method places them, with an up/down arrow giving the rank
# change against the primary. The primary is the balancing-weights model with
# Bayesian shrinkage - the ranking that feeds the main figures. Sustained and
# improvement are shown as one continuous Word table. The rank-movement helpers
# (comp_rank, move_cell, move_colour, cell_prefix/suffix) are shared in config.
TABLE_FONT_SIZE <- 7               # one place, so plain and coloured text match

# hospital display names (Title Case, corrected, code suffix stripped) come from
# the shared helper in 01_config.R, keyed by canonical site code.
hosp_name <- hospital_names()

# shared column set: display title -> estimate column, in table order. Three
# direct-standardisation generators (weighted, model-based, and none), each shown
# unshrunk and shrunk, then two funnel columns: standard, and random-effects.
method_cols <- c(
  "weighted mean + shrinkage" = "shrunk_rank",
  "weighted mean"             = "wt_mean",
  "model mean + shrinkage"    = "reg_shrunk_rank",
  "model mean"                = "reg",
  "raw mean + shrinkage"      = "raw_shrunk_rank",
  "raw mean"                  = "raw_mean",
  "funnel"                    = "funnel_z",
  "funnel RE"                 = "funnel_z_re",
  "funnel RE Wins."           = "funnel_z_re_rob"
)
base_col <- "shrunk_rank"

# which flag column carries each funnel column's significance marker
funnel_flag_of <- c("funnel"          = "flag_std",
                    "funnel RE"       = "flag_re",
                    "funnel RE Wins." = "flag_re_rob")

# format the top-n hospitals for one estimand. Returns the display text (Hospital
# name, site code, one column per method) and the matching matrix of signed moves
# used to colour the cells. The funnel significance marker is appended after the
# rank number so it stays black and the coloured move suffix is untouched.
format_block <- function(est, top_n = 20) {
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(est))
  else unname(hosp_name[as.character(est$diag_hosp)])
  est$hosp_name <- ifelse(is.na(nm) | nm == "", as.character(est$diag_hosp), nm)
  est$hosp_code <- as.character(est$diag_hosp)
  
  present <- method_cols[method_cols %in% names(est)]
  rk <- data.frame(hosp = est$hosp)
  for (col in unique(present)) rk[[col]] <- comp_rank(est[[col]])
  base_rk <- rk[[base_col]]
  keep    <- order(base_rk, est$hosp)[seq_len(min(top_n, nrow(est)))]
  
  text <- data.frame(Hospital = est$hosp_name[keep], `Hospital site code` = est$hosp_code[keep],
                     stringsAsFactors = FALSE, check.names = FALSE)
  move <- matrix(NA_real_, nrow = length(keep), ncol = length(method_cols),
                 dimnames = list(NULL, names(method_cols)))
  for (lab in names(method_cols)) {
    col <- method_cols[[lab]]
    if (!(col %in% names(est))) { text[[lab]] <- "-"; next }
    mr <- rk[[col]][keep]
    text[[lab]] <- if (col == base_col) as.character(mr) else move_cell(mr, base_rk[keep])
    move[, lab] <- if (col == base_col) 0 else base_rk[keep] - mr
    if (lab %in% names(funnel_flag_of)) {
      fl <- est[[ funnel_flag_of[[lab]] ]][keep]
      for (r in seq_along(fl)) {
        if (!is.na(fl[r]) && fl[r] != "") {
          text[[lab]][r] <- sub("^([0-9]+)", paste0("\\1", fl[r]), text[[lab]][r])
        }
      }
    }
  }
  list(text = text, move = move)
}

# improvement estimand: each generator recomputed as a change score -----------
imp_wt      <- readRDS(file.path(out_dir, "site_improve.rds"))    # weighted change (delta)
impr_shrunk <- read.csv(file.path(out_dir, "ranks_improve.csv"))  # weighted change, shrunk

di <- df %>% filter(hosp %in% imp_wt$hosp) %>%
  mutate(periodf = factor(period, levels = c("first", "second")))

# raw change and its se: later-period minus first-period mean, with the within-
# period sd pooled across hospital-periods for a stable se (as for the sustained
# raw se). se of the change adds the two period ses in quadrature.
raw_imp <- di %>% group_by(hosp, period) %>%
  summarise(m = mean(wait), sdw = sd(wait), nn = n(), .groups = "drop")
sd_pool_imp <- sqrt(weighted.mean(raw_imp$sdw^2, raw_imp$nn))
raw_change <- raw_imp %>%
  mutate(se = sd_pool_imp / sqrt(nn)) %>%
  select(hosp, period, m, se) %>%
  pivot_wider(names_from = period, values_from = c(m, se)) %>%
  transmute(hosp, raw_mean = m_second - m_first,
            se_raw = sqrt(se_first^2 + se_second^2))

# step 2 for the raw change: shrink with the same routine, centred on no change.
raw_shr_imp <- stan_shrink_rank(raw_change$raw_mean, raw_change$se_raw, mu_mean = 0)
raw_change$raw_shrunk_rank <- raw_shr_imp$exp_rank
raw_change$raw_shrunk_mean <- raw_shr_imp$post_mean

# model-based change, shrunk with the same routine (centred on no change)
reg_imp <- model_based_change(patient_data = di, continuous_covariates = cont, binary_covariates = bin_change)
reg_shr_imp <- stan_shrink_rank(reg_imp$reg_mean, reg_imp$reg_se, mu_mean = 0)
reg_imp$reg_shrunk_rank <- reg_shr_imp$exp_rank
reg_imp$reg_shrunk_mean <- reg_shr_imp$post_mean

# improvement funnel: run the indirect calculation independently in each period,
# then difference. Within a period the indirect mean is observed minus expected
# (expected from the pooled case-mix model), so the later-minus-earlier difference
# is the indirect-standardised change in days (the grand mean cancels). Each
# period's indirect mean has standard error sigma_r / sqrt(n_period), so the
# difference has se = sigma_r * sqrt(1/n1 + 1/n2) = sigma_r / sqrt(n_eff) with the
# effective volume n_eff = n1 n2 / (n1 + n2). funnel_stats then gives the standard
# and over-dispersion (random-effects and Winsorised) z-scores of that difference.
oe_period <- di %>% group_by(hosp, period) %>%
  summarise(oe = mean(wait) - mean(pred), n = n(), .groups = "drop")
imp_fun <- oe_period %>%
  pivot_wider(names_from = period, values_from = c(oe, n)) %>%
  transmute(hosp,
            indirect = oe_second - oe_first,
            n        = n_first + n_second,
            n_eff    = (n_first * n_second) / (n_first + n_second),
            se_diff  = sigma_r * sqrt(1 / n_first + 1 / n_second))
fs_imp      <- funnel_stats(imp_fun$indirect, 0, imp_fun$se_diff)
tau_imp     <- fs_imp$tau[1]
tau_imp_rob <- fs_imp$tau_robust[1]
imp_fun$funnel_z        <- fs_imp$hosp_z
imp_fun$funnel_z_re     <- fs_imp$z_re
imp_fun$funnel_z_re_rob <- fs_imp$z_re_rob
imp_fun$flag_std        <- fs_imp$flag_std
imp_fun$flag_re         <- fs_imp$flag_re
imp_fun$flag_re_rob     <- fs_imp$flag_re_rob

comp_imp <- imp_wt %>%
  transmute(hosp, diag_hosp, wt_mean = delta) %>%
  left_join(raw_change %>% select(hosp, raw_mean, raw_shrunk_rank), by = "hosp") %>%
  left_join(reg_imp %>% select(hosp, reg = reg_mean, reg_shrunk_rank), by = "hosp") %>%
  left_join(imp_fun %>% select(hosp, funnel_z, funnel_z_re, funnel_z_re_rob,
                               flag_std, flag_re, flag_re_rob), by = "hosp") %>%
  left_join(impr_shrunk %>% select(hosp, shrunk_rank = exp_rank), by = "hosp")

# per-hospital data for the method-comparison scatter figure (script 17). The
# y-axis there is the main model's shrunk posterior mean; the x-axes are the other
# generators' shrunk means (regression, no adjustment), the funnel-plot method's
# indirectly standardised mean (days), and the three funnel z-scores.
scatter_sustained <- site %>%
  select(hosp, diag_hosp, raw_shrunk_mean) %>%
  left_join(reg_sus %>% select(hosp, reg_shrunk_mean), by = "hosp") %>%
  left_join(ind %>% select(hosp, indirect, funnel_z, funnel_z_re, funnel_z_re_rob), by = "hosp") %>%
  left_join(sus %>% select(hosp, main_shrunk_mean = post_mean), by = "hosp")
write.csv(scatter_sustained, file.path(out_dir, "method_scatter_sustained.csv"), row.names = FALSE)

scatter_improve <- raw_change %>%
  select(hosp, raw_shrunk_mean) %>%
  left_join(reg_imp %>% select(hosp, reg_shrunk_mean), by = "hosp") %>%
  left_join(imp_fun %>% select(hosp, indirect, funnel_z, funnel_z_re, funnel_z_re_rob), by = "hosp") %>%
  left_join(impr_shrunk %>% select(hosp, diag_hosp, main_shrunk_mean = post_mean), by = "hosp")
write.csv(scatter_improve, file.path(out_dir, "method_scatter_improve.csv"), row.names = FALSE)

# over-dispersion summary: the estimated between-hospital sd (method of moments
# and the robust Winsorised version), the multiplicative factor for reference,
# and how many hospitals the standard and random-effects funnels flag as fast.
fun_summary <- function(fs, label) {
  data.frame(estimand = label,
             phi        = round(fs$phi[1], 2),
             tau        = round(fs$tau[1], 2),
             tau_robust = round(fs$tau_robust[1], 2),
             flagged_standard          = sum(fs$flag_std    != ""),
             flagged_random_effect     = sum(fs$flag_re     != ""),
             flagged_random_effect_rob = sum(fs$flag_re_rob != ""),
             stringsAsFactors = FALSE)
}
fun_class <- rbind(fun_summary(fs_sus, "sustained"), fun_summary(fs_imp, "improvement"))
write.csv(fun_class, file.path(out_dir, "funnel_overdispersion_summary.csv"), row.names = FALSE)
cat("\nFunnel over-dispersion summary:\n"); print(fun_class, row.names = FALSE)

# build the two blocks (top 20 each) and write the plain CSVs -----------------
sus_block <- format_block(comp,     top_n = 20)
imp_block <- format_block(comp_imp, top_n = 20)
write.csv(sus_block$text, file.path(out_dir, "movement_sustained.csv"), row.names = FALSE)
write.csv(imp_block$text, file.path(out_dir, "movement_improve.csv"),   row.names = FALSE)
cat("\nTop 20 hospitals, sustained:\n");   print(sus_block$text, row.names = FALSE)
cat("\nTop 20 hospitals, improvement:\n"); print(imp_block$text, row.names = FALSE)

# stitch into one table, all as plain body rows so the broad groupings appear
# only once: the broad group headings, a "sustained performance" banner, the
# column titles, the sustained rows, an "improvement" banner, the column titles
# again, then the improvement rows. The caption and footnote sit outside the
# table as ordinary text.
col_names   <- names(sus_block$text)                # "Hospital","Hospital site code", the methods
funnel_cols <- c("funnel", "funnel RE", "funnel RE Wins.")

blank_row <- function(first_cell = "") {
  r <- as.data.frame(as.list(rep("", length(col_names))),
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(r) <- col_names
  r[[1]] <- first_cell
  r
}
group_row <- blank_row("")
group_row[["weighted mean + shrinkage"]] <- "Balancer weighted direct standardisation"
group_row[["model mean + shrinkage"]]    <- "Regression-based direct standardisation"
group_row[["raw mean + shrinkage"]]      <- "No adjustment"
group_row[["funnel"]]                    <- "Funnel plot"

titles_top <- as.data.frame(as.list(col_names), stringsAsFactors = FALSE, check.names = FALSE)
names(titles_top) <- col_names
titles_imp <- titles_top

sus_divider <- blank_row("Sustained performance (average waiting time, 2020-2021)")
imp_divider <- blank_row("Improvement over the period (change, 2021 vs 2020)")

combined <- rbind(group_row, sus_divider, titles_top, sus_block$text,
                  imp_divider, titles_imp, imp_block$text)

n_sus <- nrow(sus_block$text); n_imp <- nrow(imp_block$text)
row_group   <- 1
row_sus_div <- 2
row_titles1 <- 3
sus_rows    <- row_titles1 + seq_len(n_sus)
row_imp_div <- max(sus_rows) + 1
row_titles2 <- row_imp_div + 1
imp_rows    <- row_titles2 + seq_len(n_imp)

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")   # every row above is already in the body

# colour only the "(arrow N)" part of each cell, one row at a time so the rank
# number (and any funnel marker) stays plain black and the same size.
colour_moves <- function(ft, block_text, block_move, rows) {
  for (lab in names(method_cols)) {
    for (r in seq_along(rows)) {
      cell <- block_text[[lab]][r]
      mv   <- block_move[r, lab]
      if (is.na(mv) || mv == 0 || !grepl("\\(", cell)) next   # "-", primary, or no change
      ft <- compose(ft, i = rows[r], j = lab, part = "body",
                    value = as_paragraph(
                      as_chunk(cell_prefix(cell),
                               props = fp_text(color = "black", font.size = TABLE_FONT_SIZE)),
                      as_chunk(cell_suffix(cell),
                               props = fp_text(color = move_colour(mv), font.size = TABLE_FONT_SIZE))
                    ))
    }
  }
  ft
}
ft <- colour_moves(ft, sus_block$text, sus_block$move, sus_rows)
ft <- colour_moves(ft, imp_block$text, imp_block$move, imp_rows)

# merge each broad heading and each banner row across the columns it covers. No
# fill colour anywhere; sections are shown with rule lines instead.
ft <- merge_at(ft, i = row_group, j = match(c("Hospital", "Hospital site code"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("weighted mean + shrinkage", "weighted mean"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("model mean + shrinkage", "model mean"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("raw mean + shrinkage", "raw mean"), col_names), part = "body")
ft <- merge_at(ft, i = row_group, j = match(funnel_cols, col_names), part = "body")
ft <- merge_at(ft, i = row_sus_div, j = seq_along(col_names), part = "body")
ft <- merge_at(ft, i = row_imp_div, j = seq_along(col_names), part = "body")

ft <- bold(ft, i = c(row_group, row_sus_div, row_titles1, row_imp_div, row_titles2), part = "body")
ft <- align(ft, part = "body", align = "center")
ft <- align(ft, j = 1:2, align = "left", part = "body")
ft <- align(ft, i = c(row_sus_div, row_imp_div), align = "left", part = "body")

# horizontal rules marking the section breaks; the outer frame already gives a
# line above the broad-headings row and below the last row. No background fill.
rule <- fp_border(color = "black", width = 1)
ft <- border_outer(ft, border = rule, part = "body")
ft <- hline(ft, i = row_group,     border = rule, part = "body")
ft <- hline(ft, i = row_sus_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles1,   border = rule, part = "body")
ft <- hline(ft, i = max(sus_rows), border = rule, part = "body")
ft <- hline(ft, i = row_imp_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles2,   border = rule, part = "body")

# very faint grey rules between individual hospitals, so rows are easy to follow
# across the width; interior rows only, since each block's last row already has a
# black rule (sustained) or the outer frame (improvement).
faint <- fp_border(color = "grey85", width = 0.5)
ft <- hline(ft, i = sus_rows[-length(sus_rows)], border = faint, part = "body")
ft <- hline(ft, i = imp_rows[-length(imp_rows)], border = faint, part = "body")

# faint dotted vertical rules splitting the broad sections (and the id columns
# from the methods); only on the structured rows, so they do not cut through the
# merged banner rows.
vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
vcols <- match(c("Hospital site code", "weighted mean", "model mean", "raw mean"), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

# column widths sized to fit landscape: narrow the ranking columns and let their
# titles wrap; keep Hospital as wide as fits so most names stay on one line.
ft <- autofit(ft)
ft <- width(ft, j = "Hospital", width = 1.9)              # +12%
ft <- width(ft, j = "Hospital site code", width = 0.63)  # +15%
ft <- width(ft, j = names(method_cols), width = 0.78)   # ranking columns, +30%
ft <- set_table_properties(ft, layout = "fixed", align = "left")

# caption above the table and a footnote below explaining the funnel columns.
caption_text <- paste(
  "Hospital rankings by method: sustained (upper block) and improvement (lower",
  "block). Each cell gives the method's rank and its change against the primary",
  "model; the change is coloured green (up) to red (down).")
footnote_text <- paste(
  "Funnel columns rank hospitals by the precision-weighted funnel z-score.",
  "'funnel' is the standard funnel; 'funnel RE' adds an additive random-effects",
  "over-dispersion correction (between-hospital sd tau by method of moments);",
  "'funnel RE Wins.' uses the same correction with tau estimated from Winsorised",
  "residuals, so a few extreme hospitals do not inflate it. The random-effects",
  "corrections reorder hospitals because they shrink high-volume z-scores more",
  "than low-volume ones. Markers flag a hospital significantly fast: ** below the",
  "99.8% funnel limit, * below the 95% limit.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.5, right = 0.5)))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "movement_methods.docx"))

# funnel plot, sustained -----------------------------------------------------
# standard 95% and 99.8% limits, plus the wider random-effects 95% band that
# adds the between-hospital sd tau to each hospital's standard error.
fp <- ind %>% mutate(
  se_i  = sigma_r / sqrt(n),
  lo95  = grand - 1.96 * se_i,
  hi95  = grand + 1.96 * se_i,
  lo998 = grand - 3.09 * se_i,
  hi998 = grand + 3.09 * se_i,
  lo95_re = grand - 1.96 * sqrt(se_i^2 + tau_sus^2),
  hi95_re = grand + 1.96 * sqrt(se_i^2 + tau_sus^2),
  lo95_rr = grand - 1.96 * sqrt(se_i^2 + tau_sus_rob^2),
  hi95_rr = grand + 1.96 * sqrt(se_i^2 + tau_sus_rob^2),
  fast = indirect < lo95_rr)
ggsave(file.path(out_dir, "funnel_sustained.pdf"),
       ggplot(fp, aes(n, indirect)) +
         geom_hline(yintercept = grand, colour = "grey50") +
         geom_line(aes(y = lo95),  linetype = 2, colour = "grey40") +
         geom_line(aes(y = hi95),  linetype = 2, colour = "grey40") +
         geom_line(aes(y = lo998), linetype = 3, colour = "grey60") +
         geom_line(aes(y = hi998), linetype = 3, colour = "grey60") +
         geom_line(aes(y = lo95_re), colour = "steelblue") +
         geom_line(aes(y = hi95_re), colour = "steelblue") +
         geom_line(aes(y = lo95_rr), colour = "darkorange") +
         geom_line(aes(y = hi95_rr), colour = "darkorange") +
         geom_point(aes(colour = fast), size = 1) +
         scale_colour_manual(values = c("black", "firebrick"), guide = "none") +
         labs(x = "Hospital volume", y = "Indirectly standardised days to DTT",
              title = "Funnel plot, sustained") +
         theme_bw(),
       width = 7, height = 5)

# funnel plot, improvement ---------------------------------------------------
fimp <- imp_fun %>% mutate(
  se_i = se_diff,
  lo95 = 0 - 1.96 * se_i,
  hi95 = 0 + 1.96 * se_i,
  lo95_re = 0 - 1.96 * sqrt(se_i^2 + tau_imp^2),
  hi95_re = 0 + 1.96 * sqrt(se_i^2 + tau_imp^2),
  lo95_rr = 0 - 1.96 * sqrt(se_i^2 + tau_imp_rob^2),
  hi95_rr = 0 + 1.96 * sqrt(se_i^2 + tau_imp_rob^2),
  fast = indirect < lo95_rr)
ggsave(file.path(out_dir, "funnel_improve.pdf"),
       ggplot(fimp, aes(n_eff, indirect)) +
         geom_hline(yintercept = 0, colour = "grey50") +
         geom_line(aes(y = lo95), linetype = 2, colour = "grey40") +
         geom_line(aes(y = hi95), linetype = 2, colour = "grey40") +
         geom_line(aes(y = lo95_re), colour = "steelblue") +
         geom_line(aes(y = hi95_re), colour = "steelblue") +
         geom_line(aes(y = lo95_rr), colour = "darkorange") +
         geom_line(aes(y = hi95_rr), colour = "darkorange") +
         geom_point(aes(colour = fast), size = 1) +
         scale_colour_manual(values = c("black", "firebrick"), guide = "none") +
         labs(x = "Effective volume, n1 n2 / (n1 + n2)",
              y = "Indirectly standardised change (days)",
              title = "Funnel plot, improvement") +
         theme_bw(),
       width = 7, height = 5)