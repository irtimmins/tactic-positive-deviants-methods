# explore_covariate_coding  choose how to code age and cci (run once, not in the
# pipeline). Two questions: (1) how should age and cci enter the outcome model of
# waiting time - is a linear age term enough or do a quadratic / spline earn their
# place, and does a finer comorbidity coding help - judged by AIC and a
# likelihood-ratio test against a simple base model; and (2) do hospital ranks
# actually move between the codings the pipeline can use. If ranks are stable,
# prefer the simpler, higher effective-sample-size coding and report robustness.

library(balancer)
library(dplyr)
library(splines)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# 1. descriptive: mean wait across the natural categories ---------------------
# uneven step changes between adjacent levels are a hint that a linear term will
# not capture the relationship.
cat("Mean wait by age band:\n")
df %>% mutate(band = cut(agediag, c(-Inf, 50, 60, 70, 80, Inf),
                         labels = c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  group_by(band) %>% summarise(n = n(), mean_wait = mean(wait), .groups = "drop") %>%
  mutate(step = mean_wait - lag(mean_wait)) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nMean wait by number of Charlson conditions:\n")
df %>% mutate(cci = pmin(cci_n_conditions, 5)) %>%
  group_by(cci) %>% summarise(n = n(), mean_wait = mean(wait), .groups = "drop") %>%
  mutate(step = mean_wait - lag(mean_wait)) %>% as.data.frame() %>% print(row.names = FALSE)

# 2. outcome-model coding: does a richer term earn its degrees of freedom? -----
# Start from a simple base model - linear age and comorbidity as 0 / 1 / 2+ - and
# ask whether each richer coding improves the fit of wait on age + cci. Every
# candidate below adds terms that CONTAIN the base as a special case (a straight
# line is a natural spline; a coarse factor is a collapsed fine factor), so a
# likelihood-ratio test against the base is valid: a negative dAIC with a small
# p-value means the extra flexibility is worth its degrees of freedom. age^2 and
# the age spline are not nested in each other, so compare those two by AIC alone.

# comorbidity as nested factors: 0/1/2+ (base), 0/1/2/3+, then the fuller 0..5+.
# pmin() collapses the tail, so each coarser factor is a special case of the finer.
df$cci012p  <- factor(pmin(df$cci_n_conditions, 2))   # 0, 1, 2+  (base)
df$cci0123p <- factor(pmin(df$cci_n_conditions, 3))   # 0, 1, 2, 3+
df$ccicat   <- factor(pmin(df$cci_n_conditions, 5))   # 0, 1, 2, 3, 4, 5+

base    <- lm(wait ~ agediag + cci012p, df)                   # linear age, cci 0/1/2+
m_age2  <- lm(wait ~ agediag + I(agediag^2) + cci012p, df)    # base + quadratic age
m_agesp <- lm(wait ~ ns(agediag, 3) + cci012p, df)            # natural spline age (3 df)
m_cci3  <- lm(wait ~ agediag + cci0123p, df)                  # base, cci 0/1/2/3+
m_ccic  <- lm(wait ~ agediag + ccicat, df)                    # base, cci 0..5+

base_name <- "base: age + cci(0/1/2+)"
models <- list(base,      m_age2,     m_agesp,             m_cci3,          m_ccic)
labels <- c(base_name, "+ age^2", "age spline (3 df)", "cci 0/1/2/3+", "cci categorical (0..5+)")

# every candidate is tested against the base by a likelihood-ratio test. The test
# degrees of freedom are just the extra parameters, i.e. the drop in residual df.
base_aic   <- AIC(base)
ll_base    <- as.numeric(logLik(base))
dfres_base <- df.residual(base)

tab <- data.frame(model = character(), params = integer(), AIC = numeric(),
                  dAIC = numeric(), LR_chisq = numeric(), LR_df = integer(),
                  p_value = numeric(), stringsAsFactors = FALSE)
for (i in seq_along(models)) {
  m   <- models[[i]]
  aic <- AIC(m)
  if (labels[i] == base_name) {
    row <- data.frame(model = labels[i], params = length(coef(m)), AIC = aic, dAIC = 0,
                      LR_chisq = NA, LR_df = NA, p_value = NA, stringsAsFactors = FALSE)
  } else {
    stat <- 2 * (as.numeric(logLik(m)) - ll_base)
    dfd  <- dfres_base - df.residual(m)
    row  <- data.frame(model = labels[i], params = length(coef(m)), AIC = aic,
                       dAIC = aic - base_aic, LR_chisq = stat, LR_df = dfd,
                       p_value = pchisq(stat, dfd, lower.tail = FALSE),
                       stringsAsFactors = FALSE)
  }
  tab <- rbind(tab, row)
}

show <- tab
show$AIC      <- round(show$AIC, 1)
show$dAIC     <- round(show$dAIC, 1)
show$LR_chisq <- round(show$LR_chisq, 2)
show$p_value  <- signif(show$p_value, 3)
cat("\nOutcome-model coding (likelihood-ratio test vs the base model):\n")
print(show, row.names = FALSE)
cat("\nRead: dAIC < 0 with a small p-value means the richer coding fits wait on\n",
    "age + cci better. age^2 and the spline are each tested against linear age; the\n",
    "cci rows against cci 0/1/2+. A term that is significant here only matters for\n",
    "the mean-balancing weights if it also moves ranks - see the check below.\n")

# 3. do the pipeline codings move the hospital estimates / ranks? -------------
# standardise under each coding the pipeline can actually use (code_covariates
# supports linear or banded age, and cci as linear / 0-1-2+ / 0-1-2-3+) and
# compare the standardised means: rank Spearman against the continuous/continuous
# reference, and the mean effective sample size.
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
res <- data.frame(age = character(), cci = character(), rank_rho = numeric(),
                  mean_eff_n = numeric(), eff_n_q25 = numeric(), eff_n_q75 = numeric(),
                  eff_n_min = numeric(), eff_n_max = numeric(),
                  median_shift = numeric(), stringsAsFactors = FALSE)
for (i in seq_len(nrow(specs))) {
  s <- site_est(specs$age[i], specs$cci[i])
  m <- inner_join(ref, s, by = "hosp", suffix = c("_ref", ""))
  q <- quantile(s$eff_n, c(0.25, 0.75))
  res <- rbind(res, data.frame(
    age = specs$age[i], cci = specs$cci[i],
    rank_rho = cor(rank(m$stand_ref), rank(m$stand), method = "spearman"),
    mean_eff_n = mean(s$eff_n), eff_n_q25 = q[[1]], eff_n_q75 = q[[2]],
    eff_n_min = min(s$eff_n), eff_n_max = max(s$eff_n),
    median_shift = median(abs(m$stand - m$stand_ref)),
    stringsAsFactors = FALSE))
}
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

cat("\nRead: mean_eff_n can be pulled up by a few well-balanced hospitals, so the\n",
    "IQR (eff_n_q25 - eff_n_q75) and range (eff_n_min - eff_n_max) are given\n",
    "alongside it - a coding that looks fine on the mean but has a low eff_n_min\n",
    "is still eroding precision at specific, likely small, hospitals. Pick the\n",
    "simplest coding whose rank_rho vs the alternatives is high (ranks robust) and\n",
    "whose eff_n spread is not badly eroded. A low rho or a large median_shift means\n",
    "the coding choice matters and should follow the AIC / LRT evidence above.\n")