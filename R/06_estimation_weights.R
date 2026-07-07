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
library(ggrepel)
library(patchwork)

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
  raw_n <- as.numeric(table(Z))
  
  # one pass per lambda, keeping the whole per-hospital effective-n vector (not
  # just its mean) so the spread across hospitals can be shown as well.
  per <- lapply(grid, function(l) {
    so <- standardize(X, rep(0, ncol(X)), Z, lambda = l, exact_global = FALSE)
    w  <- extract_weights(so)
    wm <- (t(so$weights) %*% X)                        # hospital x covariate
    wt <- as.numeric(abs(wm %*% beta))
    ne <- tapply(w, Z, function(x) sum(x)^2 / sum(x^2))
    list(summary = c(lambda       = l,
                     bias_removed = 1 - mean(wt) / mean(unw),
                     mean_eff_n   = mean(ne),
                     mean_deff    = mean(ne / raw_n)),
         ess = data.frame(lambda = l, hosp = names(ne),
                          ess = as.numeric(ne), stringsAsFactors = FALSE))
  })
  list(summary = as.data.frame(do.call(rbind, lapply(per, `[[`, "summary"))),
       ess     = do.call(rbind, lapply(per, `[[`, "ess")),
       raw     = data.frame(hosp = names(table(Z)), n = raw_n, stringsAsFactors = FALSE))
}

# main analysis: primary standardisation (age + cci + season + calendar year) -
# season and calendar year are part of the primary adjustment set; their
# inclusion, and the choice to keep age linear, are justified against outcome
# fit and effective sample size in script 15. Other patient factors (sex,
# ethnicity, deprivation, stage) remain excluded.
cv          <- code_covariates(df)
primary_bin <- c(cv$bin, season_terms, year_term)
to          <- balance_tradeoff(cv$data, cv$cont, primary_bin)
trade       <- to$summary
ess_by_hosp <- to$ess
raw_by_hosp <- to$raw
cat("balance vs effective n by lambda:\n"); print(round(trade, 3))
write.csv(trade,       file.path(out_dir, "lambda_tradeoff.csv"),        row.names = FALSE)
write.csv(ess_by_hosp, file.path(out_dir, "lambda_ess_by_hospital.csv"), row.names = FALSE)

# lambda trade-off figure. The weights shrink toward uniform as lambda grows: a
# larger lambda keeps more effective sample size but removes less case-mix
# imbalance. One point per lambda (small lambda top-left, most bias removed and
# fewest effective patients; larger lambda moves down and right), the working
# value highlighted. Same look as the other figures (theme_classic, darkblue).
axis_title_size <- 13
axis_text_size  <- 12
label_size      <- 4.2

trade_curve <- trade %>%
  mutate(pct_bias_removed = 100 * bias_removed,
         is_main = abs(lambda - lambda_main) < 1e-9)
lab_near <- trade_curve[trade_curve$lambda <  2, ]   # already well spaced: plain labels
lab_far  <- trade_curve[trade_curve$lambda >= 2, ]   # 2, 2.5, 3 bunch up: repel just these

p_trade <- ggplot(trade_curve, aes(mean_eff_n, pct_bias_removed)) +
  geom_path(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = is_main, size = is_main)) +
  geom_text(data = lab_near, aes(label = lambda), vjust = -0.9, size = label_size) +
  geom_text_repel(data = lab_far, aes(label = lambda), size = label_size,
                  direction = "y", nudge_y = 6, box.padding = 0.4,
                  min.segment.length = 0, segment.colour = "grey70",
                  segment.size = 0.3, seed = 1) +
  scale_colour_manual(values = c(`FALSE` = "darkblue", `TRUE` = "firebrick"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.8, `TRUE` = 3.8), guide = "none") +
  scale_x_continuous("Average effective sample size per hospital") +
  scale_y_continuous("Average percentage bias reduction (%)", breaks = seq(0, 100, 20)) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme_classic(base_size = 13) +
  theme(axis.title      = element_text(size = axis_title_size),
        axis.text       = element_text(size = axis_text_size),
        legend.position = "none")

ggsave(file.path(out_dir, "lambda_tradeoff.png"), p_trade,
       width = 120, height = 95, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "lambda_tradeoff.pdf"), p_trade,
       width = 120, height = 95, units = "mm", bg = "white")
cat("lambda trade-off figure written to lambda_tradeoff.png / .pdf\n")

# (b) spread of effective sample size across hospitals at each lambda, with the
# raw per-hospital counts as a reference. Each violin is the distribution of
# per-hospital effective n at one lambda (quartile lines inside); the grey box on
# the right is the raw count the weights shrink from. Effective n sits below the
# raw counts and climbs toward them as lambda grows. Same look as panel (a).
lam_lev <- as.character(lambda_grid)
ess_by_hosp$x <- factor(as.character(ess_by_hosp$lambda), levels = c(lam_lev, "Raw n"))
raw_box <- data.frame(x = factor("Raw n", levels = c(lam_lev, "Raw n")),
                      n = raw_by_hosp$n)

p_ess <- ggplot() +
  geom_hline(yintercept = median(raw_by_hosp$n), linetype = "dashed",
             colour = "grey75", linewidth = 0.4) +
  geom_violin(data = ess_by_hosp, aes(x, ess), fill = "darkblue", colour = "darkblue",
              alpha = 0.25, linewidth = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(data = raw_box, aes(x, n), width = 0.5, fill = "grey80",
               colour = "grey30", linewidth = 0.4, outlier.size = 0.6) +
  geom_vline(xintercept = length(lam_lev) + 0.5, linetype = "dotted", colour = "grey70") +
  scale_x_discrete("Balancing weight regularisation, lambda") +
  scale_y_continuous("Sample size per hospital",
                     expand = expansion(mult = c(0.02, 0.05))) +
  theme_classic(base_size = 13) +
  theme(axis.title  = element_text(size = axis_title_size),
        axis.text   = element_text(size = axis_text_size),
        axis.text.x = element_text(size = axis_text_size - 1),
        legend.position = "none")

combined <- p_trade / p_ess +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
ggsave(file.path(out_dir, "lambda_tradeoff_ess.png"), combined,
       width = 150, height = 200, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "lambda_tradeoff_ess.pdf"), combined,
       width = 150, height = 200, units = "mm", bg = "white")
cat("combined trade-off + effective-n figure written to lambda_tradeoff_ess.png / .pdf\n")

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