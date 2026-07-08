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

# main analysis fit at the working lambda; done before the figures because panel
# (b) below shows this fit's per-hospital effective sample size.
fit_main <- run_standardise(patient_data          = cv$data,
                            continuous_covariates = cv$cont,
                            binary_covariates     = primary_bin,
                            lambda                = lambda_main)
site_sustained <- fit_main$site
saveRDS(fit_main,       file.path(out_dir, "fit_primary.rds"))
saveRDS(site_sustained, file.path(out_dir, "site_sustained.rds"))

# lambda trade-off figure. The weights shrink toward uniform as lambda grows: a
# larger lambda keeps more effective sample size but removes less case-mix
# imbalance. One point per lambda (small lambda top-left, most bias removed and
# fewest effective patients; larger lambda moves down and right), the working
# value highlighted. Same look as the other figures (theme_classic, darkblue).
axis_title_size <- 13
axis_text_size  <- 12
label_size      <- 4.2
tag_size        <- 15
tag_x           <- 0      # a/b label horizontal position within each panel (0 = left)
tag_y           <- 1.03   # a/b label vertical position; raise above 1 to nudge higher
repel_spread_x  <- 2.5    # sideways fan for the crowded lambda labels (raise for more left/right)
repel_lift_y    <- 4      # small upward lift so those labels clear their points
near_nudge_x    <- 0.8    # small-lambda labels: slight rightward shift
near_nudge_y    <- 6      # small-lambda labels: slight upward shift
zero_lift_y     <- -1      # extra upward lift for the lambda = 0 label only

# plotting window. anything with effective n outside this range is cropped rather
# than trailing off the axis into the margin (clip is off for the labels and tag).
x_lo <- 60
x_hi <- 77
trade_curve <- trade %>%
  mutate(pct_bias_removed = 100 * bias_removed,
         is_main = abs(lambda - lambda_main) < 1e-9)
vis      <- trade_curve[trade_curve$mean_eff_n >= x_lo & trade_curve$mean_eff_n <= x_hi, ]
lab_near <- vis[vis$lambda <  1, ]   # well spaced along the curve
lab_far  <- vis[vis$lambda >= 1, ]   # 1 to 3 bunch up: let ggrepel place these
# fan the crowded labels out sideways from their centre, so repulsion resolves
# mostly left/right rather than stacking them vertically.
lab_far$nudge_lr <- repel_spread_x * (lab_far$mean_eff_n - mean(lab_far$mean_eff_n))
# lambda = 0 sits almost on top of 0.01, so repulsion keeps levelling the two; place
# its label at a fixed height above instead (raise zero_lift_y to lift it further).
lab_zero <- lab_near[lab_near$lambda == 0, ]
lab_near <- lab_near[lab_near$lambda != 0, ]

p_trade <- ggplot(vis, aes(mean_eff_n, pct_bias_removed)) +
  geom_path(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(colour = is_main, size = is_main)) +
  geom_text_repel(data = lab_near, aes(label = lambda), size = label_size,
                  nudge_x = near_nudge_x, nudge_y = near_nudge_y,
                  direction = "both", force = 0.5, box.padding = 0.3,
                  point.padding = 0.2, max.overlaps = Inf,
                  segment.colour = "grey70", segment.size = 0.3, seed = 1) +
  geom_text(data = lab_zero, aes(label = lambda), size = label_size,
            nudge_x = near_nudge_x, nudge_y = near_nudge_y + zero_lift_y) +
  geom_text_repel(data = lab_far, aes(label = lambda), size = label_size,
                  nudge_x = lab_far$nudge_lr, nudge_y = repel_lift_y,
                  direction = "both", box.padding = 0.6, point.padding = 0.4,
                  force = 2, max.overlaps = Inf, min.segment.length = 0,
                  segment.colour = "grey70", segment.size = 0.3, seed = 1) +
  scale_colour_manual(values = c(`FALSE` = "darkblue", `TRUE` = "firebrick"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.8, `TRUE` = 3.8), guide = "none") +
  scale_x_continuous("Average effective sample size per hospital",
                     breaks = c(60, 65, 70, 75)) +
  scale_y_continuous("Average percentage bias reduction (%)", breaks = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0, 100), clip = "off") +
  theme_classic(base_size = 13) +
  theme(axis.title      = element_text(size = axis_title_size),
        axis.text       = element_text(size = axis_text_size),
        legend.position = "none")

ggsave(file.path(out_dir, "lambda_tradeoff.png"), p_trade,
       width = 120, height = 95, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "lambda_tradeoff.pdf"), p_trade,
       width = 120, height = 95, units = "mm", bg = "white")
cat("lambda trade-off figure written to lambda_tradeoff.png / .pdf\n")

# (b) overlap check at the working lambda: each hospital's effective sample size
# after weighting against its raw sample size, point area proportional to the raw
# count. Effective n can only sit at or below the raw count (dashed y = x line),
# and falls further below it where a hospital's case-mix is further from the
# reference. Same look as panel (a).
nmax        <- max(site_sustained$n)
size_breaks <- c(25, 50, 100, 200)

p_scatter <- ggplot(site_sustained, aes(n, n_eff)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.4) +
  geom_point(aes(size = n), colour = "darkblue", alpha = 0.5) +
  scale_size_area("Hospital sample size", max_size = 6, breaks = size_breaks) +
  scale_x_continuous("Hospital sample size", limits = c(0, nmax),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_y_continuous("Hospital effective sample size", limits = c(0, nmax),
                     expand = expansion(mult = c(0, 0.04))) +
  coord_fixed(ratio = 1) +                       # equal units both axes: y = x is a true 45 degrees
  theme_classic(base_size = 13) +
  theme(axis.title           = element_text(size = axis_title_size),
        axis.text            = element_text(size = axis_text_size),
        legend.position      = c(0.03, 0.97),    # inside the empty upper-left corner
        legend.justification = c(0, 1),
        legend.background    = element_rect(fill = "white", colour = NA),
        legend.key           = element_blank(),
        legend.title         = element_text(size = axis_text_size - 1),
        legend.text          = element_text(size = axis_text_size - 1))

combined <- p_trade / p_scatter +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag          = element_text(size = tag_size, face = "bold"),
        plot.tag.position = c(tag_x, tag_y),
        plot.margin       = margin(t = 12, r = 6, b = 6, l = 6))
ggsave(file.path(out_dir, "lambda_tradeoff_ess.png"), combined,
       width = 125, height = 220, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "lambda_tradeoff_ess.pdf"), combined,
       width = 125, height = 220, units = "mm", bg = "white")
cat("combined trade-off + effective-n figure written to lambda_tradeoff_ess.png / .pdf\n")

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