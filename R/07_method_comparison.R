# 07  ranking-method comparison and funnel plots
# -----------------------------------------------------------------------------
# Compares hospital rankings under different case-mix-adjustment methods on the
# sustained estimand (the basis for Table 3), and draws funnel plots. Methods:
#   - raw mean and median (no adjustment)
#   - balancing-weights direct standardisation, mean and median, on the age +
#     comorbidity set and on the full patient mix
#   - regression direct standardisation (fixed hospital effects, g-computation)
#   - regression direct standardisation with shrinkage (random hospital effects
#     via a linear mixed model, g-computed over the patient mix) - the regression
#     analogue of the headline balancing-weights-plus-Bayesian-shrinkage method
#   - indirect standardisation (observed minus expected)
#   - the shrunk expected rank from the Bayesian fit (the headline method)

library(dplyr)
library(ggplot2)
library(lme4)

source("R/01_config.R")
df    <- readRDS(file.path(out_dir, "analysis_data.rds"))
fit_p <- readRDS(file.path(out_dir, "fit_primary.rds"))
fit_f <- readRDS(file.path(out_dir, "fit_full.rds"))
sus   <- read.csv(file.path(out_dir, "ranks_sustained.csv"))

# regression direct standardisation (fixed hospital effects): fit hospital as a
# fixed effect, predict every hospital over the whole patient mix, and average.
reg_std <- function(d, cont, bin) {
  d <- d %>% mutate(hospf = factor(hosp))
  form <- as.formula(paste("wait ~ hospf +", paste(c(cont, bin), collapse = " + ")))
  m <- lm(form, data = d)
  hs <- levels(d$hospf)
  vapply(hs, function(h) { nd <- d; nd$hospf <- factor(h, levels = hs); mean(predict(m, nd)) },
         numeric(1))
}
g_primary <- reg_std(df, cont_vars, bin_primary)
g_full    <- reg_std(df, cont_vars, bin_full)

# regression direct standardisation with shrinkage: hospital as a random effect.
# The g-computation average over the whole patient mix equals the population
# fixed-effect prediction plus the hospital's shrunken effect (its BLUP), so the
# hospital estimates are pulled toward the average exactly as in the Bayesian
# headline method, but fitted by a mixed model rather than balancing weights.
reg_std_shrunk <- function(d, cont, bin) {
  d <- d %>% mutate(hospf = factor(hosp))
  form <- as.formula(paste("wait ~", paste(c(cont, bin), collapse = " + "), "+ (1 | hospf)"))
  m <- lmer(form, data = d, REML = TRUE)
  pop_fixed <- mean(predict(m, re.form = NA))           # population fixed-effect mean
  re <- ranef(m)$hospf[, "(Intercept)"]                 # shrunken hospital effects
  names(re) <- rownames(ranef(m)$hospf)
  hs <- levels(d$hospf)
  setNames(pop_fixed + re[hs], hs)
}
g_shrunk_primary <- reg_std_shrunk(df, cont_vars, bin_primary)
g_shrunk_full    <- reg_std_shrunk(df, cont_vars, bin_full)

# indirect standardisation: expected wait from a pooled model on each hospital's
# own patients, then observed - expected + grand mean.
pm <- lm(as.formula(paste("wait ~", paste(c(cont_vars, bin_primary), collapse = " + "))),
         data = df)
df$pred <- predict(pm)
grand   <- mean(df$wait)
ind <- df %>% group_by(hosp) %>%
  summarise(n = n(), obs = mean(wait), exp = mean(pred), .groups = "drop") %>%
  mutate(indirect = obs - exp + grand)

# assemble point estimates per hospital --------------------------------------
comp <- fit_p$site %>%
  transmute(hosp, diag_hosp,
            raw_mean, raw_median,
            wt_ac_mean = stand, wt_ac_med = stand_med) %>%
  left_join(fit_f$site %>% transmute(hosp, wt_full_mean = stand, wt_full_med = stand_med),
            by = "hosp") %>%
  left_join(tibble(hosp = as.integer(names(g_primary)),        reg_ac        = g_primary),        by = "hosp") %>%
  left_join(tibble(hosp = as.integer(names(g_full)),           reg_full      = g_full),           by = "hosp") %>%
  left_join(tibble(hosp = as.integer(names(g_shrunk_primary)), reg_shrunk_ac = g_shrunk_primary), by = "hosp") %>%
  left_join(tibble(hosp = as.integer(names(g_shrunk_full)),    reg_shrunk_full = g_shrunk_full),  by = "hosp") %>%
  left_join(ind %>% select(hosp, indirect), by = "hosp") %>%
  left_join(sus %>% select(hosp, shrunk_rank = exp_rank), by = "hosp")

# turn each estimate into a rank (1 = fastest); shrunk_rank is already a rank
methods <- c("raw_mean","raw_median","wt_ac_mean","wt_ac_med",
             "wt_full_mean","wt_full_med","reg_ac","reg_full",
             "reg_shrunk_ac","reg_shrunk_full","indirect")
ranks <- comp %>%
  mutate(across(all_of(methods), ~ rank(.x, ties.method = "average"),
                .names = "rank_{.col}")) %>%
  mutate(rank_shrunk = rank(shrunk_rank, ties.method = "average"))

rank_cols <- c(paste0("rank_", methods), "rank_shrunk")
rho <- cor(ranks[rank_cols], method = "spearman")
cat("Spearman rank correlation between methods:\n"); print(round(rho, 2))

write.csv(comp, file.path(out_dir, "method_comparison_estimates.csv"), row.names = FALSE)
write.csv(ranks %>% select(hosp, diag_hosp, all_of(rank_cols)),
          file.path(out_dir, "method_comparison_ranks.csv"), row.names = FALSE)
write.csv(round(rho, 3), file.path(out_dir, "method_comparison_rho.csv"))

# funnel plot, sustained -----------------------------------------------------
sigma_r <- sd(resid(pm))
fp <- ind %>% mutate(
  lo95 = grand - 1.96  * sigma_r / sqrt(n),
  hi95 = grand + 1.96  * sigma_r / sqrt(n),
  lo998 = grand - 3.09 * sigma_r / sqrt(n),
  hi998 = grand + 3.09 * sigma_r / sqrt(n),
  fast = indirect < lo95)
ggsave(file.path(out_dir, "funnel_sustained.pdf"),
       ggplot(fp, aes(n, indirect)) +
         geom_hline(yintercept = grand, colour = "grey50") +
         geom_line(aes(y = lo95),  linetype = 2, colour = "grey40") +
         geom_line(aes(y = hi95),  linetype = 2, colour = "grey40") +
         geom_line(aes(y = lo998), linetype = 3, colour = "grey60") +
         geom_line(aes(y = hi998), linetype = 3, colour = "grey60") +
         geom_point(aes(colour = fast), size = 1) +
         scale_colour_manual(values = c("black", "firebrick"), guide = "none") +
         labs(x = "Hospital volume", y = "Indirectly standardised days to DTT",
              title = "Funnel plot, sustained") +
         theme_bw(),
       width = 7, height = 5)

# funnel plot, improvement ---------------------------------------------------
imp <- readRDS(file.path(out_dir, "site_improve.rds")) %>%
  mutate(n = n1 + n2)
sd_delta <- median(imp$se_delta * sqrt(imp$n))
center   <- 0
fimp <- imp %>% mutate(
  lo95 = center - 1.96 * sd_delta / sqrt(n),
  hi95 = center + 1.96 * sd_delta / sqrt(n),
  fast = delta < lo95)
ggsave(file.path(out_dir, "funnel_improve.pdf"),
       ggplot(fimp, aes(n, delta)) +
         geom_hline(yintercept = 0, colour = "grey50") +
         geom_line(aes(y = lo95), linetype = 2, colour = "grey40") +
         geom_line(aes(y = hi95), linetype = 2, colour = "grey40") +
         geom_point(aes(colour = fast), size = 1) +
         scale_colour_manual(values = c("black", "firebrick"), guide = "none") +
         labs(x = "Hospital volume (both halves)",
              y = "Change in standardised days", title = "Funnel plot, improvement") +
         theme_bw(),
       width = 7, height = 5)
