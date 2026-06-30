# Consistency of positive-deviant hospitals across comorbidity strata.
# Question: do the hospitals that look fast in the overall (whole-cohort)
# sustained analysis also look fast within each comorbidity stratum? The overall
# result is the reference on both comparisons; the two strata are not compared
# with each other. Three views are produced and the per-stratum caterpillars in
# 05 are left as they are.
#  1. "shake" caterpillars: each series placed in the overall ranking order, so
#     departures from a rising line show how much the ranking moves within strata
#  2. two concordance scatter plots, overall estimate (y) vs each stratum (x)
#  3. a concordance table and rank correlations

library(rstan)
library(dplyr)
library(ggplot2)

source("R/01_config.R")
prob_cut <- 0.80   # candidate threshold, matching 05

process_fit <- function(obj) {
  draws <- rstan::extract(obj$fit, pars = "y_site_true")$y_site_true
  obj$site %>%
    mutate(post_mean = colMeans(draws),
           post_sd   = apply(draws, 2, sd),
           ci_lo = post_mean - 1.96 * post_sd,
           ci_hi = post_mean + 1.96 * post_sd) %>%
    bind_cols(rank_metrics(draws))
}

overall <- process_fit(readRDS(file.path(out_dir, "stan_sustained.rds")))
s01     <- process_fit(readRDS(file.path(out_dir, "stan_strata_01.rds")))
s2      <- process_fit(readRDS(file.path(out_dir, "stan_strata_2.rds")))

# overall positive-deviant set and the common ordering (fastest first)
cand_set <- overall %>% filter(p_top20 >= prob_cut) %>% pull(diag_hosp)
ord <- overall %>% arrange(post_mean) %>%
  transmute(diag_hosp, rank_order = row_number())

# 1. shake caterpillars ------------------------------------------------------
cat_long <- bind_rows(
  overall %>% transmute(diag_hosp, series = "Overall", post_mean, ci_lo, ci_hi),
  s01     %>% transmute(diag_hosp, series = "CCI 0-1", post_mean, ci_lo, ci_hi),
  s2      %>% transmute(diag_hosp, series = "CCI 2+",  post_mean, ci_lo, ci_hi)
) %>%
  inner_join(ord, by = "diag_hosp") %>%
  mutate(series = factor(series, levels = c("Overall", "CCI 0-1", "CCI 2+")),
         cand   = diag_hosp %in% cand_set)

p_cat <- ggplot(cat_long, aes(rank_order, post_mean)) +
  geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), colour = "grey78") +
  geom_point(aes(colour = cand), size = 0.8) +
  scale_colour_manual(values = c("grey25", "firebrick"),
                      labels = c("other", "overall top 20%"), name = NULL) +
  facet_grid(series ~ ., scales = "free_y") +
  labs(x = "Hospital (ordered by overall sustained estimate, fastest first)",
       y = "Standardised days to decision-to-treat",
       title = "Stratum estimates in the overall ranking order") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom")
ggsave(file.path(out_dir, "strata_shake_caterpillar.pdf"), p_cat, width = 7, height = 8)

# 2. concordance scatter plots -----------------------------------------------
base <- overall %>%
  transmute(diag_hosp, overall = post_mean, o_lo = ci_lo, o_hi = ci_hi, n,
            cand = diag_hosp %in% cand_set)

conc_plot <- function(strat_df, xlab, title) {
  d <- base %>%
    inner_join(strat_df %>% transmute(diag_hosp, strat = post_mean,
                                      s_lo = ci_lo, s_hi = ci_hi),
               by = "diag_hosp")
  r   <- cor(d$overall, d$strat)
  rho <- cor(d$overall, d$strat, method = "spearman")
  ggplot(d, aes(strat, overall)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_errorbar(aes(ymin = o_lo, ymax = o_hi), colour = "grey85", width = 0) +
    geom_errorbarh(aes(xmin = s_lo, xmax = s_hi), colour = "grey85", height = 0) +
    geom_point(aes(size = n, colour = cand), alpha = 0.6) +
    scale_colour_manual(values = c("grey25", "firebrick"),
                        labels = c("other", "overall top 20%"), name = NULL) +
    scale_size_area(guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f", r, rho)) +
    labs(x = xlab, y = "Overall sustained estimate (days)", title = title) +
    theme_bw() + theme(legend.position = "bottom")
}

ggsave(file.path(out_dir, "concordance_overall_vs_cci01.pdf"),
       conc_plot(s01, "CCI 0-1 sustained estimate (days)", "(a) Overall vs CCI 0-1"),
       width = 6, height = 6)
ggsave(file.path(out_dir, "concordance_overall_vs_cci2.pdf"),
       conc_plot(s2, "CCI 2+ sustained estimate (days)", "(b) Overall vs CCI 2+"),
       width = 6, height = 6)

# 3. concordance table and correlations --------------------------------------
conc_tab <- overall %>%
  transmute(diag_hosp, n,
            overall_est = post_mean, overall_ptop20 = p_top20,
            overall_rank = exp_rank, overall_cand = diag_hosp %in% cand_set) %>%
  left_join(s01 %>% transmute(diag_hosp, cci01_est = post_mean,
                              cci01_ptop20 = p_top20, cci01_rank = exp_rank),
            by = "diag_hosp") %>%
  left_join(s2 %>% transmute(diag_hosp, cci2_est = post_mean,
                             cci2_ptop20 = p_top20, cci2_rank = exp_rank),
            by = "diag_hosp") %>%
  arrange(desc(overall_ptop20))
write.csv(conc_tab, file.path(out_dir, "strata_concordance_table.csv"), row.names = FALSE)

sp01 <- conc_tab %>% filter(!is.na(cci01_rank))
sp2  <- conc_tab %>% filter(!is.na(cci2_rank))
cat(sprintf("overall vs CCI 0-1: Spearman %.3f (n = %d)\n",
            cor(sp01$overall_rank, sp01$cci01_rank, method = "spearman"), nrow(sp01)))
cat(sprintf("overall vs CCI 2+ : Spearman %.3f (n = %d)\n",
            cor(sp2$overall_rank,  sp2$cci2_rank,  method = "spearman"), nrow(sp2)))

# of the overall candidates, how many also reach the top 20% within each stratum
cat(sprintf("overall candidates: %d\n", length(cand_set)))
cat(sprintf("  also top 20%% in CCI 0-1: %d\n",
            sum(conc_tab$overall_cand & conc_tab$cci01_ptop20 >= prob_cut, na.rm = TRUE)))
cat(sprintf("  also top 20%% in CCI 2+ : %d\n",
            sum(conc_tab$overall_cand & conc_tab$cci2_ptop20 >= prob_cut, na.rm = TRUE)))

