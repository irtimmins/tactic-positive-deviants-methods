# Posterior processing: shrunk hospital means with credible intervals,
# expected ranks, probability of being in the top X%, caterpillar plots, and
# selection of sustained / improved positive-deviant candidates. The rank logic
# follows the template (rank within each posterior draw), with shorter waits as
# better so rank 1 is the best performer.

library(rstan)
library(dplyr)
library(ggplot2)

source("00_config_funcs.R")

# optional: a diag_hosp code to highlight in the caterpillar plots
highlight_hosp <- NA

process_fit <- function(obj, id_col = "diag_hosp") {
  draws <- rstan::extract(obj$fit, pars = "y_site_true")$y_site_true  # draws x J
  site  <- obj$site
  rm <- rank_metrics(draws)
  site %>%
    mutate(
      post_mean = colMeans(draws),
      post_sd   = apply(draws, 2, sd),
      ci_lo     = post_mean - 1.96 * post_sd,
      ci_hi     = post_mean + 1.96 * post_sd
    ) %>%
    bind_cols(rm)
}

caterpillar <- function(d, ylab, title, highlight = NA) {
  d <- d %>% arrange(post_mean) %>% mutate(rank_order = row_number())
  hl <- if (!is.na(highlight)) d$diag_hosp == highlight else rep(FALSE, nrow(d))
  ggplot(d, aes(rank_order, post_mean)) +
    geom_hline(yintercept = mean(d$post_mean), colour = "grey50") +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), colour = "grey70") +
    geom_point(size = 0.9) +
    { if (any(hl)) geom_point(data = d[hl, ], colour = "firebrick", size = 2) } +
    labs(x = "Hospital (ordered)", y = ylab, title = title) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# sustained ------------------------------------------------------------------
sus <- process_fit(readRDS(file.path(out_dir, "stan_sustained.rds")))
write.csv(sus, file.path(out_dir, "ranks_sustained.csv"), row.names = FALSE)

ggsave(file.path(out_dir, "caterpillar_sustained.pdf"),
       caterpillar(sus, "Standardised days to decision-to-treat",
                   "Sustained performance", highlight_hosp),
       width = 7, height = 6)

# improvement ----------------------------------------------------------------
imp <- process_fit(readRDS(file.path(out_dir, "stan_improve.rds")))
write.csv(imp, file.path(out_dir, "ranks_improve.csv"), row.names = FALSE)

ggsave(file.path(out_dir, "caterpillar_improve.pdf"),
       caterpillar(imp, "Change in standardised days (second half - first)",
                   "Improvement", highlight_hosp),
       width = 7, height = 6)

# strata caterpillars --------------------------------------------------------
for (tag in c("01", "2")) {
  f <- file.path(out_dir, sprintf("stan_strata_%s.rds", tag))
  if (!file.exists(f)) next
  st <- process_fit(readRDS(f))
  lab <- ifelse(tag == "01", "CCI 0-1", "CCI 2+")
  write.csv(st, file.path(out_dir, sprintf("ranks_strata_%s.csv", tag)), row.names = FALSE)
  ggsave(file.path(out_dir, sprintf("caterpillar_strata_%s.pdf", tag)),
         caterpillar(st, "Standardised days to decision-to-treat",
                     paste("Sustained performance,", lab), highlight_hosp),
         width = 7, height = 6)
}

# candidate selection --------------------------------------------------------
# a sustained candidate has a high posterior probability of sitting in the
# fastest 20%; an improved candidate likewise on the change estimate.
prob_cut <- 0.80
sus_cand <- sus %>% filter(p_top20 >= prob_cut) %>% pull(diag_hosp)
imp_cand <- imp %>% filter(p_top20 >= prob_cut) %>% pull(diag_hosp)

cat(sprintf("sustained candidates (P(top 20%%) >= %.2f): %d\n", prob_cut, length(sus_cand)))
cat(sprintf("improved candidates  (P(top 20%%) >= %.2f): %d\n", prob_cut, length(imp_cand)))
cat("in both:", length(intersect(sus_cand, imp_cand)), "\n")

candidates <- tibble(
  diag_hosp = union(sus_cand, imp_cand)) %>%
  mutate(sustained = diag_hosp %in% sus_cand,
         improved  = diag_hosp %in% imp_cand,
         both      = sustained & improved)
write.csv(candidates, file.path(out_dir, "candidates.csv"), row.names = FALSE)

# overlap of the two candidate sets (simple two-set diagram)
n_s <- length(setdiff(sus_cand, imp_cand))
n_i <- length(setdiff(imp_cand, sus_cand))
n_b <- length(intersect(sus_cand, imp_cand))
pdf(file.path(out_dir, "candidate_overlap.pdf"), width = 6, height = 4)
plot.new(); plot.window(c(0, 10), c(0, 6), asp = 1)
symbols(c(4, 6), c(3, 3), circles = c(2, 2), inches = FALSE, add = TRUE)
text(3, 3, n_s); text(7, 3, n_i); text(5, 3, n_b)
text(3, 5.4, "Sustained"); text(7, 5.4, "Improved")
dev.off()