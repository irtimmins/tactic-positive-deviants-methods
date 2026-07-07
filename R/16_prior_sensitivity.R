# 16  prior sensitivity of the shrinkage step
# -----------------------------------------------------------------------------
# The headline ranking (07) shrinks each hospital's estimate with a normal-normal
# model: the hospital true means are y_site_true ~ Normal(mu, tau), the observed
# estimates are y_obs ~ Normal(y_site_true, se), and the between-hospital sd has
# the prior tau ~ half-Cauchy(0, prior_tau_scale) (= half-Cauchy(0, 10) days).
#
# This script re-fits the same model under alternative priors on tau and checks
# how far the ranking moves. ONLY the tau prior changes. The prior on the grand
# mean is UNCHANGED throughout: mu ~ Normal(prior_mu_mean, prior_mu_sd), with
# prior_mu_sd = 50 days and prior_mu_mean the data mean for the sustained
# estimand and 0 for the change score (exactly as in the main fits). The methods:
#   main   tau ~ half-Cauchy(0, 10)          the headline prior (reference column)
#   1      tau ~ half-Normal(0, 25)          vague, light tails
#   2      tau ~ Uniform(0, 100)             flat over a wide range
#   3      tau ~ half-Student-t(3, 0, 10)    heavier tails than normal, lighter
#                                            than the Cauchy
#   DL     empirical Bayes: tau^2 by the DerSimonian-Laird method of moments, then
#          the analytic shrinkage mu + B (y - mu), B = tau^2 / (tau^2 + se^2)
# Every method is a shrinkage method; only the prior on tau (or, for DL, the tau
# estimate) differs. Each column title spells its tau prior out in full.
#
# For the sustained and improvement estimands it builds a rank-movement table in
# the style of movement_methods.docx: the main-model rank, then the rank change
# under each of the four sensitivity methods (green up / red down).

library(rstan)
library(dplyr)
library(flextable)
library(officer)

source("R/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

priors_stan <- file.path(stan_dir, "dp_normal_priors.stan")

# the four sensitivity priors on the between-hospital sd (tau). The main model
# uses tau ~ half-Cauchy(0, prior_tau_scale); the alternatives are deliberately
# different vague choices plus a heavier-tailed Student-t. The codes match
# dp_normal_priors.stan (1 half-normal, 2 uniform, 3 half-student_t). The uniform
# needs a finite upper bound; the others get a large one so it never binds.
big_upper <- 500
prior_specs <- list(
  list(label = "Vague half-normal", code = 1, scale = 25, df = 3, upper = big_upper),
  list(label = "Vague uniform",     code = 2, scale = 10, df = 3, upper = 100),
  list(label = "Student-t",         code = 3, scale = 10, df = 3, upper = big_upper))

# write each prior's tau distribution out in full, from its own parameters so the
# printed title cannot drift from what is fitted. These become the column titles.
tau_string <- function(spec) switch(spec$code,
                                    sprintf("tau ~ half-Normal(0, %g)", spec$scale),                    # 1
                                    sprintf("tau ~ Uniform(0, %g)", spec$upper),                        # 2
                                    sprintf("tau ~ half-Student-t(%g, 0, %g)", spec$df, spec$scale))    # 3

main_title  <- sprintf("Main model: tau ~ half-Cauchy(0, %g)", prior_tau_scale)
sens_titles <- vapply(prior_specs,
                      function(s) sprintf("%s: %s", s$label, tau_string(s)), character(1))
dl_title    <- "DerSimonian-Laird: empirical Bayes, plug-in tau"
sens_cols   <- c(sens_titles, dl_title)

# fit one Bayesian sensitivity prior; return the expected posterior rank per
# hospital, in the input order.
fit_prior <- function(y, se, spec, mu_mean) {
  dat <- list(J = length(y), y_site_obs = y, sigma_site_obs = se,
              prior_mu_mean = mu_mean, prior_mu_sd = prior_mu_sd,
              tau_prior = spec$code, tau_scale = spec$scale,
              tau_upper = spec$upper, tau_df = spec$df)
  fit <- rstan::stan(priors_stan, data = dat, seed = 8675309,
                     chains = 4, iter = 4000, warmup = 2000, refresh = 0,
                     control = list(adapt_delta = 0.95, max_treedepth = 12))
  rank_metrics(rstan::extract(fit, pars = "y_site_true")$y_site_true)$exp_rank
}

# empirical-Bayes DerSimonian-Laird shrinkage: between-hospital variance tau^2 by
# the method of moments about the precision-weighted mean (held at zero below its
# null expectation), giving the analytic posterior mean mu + B (y - mu) and its
# posterior sd sqrt(B se^2), with B = tau^2 / (tau^2 + se^2). The expected rank is
# then read off simulated draws from that Gaussian empirical-Bayes posterior, the
# same way the Bayesian methods' expected ranks are formed.
dl_rank <- function(y, se, n_draws = 4000) {
  w     <- 1 / se^2
  mu    <- sum(w * y) / sum(w)
  Q     <- sum(w * (y - mu)^2)
  denom <- sum(w) - sum(w^2) / sum(w)
  tau2  <- max((Q - (length(y) - 1)) / denom, 0)
  B     <- tau2 / (tau2 + se^2)
  shrunk_mean <- mu + B * (y - mu)
  shrunk_sd   <- sqrt(B * se^2)
  set.seed(8675309)
  draws <- matrix(rnorm(n_draws * length(y),
                        mean = rep(shrunk_mean, each = n_draws),
                        sd   = rep(shrunk_sd,   each = n_draws)),
                  nrow = n_draws)
  rank_metrics(draws)$exp_rank
}

# per-hospital expected posterior rank for one estimand: the main model (read from
# the saved 07 fit) plus the four sensitivity methods, all in hospital order.
sensitivity_estimates <- function(site_file, main_fit_file, y_col, se_col, mu0 = NULL) {
  site <- readRDS(file.path(out_dir, site_file)) %>% arrange(hosp)
  y  <- site[[y_col]]; se <- site[[se_col]]
  mu_mean <- if (is.null(mu0)) mean(y) else mu0
  
  main <- readRDS(file.path(out_dir, main_fit_file))
  est  <- data.frame(hosp = site$hosp, diag_hosp = site$diag_hosp,
                     stringsAsFactors = FALSE)
  est[[main_title]] <- rank_metrics(rstan::extract(main$fit, pars = "y_site_true")$y_site_true)$exp_rank
  for (i in seq_along(prior_specs))
    est[[sens_titles[i]]] <- fit_prior(y, se, prior_specs[[i]], mu_mean)
  est[[dl_title]] <- dl_rank(y, se)
  est
}

sus_est <- sensitivity_estimates("site_sustained.rds", "stan_sustained.rds",
                                 "stand_adj", "se_adj_pool")            # shrink to the mean
imp_est <- sensitivity_estimates("site_improve.rds",   "stan_improve.rds",
                                 "delta",     "se_delta",     mu0 = 0)  # change, target 0

# -----------------------------------------------------------------------------
hosp_name <- hospital_names()

# top-20 (by the main model) rank-movement block for one estimand. The main-model
# rank is plain; each sensitivity column shows its rank and the change against the
# main model, reusing the shared movement helpers (1 = fastest / most improved).
format_block <- function(est, top_n = 20) {
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(est))
  else unname(hosp_name[as.character(est$diag_hosp)])
  est$hosp_name <- ifelse(is.na(nm) | nm == "", as.character(est$diag_hosp), nm)
  
  rk <- data.frame(main = comp_rank(est[[main_title]]))
  for (col in sens_cols) rk[[col]] <- comp_rank(est[[col]])
  base_rk <- rk$main
  keep    <- order(base_rk, est$hosp)[seq_len(min(top_n, nrow(est)))]
  
  text <- data.frame(Hospital = est$hosp_name[keep],
                     `Hospital site code` = as.character(est$diag_hosp[keep]),
                     check.names = FALSE, stringsAsFactors = FALSE)
  text[[main_title]] <- as.character(base_rk[keep])
  move <- matrix(0, nrow = length(keep), ncol = length(sens_cols),
                 dimnames = list(NULL, sens_cols))
  for (col in sens_cols) {
    mr <- rk[[col]][keep]
    text[[col]] <- move_cell(mr, base_rk[keep])
    move[, col]  <- base_rk[keep] - mr
  }
  list(text = text, move = move)
}

sus_block <- format_block(sus_est)
imp_block <- format_block(imp_est)
cat("\nTop 20 sustained (main-model rank):\n");   print(sus_block$text, row.names = FALSE)
cat("\nTop 20 improvement (main-model rank):\n"); print(imp_block$text, row.names = FALSE)
write.csv(sus_block$text, file.path(out_dir, "prior_sensitivity_sustained.csv"), row.names = FALSE)
write.csv(imp_block$text, file.path(out_dir, "prior_sensitivity_improve.csv"),   row.names = FALSE)

# assemble one Word table: a sustained banner, column titles, the sustained rows,
# then an improvement banner, titles again, and the improvement rows ----------
TABLE_FONT_SIZE <- 8
col_names <- names(sus_block$text)     # Hospital, Hospital site code, Main model, the 4 priors

blank_row <- function(first_cell = "") {
  r <- as.data.frame(as.list(rep("", length(col_names))),
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(r) <- col_names; r[[1]] <- first_cell; r
}
titles <- as.data.frame(as.list(col_names), stringsAsFactors = FALSE, check.names = FALSE)
names(titles) <- col_names

sus_div  <- blank_row("Sustained performance")
imp_div  <- blank_row("Improvement over the period")
combined <- rbind(sus_div, titles, sus_block$text, imp_div, titles, imp_block$text)

n_sus <- nrow(sus_block$text); n_imp <- nrow(imp_block$text)
row_sus_div <- 1
row_titles1 <- 2
sus_rows    <- row_titles1 + seq_len(n_sus)
row_imp_div <- max(sus_rows) + 1
row_titles2 <- row_imp_div + 1
imp_rows    <- row_titles2 + seq_len(n_imp)

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")

# colour only the "(arrow N)" suffix of each sensitivity cell; the rank stays black
colour_moves <- function(ft, block_text, block_move, rows) {
  for (col in sens_cols) {
    for (r in seq_along(rows)) {
      cell <- block_text[[col]][r]; mv <- block_move[r, col]
      if (is.na(mv) || mv == 0 || !grepl("\\(", cell)) next
      ft <- compose(ft, i = rows[r], j = col, part = "body",
                    value = as_paragraph(
                      as_chunk(cell_prefix(cell),
                               props = fp_text(color = "black", font.size = TABLE_FONT_SIZE)),
                      as_chunk(cell_suffix(cell),
                               props = fp_text(color = move_colour(mv), font.size = TABLE_FONT_SIZE))))
    }
  }
  ft
}
ft <- colour_moves(ft, sus_block$text, sus_block$move, sus_rows)
ft <- colour_moves(ft, imp_block$text, imp_block$move, imp_rows)

ft <- merge_at(ft, i = row_sus_div, j = seq_along(col_names), part = "body")
ft <- merge_at(ft, i = row_imp_div, j = seq_along(col_names), part = "body")
ft <- bold(ft, i = c(row_sus_div, row_titles1, row_imp_div, row_titles2), part = "body")
ft <- align(ft, part = "body", align = "center")
ft <- align(ft, j = 1:2, align = "left", part = "body")
ft <- align(ft, i = c(row_sus_div, row_imp_div), align = "left", part = "body")

rule <- fp_border(color = "black", width = 1)
ft <- border_outer(ft, border = rule, part = "body")
ft <- hline(ft, i = row_sus_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles1,   border = rule, part = "body")
ft <- hline(ft, i = max(sus_rows), border = rule, part = "body")
ft <- hline(ft, i = row_imp_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles2,   border = rule, part = "body")

faint <- fp_border(color = "grey85", width = 0.5)
ft <- hline(ft, i = sus_rows[-length(sus_rows)], border = faint, part = "body")
ft <- hline(ft, i = imp_rows[-length(imp_rows)], border = faint, part = "body")

# dotted verticals separating the id columns and the main-model rank from the priors
vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
vcols <- match(c("Hospital site code", main_title), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

ft <- autofit(ft)
ft <- width(ft, j = "Hospital",           width = 1.95)
ft <- width(ft, j = "Hospital site code", width = 0.7)
ft <- width(ft, j = main_title,           width = 1.3)   # titles wrap to a few lines
ft <- width(ft, j = sens_cols,            width = 1.45)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Prior-sensitivity of the shrinkage ranking: sustained (upper block) and",
  "improvement (lower block). The main-model column gives each hospital's rank",
  "under the headline shrinkage prior; each sensitivity column gives its rank and",
  "the change against the main model, coloured green (up) to red (down). All",
  "columns are shrinkage estimates; only the prior on the between-hospital sd (tau)",
  "differs, and is written out in each column title.")
footnote_text <- paste0(
  "All columns use the same normal-normal shrinkage model and the same prior on ",
  "the grand mean, which is unchanged throughout: mu ~ Normal(prior mean, ",
  prior_mu_sd, ") in days, where the prior mean is the data mean for the sustained ",
  "estimand and 0 for the change. Only the prior on the between-hospital sd (tau) ",
  "differs, as given in each column title. DerSimonian-Laird is not Bayesian: it ",
  "plugs in the method-of-moments estimate of tau and shrinks analytically, with ",
  "its expected rank taken from simulated draws of that Gaussian posterior. Ranks ",
  "are by expected posterior rank (1 = fastest).")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.5, right = 0.5)))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "prior_sensitivity_methods.docx"))
cat("prior-sensitivity table written to prior_sensitivity_methods.docx\n")