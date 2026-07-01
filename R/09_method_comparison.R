# 09  ranking-method comparison and funnel plots
# -----------------------------------------------------------------------------
# Compares hospital rankings on the sustained and improvement estimands (the
# basis for Table 3), draws funnel plots, and builds rank-movement tables.
#
# The comparison is deliberately built from two clean steps, so every method is
# on the same footing:
#   step 1 - generate a per-hospital mean waiting time, either case-mix adjusted
#            (balancing-weights direct standardisation, the main model) or not
#            adjusted at all (raw mean);
#   step 2 - optionally shrink those means with the one Bayesian normal-normal
#            routine (fit_shrink, in 01_config.R) - the same routine that feeds
#            the headline figures.
# So each generator appears twice, unshrunk and shrunk, and the shrinkage is
# identical across them. Indirect standardisation (observed minus expected) is
# also shown, as the basis for the funnel plots.
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

cv   <- code_covariates(df)          # the main model case-mix (age + cci)
df   <- cv$data
cont <- cv$cont; bin <- cv$bin

# model-based direct standardisation (g-computation). Fit a fixed-effect model
# with hospital coded outright (no intercept, one term per hospital) plus the
# same age + cci covariates the weighting uses. The direct-standardised mean for
# a hospital is the average predicted wait if the whole sample were treated at
# that hospital, which for this linear, additive model is that hospital's own
# effect plus the population-mean covariate vector times the covariate
# coefficients - a linear contrast of the fitted coefficients. Its standard error
# is that contrast's se from the model covariance (L V L'), so the estimate and se
# feed the same Bayesian shrinkage as the weighted estimate.
regstd_gcomp <- function(d, cont, bin) {
  d <- d %>% mutate(hospf = factor(hosp))
  covs <- c(cont, bin)
  form <- as.formula(paste("wait ~ 0 + hospf +", paste(covs, collapse = " + ")))
  m <- lm(form, data = d)
  b <- coef(m); V <- vcov(m)
  hs   <- levels(d$hospf)
  xbar <- colMeans(d[covs])                       # population mean covariate vector
  L <- matrix(0, nrow = length(hs), ncol = length(b), dimnames = list(hs, names(b)))
  for (i in seq_along(hs)) {
    L[i, paste0("hospf", hs[i])] <- 1             # this hospital's own effect
    L[i, names(xbar)]            <- xbar          # averaged over the whole case-mix
  }
  data.frame(hosp = as.integer(hs),
             reg_mean = as.numeric(L %*% b),
             reg_se   = sqrt(diag(L %*% V %*% t(L))))
}

# the improvement analogue: the change in the direct-standardised mean between
# the two periods. With additive covariates the covariate part cancels in the
# difference, so the change is the contrast (second - first) of each hospital's
# period-specific effect, again with its se from the model covariance.
regstd_change <- function(d, cont, bin) {
  d <- d %>% mutate(hp = factor(paste(hosp, period, sep = "_")))
  covs <- c(cont, bin)
  form <- as.formula(paste("wait ~ 0 + hp +", paste(covs, collapse = " + ")))
  m <- lm(form, data = d)
  b <- coef(m); V <- vcov(m)
  hosps <- sort(unique(d$hosp))
  L <- matrix(0, nrow = length(hosps), ncol = length(b), dimnames = list(hosps, names(b)))
  for (i in seq_along(hosps)) {
    L[i, paste0("hp", hosps[i], "_second")] <-  1
    L[i, paste0("hp", hosps[i], "_first")]  <- -1
  }
  data.frame(hosp = hosps,
             reg_mean = as.numeric(L %*% b),
             reg_se   = sqrt(diag(L %*% V %*% t(L))))
}

# indirect standardisation: expected wait from a pooled case-mix model, then
# observed - expected + grand mean. Also the basis for the sustained funnel plot.
pm <- lm(as.formula(paste("wait ~", paste(c(cont, bin), collapse = " + "))),
         data = df)
df$pred <- predict(pm)
grand   <- mean(df$wait)
ind <- df %>% group_by(hosp) %>%
  summarise(n = n(), obs = mean(wait), exp = mean(pred), .groups = "drop") %>%
  mutate(indirect = obs - exp + grand)

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
# headline. Returns the shrunk posterior mean rank (1 = fastest).
raw_shr <- stan_shrink_rank(site$raw_mean, site$se_raw)
site$raw_shrunk_rank <- raw_shr$exp_rank

# model-based direct standardisation and its shrunk rank (same routine again)
reg_sus <- regstd_gcomp(df, cont, bin)
reg_sus$reg_shrunk_rank <- stan_shrink_rank(reg_sus$reg_mean, reg_sus$reg_se)$exp_rank

comp <- site %>%
  transmute(hosp, diag_hosp,
            wt_mean = stand_adj,       # weighted direct-standardised mean (unshrunk)
            raw_mean,                  # raw mean (unshrunk)
            raw_shrunk_rank) %>%       # raw mean, shrunk
  left_join(reg_sus %>% select(hosp, reg = reg_mean, reg_shrunk_rank), by = "hosp") %>%
  left_join(ind %>% select(hosp, indirect), by = "hosp") %>%
  left_join(sus %>% select(hosp, shrunk_rank = exp_rank), by = "hosp")   # weighted, shrunk (headline)

# rank agreement between methods (1 = fastest); the shrunk columns are already
# posterior ranks, so rank() just re-expresses them on the same scale.
methods <- c("raw_mean", "wt_mean", "reg", "indirect")
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
# change against the primary. The primary is the balancing-weights (age + cci)
# model with Bayesian shrinkage - the ranking that feeds the main figures.
# Sustained and improvement are shown as one continuous Word table: the broad
# method groupings appear once at the top, the individual column titles repeat
# before the improvement block, and only the "(arrow N)" part of each cell is
# coloured (green up, red down); the rank number itself stays plain black.
TABLE_FONT_SIZE <- 7               # one place, so plain and coloured text match

up_arrow   <- intToUtf8(8593)      # built at run time so the script stays ASCII
down_arrow <- intToUtf8(8595)

comp_rank <- function(x) rank(x, ties.method = "min")   # 1 = best (lowest value)

# format one method cell: the method's rank, and its move relative to the base.
# move > 0 means the method ranks the hospital nearer the top than the primary.
move_cell <- function(method_rank, base_rank) {
  move <- base_rank - method_rank
  out  <- sprintf("%d (=)", method_rank)
  up <- move > 0; down <- move < 0
  out[up]   <- sprintf("%d (%s%d)", method_rank[up],   up_arrow,    move[up])
  out[down] <- sprintf("%d (%s%d)", method_rank[down], down_arrow, -move[down])
  out
}

# green (moved up) to red (moved down) colour for a signed rank move, used as a
# font colour on just the "(arrow N)" part of a cell. Six colours plus black for
# no change, so the size of the move is visible at a glance.
move_colour <- function(m) {
  pick <- function(x) {
    if (is.na(x))  return("black")
    if (x >=  6)   return("#1A9850")   # up 6+   strong green
    if (x >=  3)   return("#66BD63")   # up 3-5  green
    if (x >=  1)   return("#4D9221")   # up 1-2  darker light-green (legible as text)
    if (x ==  0)   return("black")     # no change
    if (x >= -2)   return("#B35806")   # down 1-2 darker light-orange (legible as text)
    if (x >= -5)   return("#F46D43")   # down 3-5 orange-red
    "#D73027"                          # down 6+  red
  }
  vapply(m, pick, character(1))
}

# split a formatted cell like "2 (up_arrow1)" into the plain rank prefix "2 " and
# the parenthetical suffix "(up_arrow1)"; used to colour only the suffix.
cell_prefix <- function(x) sub("\\(.*$", "", x)
cell_suffix <- function(x) sub("^[^(]*", "", x)

# hospital display names (Title Case, corrected, code suffix stripped) come from
# the shared helper in 01_config.R, keyed by canonical site code.
hosp_name <- hospital_names()

# shared column set: display title -> estimate column, in table order. Three
# direct-standardisation generators (weighted, model-based, and none), each shown
# unshrunk and shrunk, then the indirect (funnel) estimate.
method_cols <- c(
  "weighted mean + shrinkage" = "shrunk_rank",
  "weighted mean"             = "wt_mean",
  "model mean + shrinkage"    = "reg_shrunk_rank",
  "model mean"                = "reg",
  "raw mean + shrinkage"      = "raw_shrunk_rank",
  "raw mean"                  = "raw_mean",
  "indirect"                  = "indirect"
)
base_col <- "shrunk_rank"

# format the top-n hospitals for one estimand. Returns the display text (Hospital
# name, site code, one column per method) and the matching matrix of signed moves
# used to colour the cells.
format_block <- function(est, top_n = 20) {
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(est))
  else unname(hosp_name[as.character(est$diag_hosp)])
  est$hosp_name <- ifelse(is.na(nm) | nm == "", as.character(est$diag_hosp), nm)
  est$hosp_code <- as.character(est$diag_hosp)
  
  present <- method_cols[method_cols %in% names(est)]
  rk <- data.frame(hosp = est$hosp)
  for (col in present) rk[[col]] <- comp_rank(est[[col]])
  base_rk <- rk[[base_col]]
  keep    <- order(base_rk, est$hosp)[seq_len(min(top_n, nrow(est)))]
  
  text <- data.frame(Hospital = est$hosp_name[keep], `Site code` = est$hosp_code[keep],
                     stringsAsFactors = FALSE, check.names = FALSE)
  move <- matrix(NA_real_, nrow = length(keep), ncol = length(method_cols),
                 dimnames = list(NULL, names(method_cols)))
  for (lab in names(method_cols)) {
    col <- method_cols[[lab]]
    if (!(col %in% names(est))) { text[[lab]] <- "-"; next }
    mr <- rk[[col]][keep]
    text[[lab]] <- if (col == base_col) as.character(mr) else move_cell(mr, base_rk[keep])
    move[, lab] <- if (col == base_col) 0 else base_rk[keep] - mr
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

# model-based change and its shrunk rank (same routine, centred on no change)
reg_imp <- regstd_change(di, cont, bin)
reg_imp$reg_shrunk_rank <- stan_shrink_rank(reg_imp$reg_mean, reg_imp$reg_se, mu_mean = 0)$exp_rank

# indirect change (best guess): an indirect standardisation is not uniquely
# defined for a change score, so take the difference-in-differences of the
# observed-minus-expected gap - each period's observed mean minus what the
# pooled case-mix model expects, then later minus first.
di$pred <- predict(pm, newdata = di)
ind_chg <- di %>%
  group_by(hosp, period) %>%
  summarise(gap = mean(wait) - mean(pred), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = gap) %>%
  transmute(hosp, indirect = second - first)

comp_imp <- imp_wt %>%
  transmute(hosp, diag_hosp, wt_mean = delta) %>%
  left_join(raw_change %>% select(hosp, raw_mean, raw_shrunk_rank), by = "hosp") %>%
  left_join(reg_imp %>% select(hosp, reg = reg_mean, reg_shrunk_rank), by = "hosp") %>%
  left_join(ind_chg, by = "hosp") %>%
  left_join(impr_shrunk %>% select(hosp, shrunk_rank = exp_rank), by = "hosp")

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
# again (with * on the improvement approximation), then the improvement rows.
# The caption and the footnote sit outside the table as ordinary text.
col_names <- names(sus_block$text)                  # "Hospital","Site code", the 5 methods

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
group_row[["indirect"]]                  <- "Funnel plot (indirect standardisation)"

titles_top <- as.data.frame(as.list(col_names), stringsAsFactors = FALSE, check.names = FALSE)
names(titles_top) <- col_names
titles_imp <- titles_top
titles_imp[["indirect"]] <- "indirect*"

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
# number itself stays plain black and the same size as the coloured text.
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
# fill colour anywhere; sections are shown with rule lines instead. "Funnel plot"
# is a single column, so no merge is needed for it.
ft <- merge_at(ft, i = row_group, j = match(c("Hospital", "Site code"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("weighted mean + shrinkage", "weighted mean"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("model mean + shrinkage", "model mean"), col_names), part = "body")
ft <- merge_at(ft, i = row_group,
               j = match(c("raw mean + shrinkage", "raw mean"), col_names), part = "body")
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

# faint dotted vertical rules splitting the broad sections (and the id columns
# from the methods); only on the structured rows, so they do not cut through the
# merged banner rows.
vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
vcols <- match(c("Site code", "weighted mean", "model mean", "raw mean"), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

# column widths sized to fit landscape: narrow the ranking columns and let their
# titles wrap; keep Hospital as wide as fits so most names stay on one line.
# autofit() first, then override, then fix the layout so Word keeps these widths.
ft <- autofit(ft)
ft <- width(ft, j = "Hospital", width = 1.8)
ft <- width(ft, j = "Site code", width = 0.6)
ft <- width(ft, j = names(method_cols), width = 0.62)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

# caption above the table (plain document text, not the table's small font) and
# a footnote below explaining the one improvement approximation. Built with
# read_docx() rather than save_as_docx() so the caption and footnote can sit
# outside the flextable as ordinary paragraphs.
caption_text <- paste(
  "Hospital rankings by method: sustained (upper block) and improvement (lower",
  "block). Each cell gives the method's rank and its change against the primary",
  "model; the change is coloured green (up) to red (down).")
footnote_text <- paste(
  "* For the improvement estimand the indirect column is an approximation: an",
  "indirect standardisation of a change score is not uniquely defined. Shown as",
  "the change in the observed-minus-expected gap between the two periods.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "landscape")))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "movement_methods.docx"))

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