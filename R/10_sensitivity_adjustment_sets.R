# 10  sensitivity to the standardisation adjustment set
# -----------------------------------------------------------------------------
# The main model standardises for age and comorbidity only. This script asks
# whether adding season (quarter of diagnosis) and calendar year to the direct
# standardisation changes the picture, in three ways:
#   Part A - does adding season and calendar year to the outcome model improve
#            fit, judged by AIC and a likelihood-ratio test against the base;
#   Part B - a Word table (balancing-weights direct standardisation, shrunk) for
#            the sustained and improvement estimands, giving for the base case
#            (age + cci) and for the base plus calendar year and season each
#            hospital's shrunk rank and its shrunk mean waiting time (s.e.), with
#            the year+season ranks coloured by their change against the base case;
#   Part C - the rank correlation across the two sets.
# Season and calendar year are held OUT of the main model; this is where their
# effect on the ranking is examined.

library(balancer)
library(dplyr)
library(rstan)
library(flextable)
library(officer)

source("R/01_config.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

df <- readRDS(file.path(out_dir, "analysis_data.rds"))
cv <- code_covariates(df); df <- cv$data
cont <- cv$cont; bin <- cv$bin
season <- c("q2", "q3", "q4")     # quarter-of-diagnosis dummies
year   <- c("yr_late")            # later calendar year of the window

# Part A  outcome-model fit: does season / year earn its degrees of freedom? ---
# Base is the main model (age + cci). Each candidate adds terms that CONTAIN the
# base as a special case, so a likelihood-ratio test against the base is valid; a
# negative dAIC with a small p-value means the extra terms improve the fit of
# wait on the adjustment set.
base  <- lm(reformulate(c(cont, bin),               "wait"), df)
m_sea <- lm(reformulate(c(cont, bin, season),       "wait"), df)
m_yr  <- lm(reformulate(c(cont, bin, year),         "wait"), df)
m_all <- lm(reformulate(c(cont, bin, season, year), "wait"), df)

models <- list(base, m_sea, m_yr, m_all)
labels <- c("base: age + cci", "+ season", "+ calendar year", "+ season + year")

base_aic   <- AIC(base)
ll_base    <- as.numeric(logLik(base))
dfres_base <- df.residual(base)

aic_tab <- data.frame(model = character(), params = integer(), AIC = numeric(),
                      dAIC = numeric(), LR_chisq = numeric(), LR_df = integer(),
                      p_value = numeric(), stringsAsFactors = FALSE)
for (i in seq_along(models)) {
  m <- models[[i]]; aic <- AIC(m)
  if (labels[i] == "base: age + cci") {
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
  aic_tab <- rbind(aic_tab, row)
}
show <- aic_tab
show$AIC <- round(show$AIC, 1); show$dAIC <- round(show$dAIC, 1)
show$LR_chisq <- round(show$LR_chisq, 2); show$p_value <- signif(show$p_value, 3)
cat("Adjustment-set fit (likelihood-ratio test vs age + cci base):\n")
print(show, row.names = FALSE)
write.csv(aic_tab, file.path(out_dir, "adjustment_set_aic_lrt.csv"), row.names = FALSE)

# Part B  shrunk ranks and means under two adjustment sets ---------------------
# base case: reuse the headline shrinkage outputs (age + cci). year + season:
# standardise afresh with the extra covariates, then shrink with the same routine.
hosp_name <- hospital_names()
base_sus <- read.csv(file.path(out_dir, "ranks_sustained.csv"))  # post_mean, post_sd, exp_rank
base_imp <- read.csv(file.path(out_dir, "ranks_improve.csv"))

# sustained, adding season and calendar year to the age + comorbidity covariates
ys_sus <- run_standardise(patient_data          = df,
                          continuous_covariates = cont,
                          binary_covariates     = c(bin, season, year))$site
ys_sus <- bind_cols(ys_sus, stan_shrink_rank(ys_sus$stand_adj, ys_sus$se_adj_pool))

# improvement, adding season and calendar year (the later-year indicator is
# constant within a half and is dropped, so this is effectively the + season change)
ys_imp <- standardise_change(patient_data          = df,
                             continuous_covariates = cont,
                             binary_covariates     = c(bin, season, year))
ys_imp <- bind_cols(ys_imp, stan_shrink_rank(ys_imp$delta, ys_imp$se_delta, mu_mean = 0))

# join the base and year+season sets for one estimand, and add competition ranks
join_sets <- function(bd, yd) {
  d <- bd %>%
    transmute(hosp, diag_hosp, base_mean = post_mean, base_sd = post_sd, base_er = exp_rank) %>%
    inner_join(yd %>% select(hosp, ys_mean = post_mean, ys_sd = post_sd, ys_er = exp_rank),
               by = "hosp")
  d$r_base <- comp_rank(d$base_er)
  d$r_ys   <- comp_rank(d$ys_er)
  d
}
sus_full <- join_sets(base_sus, ys_sus)
imp_full <- join_sets(base_imp, ys_imp)

# the top-20 display rows for one estimand: base rank (reference, plain) and its
# shrunk mean (s.e.); year+season rank (move vs base, coloured) and its mean (s.e.)
block_disp <- function(d, top_n = 20) {
  d <- d[order(d$r_base), ]
  d <- d[seq_len(min(top_n, nrow(d))), ]
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(d))
  else unname(hosp_name[as.character(d$diag_hosp)])
  disp <- data.frame(
    Hospital  = ifelse(is.na(nm) | nm == "", as.character(d$diag_hosp), nm),
    SiteCode  = as.character(d$diag_hosp),
    base_rank = as.character(d$r_base),
    base_ms   = sprintf("%.1f (%.1f)", d$base_mean, d$base_sd),
    ys_rank   = move_cell(d$r_ys, d$r_base),
    ys_ms     = sprintf("%.1f (%.1f)", d$ys_mean, d$ys_sd),
    stringsAsFactors = FALSE, check.names = FALSE)
  list(disp = disp, ys_move = d$r_base - d$r_ys)
}
sus_block <- block_disp(sus_full)
imp_block <- block_disp(imp_full)
write.csv(sus_block$disp, file.path(out_dir, "adjustment_set_sustained.csv"), row.names = FALSE)
write.csv(imp_block$disp, file.path(out_dir, "adjustment_set_improve.csv"),   row.names = FALSE)

# assemble one table: broad headings once, column titles repeated before each
# block, sustained then improvement, with sustained/improvement banner rows.
TABLE_FONT_SIZE <- 8
col_names <- names(sus_block$disp)   # 6 columns

blank_row <- function(first_cell = "") {
  r <- as.data.frame(as.list(rep("", length(col_names))), stringsAsFactors = FALSE)
  names(r) <- col_names
  r[[1]] <- first_cell
  r
}
group_row <- blank_row()
group_row[["base_rank"]] <- "Base: age + comorbidity"
group_row[["ys_rank"]]   <- "+ calendar year and season"

sub <- c("rank", "mean (s.e.), days")
titles_row <- as.data.frame(as.list(c("Hospital", "Site code", rep(sub, 2))),
                            stringsAsFactors = FALSE)
names(titles_row) <- col_names

sus_divider <- blank_row("Sustained performance (average waiting time, 2020-2021)")
imp_divider <- blank_row("Improvement over the period (change, 2021 vs 2020)")

combined <- rbind(group_row, sus_divider, titles_row, sus_block$disp,
                  imp_divider, titles_row, imp_block$disp)

n_sus <- nrow(sus_block$disp); n_imp <- nrow(imp_block$disp)
row_group   <- 1
row_sus_div <- 2
row_titles1 <- 3
sus_rows    <- row_titles1 + seq_len(n_sus)
row_imp_div <- max(sus_rows) + 1
row_titles2 <- row_imp_div + 1
imp_rows    <- row_titles2 + seq_len(n_imp)

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")

# colour only the parenthetical move on the year+season rank column
colour_ys <- function(ft, disp, move, rows) {
  for (r in seq_along(rows)) {
    cell <- disp$ys_rank[r]; mv <- move[r]
    if (is.na(mv) || mv == 0 || !grepl("\\(", cell)) next
    ft <- compose(ft, i = rows[r], j = "ys_rank", part = "body",
                  value = as_paragraph(
                    as_chunk(cell_prefix(cell),
                             props = fp_text(color = "black", font.size = TABLE_FONT_SIZE)),
                    as_chunk(cell_suffix(cell),
                             props = fp_text(color = move_colour(mv), font.size = TABLE_FONT_SIZE))))
  }
  ft
}
ft <- colour_ys(ft, sus_block$disp, sus_block$ys_move, sus_rows)
ft <- colour_ys(ft, imp_block$disp, imp_block$ys_move, imp_rows)

# merge the two broad headings and the two banners
ft <- merge_at(ft, i = row_group, j = match(c("base_rank", "base_ms"), col_names), part = "body")
ft <- merge_at(ft, i = row_group, j = match(c("ys_rank", "ys_ms"), col_names), part = "body")
ft <- merge_at(ft, i = row_sus_div, j = seq_along(col_names), part = "body")
ft <- merge_at(ft, i = row_imp_div, j = seq_along(col_names), part = "body")

ft <- bold(ft, i = c(row_group, row_sus_div, row_titles1, row_imp_div, row_titles2), part = "body")
ft <- align(ft, part = "body", align = "center")
ft <- align(ft, j = match(c("Hospital", "SiteCode"), col_names), align = "left", part = "body")
ft <- align(ft, i = c(row_sus_div, row_imp_div), align = "left", part = "body")

rule <- fp_border(color = "black", width = 1)
ft <- border_outer(ft, border = rule, part = "body")
ft <- hline(ft, i = row_group,     border = rule, part = "body")
ft <- hline(ft, i = row_sus_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles1,   border = rule, part = "body")
ft <- hline(ft, i = max(sus_rows), border = rule, part = "body")
ft <- hline(ft, i = row_imp_div,   border = rule, part = "body")
ft <- hline(ft, i = row_titles2,   border = rule, part = "body")

# faint dotted rule between the id columns and each group
vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
vcols <- match(c("SiteCode", "base_ms"), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_group, row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

ft <- autofit(ft)
ft <- width(ft, j = "Hospital", width = 1.8)
ft <- width(ft, j = "SiteCode", width = 0.55)
ft <- width(ft, j = c("base_rank", "ys_rank"), width = 0.7)
ft <- width(ft, j = c("base_ms", "ys_ms"), width = 1.0)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Sensitivity to the standardisation adjustment set (balancing-weights direct",
  "standardisation, shrunk): sustained (upper block) and improvement (lower",
  "block), top 20 by the base-case rank. For the base case (age + comorbidity)",
  "and for the base plus calendar year and season, each hospital's shrunk rank",
  "and shrunk mean waiting time (s.e.) in days are shown; the year+season ranks",
  "are coloured green (up) to red (down) by their change against the base case.")
footnote_text <- paste(
  "Ranks are competition ranks (1 = fastest); the base-case column is the",
  "reference. For the improvement block the later-year indicator is constant",
  "within each half and is dropped, so that column is effectively the change",
  "with season added.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "portrait")))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "sensitivity_adjustment_sets.docx"))

# Part C  rank agreement between the two sets ---------------------------------
cat("\nSpearman rank correlation, base vs + calendar year + season (shrunk):\n")
cat(sprintf("  sustained:   %.3f\n", cor(sus_full$r_base, sus_full$r_ys, method = "spearman")))
cat(sprintf("  improvement: %.3f\n", cor(imp_full$r_base, imp_full$r_ys, method = "spearman")))

# fastest-quintile stability: is the fastest 20% under the base case still the
# fastest 20% once calendar year and season are added?
in_q <- function(r, n) r <= ceiling(0.20 * n)
retain <- function(d) {
  b <- in_q(d$r_base, nrow(d)); y <- in_q(d$r_ys, nrow(d))
  100 * sum(b & y) / sum(b)
}
cat(sprintf("\nFastest-quintile retained after adding year + season: sustained %.0f%%, improvement %.0f%%\n",
            retain(sus_full), retain(imp_full)))
cat("\nadjustment-set sensitivity outputs written.\n")