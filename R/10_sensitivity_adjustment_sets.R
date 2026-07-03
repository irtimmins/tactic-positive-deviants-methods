# 10  sensitivity of the ranking to the calendar adjustment
# -----------------------------------------------------------------------------
# The primary model standardises for age, comorbidity, season and calendar year.
# This script asks how much the hospital ranking depends on the calendar part of
# that set: it re-standardises with age and comorbidity only (season and calendar
# year removed) and compares. The case for including season and calendar year, on
# outcome fit and effective sample size, is made separately in script 15; here
# the question is purely what dropping them does to the ranks.
#   Part A - a Word table (balancing-weights direct standardisation, shrunk) for
#            the sustained and improvement estimands: for the primary set and for
#            age + comorbidity only, each hospital's shrunk rank and shrunk mean
#            waiting time (s.e.), with the age + comorbidity ranks coloured by
#            their change against the primary;
#   Part B - the rank correlation between the two sets.

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

# Part A  shrunk ranks and means under the primary and reduced sets -----------
# primary: reuse the headline shrinkage outputs (age + cci + season + calendar
# year). reduced: re-standardise with age and comorbidity only, then shrink with
# the same routine.
hosp_name <- hospital_names()
prim_sus <- read.csv(file.path(out_dir, "ranks_sustained.csv"))  # post_mean, post_sd, exp_rank
prim_imp <- read.csv(file.path(out_dir, "ranks_improve.csv"))

# sustained, age and comorbidity only (season and calendar year removed)
red_sus <- run_standardise(patient_data          = df,
                           continuous_covariates = cont,
                           binary_covariates     = bin)$site
red_sus <- bind_cols(red_sus, stan_shrink_rank(red_sus$stand_adj, red_sus$se_adj_pool))

# improvement, age and comorbidity only
red_imp <- standardise_change(patient_data          = df,
                              continuous_covariates = cont,
                              binary_covariates     = bin)
red_imp <- bind_cols(red_imp, stan_shrink_rank(red_imp$delta, red_imp$se_delta, mu_mean = 0))

# join the primary (reference) and reduced sets for one estimand, and add
# competition ranks. r_prim is the reference rank; r_red is age + comorbidity only.
join_sets <- function(reference, reduced) {
  d <- reference %>%
    transmute(hosp, diag_hosp, prim_mean = post_mean, prim_sd = post_sd, prim_er = exp_rank) %>%
    inner_join(reduced %>% select(hosp, red_mean = post_mean, red_sd = post_sd, red_er = exp_rank),
               by = "hosp")
  d$r_prim <- comp_rank(d$prim_er)
  d$r_red  <- comp_rank(d$red_er)
  d
}
sus_full <- join_sets(prim_sus, red_sus)
imp_full <- join_sets(prim_imp, red_imp)

# the top-20 display rows for one estimand: primary rank (reference, plain) and
# its shrunk mean (s.e.); age + comorbidity rank (move vs primary, coloured) and
# its mean (s.e.)
block_disp <- function(d, top_n = 20) {
  d <- d[order(d$r_prim), ]
  d <- d[seq_len(min(top_n, nrow(d))), ]
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(d))
  else unname(hosp_name[as.character(d$diag_hosp)])
  disp <- data.frame(
    Hospital  = ifelse(is.na(nm) | nm == "", as.character(d$diag_hosp), nm),
    SiteCode  = as.character(d$diag_hosp),
    prim_rank = as.character(d$r_prim),
    prim_ms   = sprintf("%.1f (%.1f)", d$prim_mean, d$prim_sd),
    red_rank  = move_cell(d$r_red, d$r_prim),
    red_ms    = sprintf("%.1f (%.1f)", d$red_mean, d$red_sd),
    stringsAsFactors = FALSE, check.names = FALSE)
  list(disp = disp, red_move = d$r_prim - d$r_red)
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
group_row[["prim_rank"]] <- "Primary: age + comorbidity + season + calendar year"
group_row[["red_rank"]]  <- "Age + comorbidity only"

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

# colour only the parenthetical move on the age + comorbidity rank column
colour_red <- function(ft, disp, move, rows) {
  for (r in seq_along(rows)) {
    cell <- disp$red_rank[r]; mv <- move[r]
    if (is.na(mv) || mv == 0 || !grepl("\\(", cell)) next
    ft <- compose(ft, i = rows[r], j = "red_rank", part = "body",
                  value = as_paragraph(
                    as_chunk(cell_prefix(cell),
                             props = fp_text(color = "black", font.size = TABLE_FONT_SIZE)),
                    as_chunk(cell_suffix(cell),
                             props = fp_text(color = move_colour(mv), font.size = TABLE_FONT_SIZE))))
  }
  ft
}
ft <- colour_red(ft, sus_block$disp, sus_block$red_move, sus_rows)
ft <- colour_red(ft, imp_block$disp, imp_block$red_move, imp_rows)

# merge the two broad headings and the two banners
ft <- merge_at(ft, i = row_group, j = match(c("prim_rank", "prim_ms"), col_names), part = "body")
ft <- merge_at(ft, i = row_group, j = match(c("red_rank", "red_ms"), col_names), part = "body")
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
vcols <- match(c("SiteCode", "prim_ms"), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_group, row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

ft <- autofit(ft)
ft <- width(ft, j = "Hospital", width = 1.8)
ft <- width(ft, j = "SiteCode", width = 0.55)
ft <- width(ft, j = c("prim_rank", "red_rank"), width = 0.7)
ft <- width(ft, j = c("prim_ms", "red_ms"), width = 1.0)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Sensitivity of the ranking to the calendar adjustment (balancing-weights",
  "direct standardisation, shrunk): sustained (upper block) and improvement",
  "(lower block), top 20 by the primary rank. For the primary set (age +",
  "comorbidity + season + calendar year) and for age + comorbidity only, each",
  "hospital's shrunk rank and shrunk mean waiting time (s.e.) in days are shown;",
  "the age + comorbidity ranks are coloured green (up) to red (down) by their",
  "change against the primary.")
footnote_text <- paste(
  "Ranks are competition ranks (1 = fastest); the primary column is the",
  "reference. For the improvement block the later-year indicator is constant",
  "within each half and is dropped from the primary set, so the primary and the",
  "age + comorbidity columns there differ only by season.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "portrait")))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "sensitivity_adjustment_sets.docx"))

# Part B  rank agreement between the two sets ---------------------------------
cat("\nSpearman rank correlation, primary vs age + comorbidity only (shrunk):\n")
cat(sprintf("  sustained:   %.3f\n", cor(sus_full$r_prim, sus_full$r_red, method = "spearman")))
cat(sprintf("  improvement: %.3f\n", cor(imp_full$r_prim, imp_full$r_red, method = "spearman")))

# fastest-quintile stability: is the fastest 20% under the primary set still the
# fastest 20% once season and calendar year are removed?
in_q <- function(r, n) r <= ceiling(0.20 * n)
retain <- function(d) {
  p <- in_q(d$r_prim, nrow(d)); r <- in_q(d$r_red, nrow(d))
  100 * sum(p & r) / sum(p)
}
cat(sprintf("\nFastest-quintile retained after removing season + year: sustained %.0f%%, improvement %.0f%%\n",
            retain(sus_full), retain(imp_full)))
cat("\ncalendar-adjustment sensitivity outputs written.\n")