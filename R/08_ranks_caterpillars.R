# Posterior processing: shrunk hospital means with credible intervals,
# expected ranks, probability of being in the top X%, caterpillar plots, and
# selection of sustained / improved positive-deviant candidates. The rank logic
# follows the template (rank within each posterior draw), with shorter waits as
# better so rank 1 is the best performer.

library(rstan)
library(dplyr)
library(ggplot2)
library(flextable)
library(officer)

source("R/01_config.R")

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

# top-20 Word table: sustained and improvement in one document ----------------
# One row per hospital for the fastest / most-improved 20, with the shrinkage
# model estimate (standardised waiting days) and the posterior probability of
# sitting in the top 10 / 20 / 25 percent. Built like the comparison table in 09:
# the probability heading appears once, the column titles repeat before each
# block, and the caption and footnote sit outside the table as ordinary text.
hosp_name <- hospital_names()

top_table <- function(d, top_n = 20) {
  d <- d %>% arrange(post_mean) %>% mutate(rank = row_number())
  d <- d[seq_len(min(top_n, nrow(d))), ]
  nm <- if (is.null(hosp_name)) rep(NA_character_, nrow(d))
  else unname(hosp_name[as.character(d$diag_hosp)])
  data.frame(
    Rank     = as.character(d$rank),
    Hospital = ifelse(is.na(nm) | nm == "", as.character(d$diag_hosp), nm),
    SiteCode = as.character(d$diag_hosp),
    Patients = as.character(d$patients),
    Estimate = sprintf("%.1f", d$post_mean),
    p10      = sprintf("%.2f", d$p_top10),
    p20      = sprintf("%.2f", d$p_top20),
    p25      = sprintf("%.2f", d$p_top25),
    stringsAsFactors = FALSE, check.names = FALSE)
}

sus$patients <- sus$n            # patients over the window
imp$patients <- imp$n1 + imp$n2  # patients over both halves
sus_tab <- top_table(sus)
imp_tab <- top_table(imp)

col_names <- c("Rank", "Hospital", "SiteCode", "Patients", "Estimate", "p10", "p20", "p25")
blank_row <- function(first_cell = "") {
  r <- as.data.frame(as.list(rep("", length(col_names))), stringsAsFactors = FALSE)
  names(r) <- col_names
  r[[1]] <- first_cell
  r
}

# broad heading over the three probability columns; the rest of the row is blank
group_row <- blank_row("")
group_row[["p10"]] <- "Probability in the top X% of performers"

# the column titles; the estimate label differs between the two blocks
titles_row <- function(est_label) {
  r <- as.data.frame(as.list(c("Rank", "Hospital", "Site code", "Patients (n)",
                               est_label, "10%", "20%", "25%")), stringsAsFactors = FALSE)
  names(r) <- col_names
  r
}
titles_sus <- titles_row("Standardised wait (days)")
titles_imp <- titles_row("Change in wait (days)")

sus_divider <- blank_row("Sustained performance (average waiting time, 2020-2021)")
imp_divider <- blank_row("Improvement over the period (change, 2021 vs 2020)")

combined <- rbind(group_row, sus_divider, titles_sus, sus_tab,
                  imp_divider, titles_imp, imp_tab)

n_sus <- nrow(sus_tab); n_imp <- nrow(imp_tab)
row_group   <- 1
row_sus_div <- 2
row_titles1 <- 3
sus_rows    <- row_titles1 + seq_len(n_sus)
row_imp_div <- max(sus_rows) + 1
row_titles2 <- row_imp_div + 1
imp_rows    <- row_titles2 + seq_len(n_imp)

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")   # every row above is already in the body

# merge the probability heading over its three columns, the blank left part of
# the heading row, and each banner across the full width
ft <- merge_at(ft, i = row_group,
               j = match(c("Rank", "Hospital", "SiteCode", "Patients", "Estimate"), col_names),
               part = "body")
ft <- merge_at(ft, i = row_group, j = match(c("p10", "p20", "p25"), col_names), part = "body")
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

# faint dotted rules separating identifiers, estimate and probabilities; only on
# the title and data rows, so they do not cut through the merged heading / banners
vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
vcols <- match(c("Patients", "Estimate"), col_names)
structured <- setdiff(seq_len(nrow(combined)), c(row_group, row_sus_div, row_imp_div))
ft <- vline(ft, i = structured, j = vcols, border = vrule, part = "body")

ft <- fontsize(ft, size = 9, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

ft <- autofit(ft)
ft <- width(ft, j = "Rank", width = 0.45)
ft <- width(ft, j = "Hospital", width = 1.75)
ft <- width(ft, j = "SiteCode", width = 0.55)
ft <- width(ft, j = "Patients", width = 0.6)
ft <- width(ft, j = "Estimate", width = 1.0)
ft <- width(ft, j = c("p10", "p20", "p25"), width = 0.5)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Top 20 hospitals by the shrinkage model: sustained (upper block) and",
  "improvement (lower block). The estimate is the age- and comorbidity-",
  "standardised waiting time from the Bayesian shrinkage model; the probabilities",
  "are the posterior probability the hospital sits in the top X% of performers.")
footnote_text <- paste(
  "For the improvement block the estimate is the change in standardised waiting",
  "time (a negative value = faster over the period) and 'top X%' means the most",
  "improved X% of hospitals.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "portrait")))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "top20_ranking.docx"))
cat("top-20 ranking table saved.\n")

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