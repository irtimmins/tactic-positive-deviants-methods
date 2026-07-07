# 15  justifying the primary adjustment set
# -----------------------------------------------------------------------------
# The primary model adjusts for age (linear), comorbidity (0/1/2+), season and
# calendar year. This script justifies two choices behind that set and writes a
# supplementary Word table, then draws the effective-sample-size scatter for the
# primary model.
#
# Every candidate is judged two ways at once:
#   fit  - does adding the terms improve the outcome model, by AIC and by a
#          likelihood-ratio test against the age + comorbidity base;
#   cost - what does balancing on the larger set do to the effective sample size
#          (median and minimum across hospitals, and the mean percent reduction).
# Starting from age + comorbidity we test, separately, whether age needs a
# quadratic or cubic term (non-linearity), whether season and calendar year earn
# their place, and whether a finer comorbidity coding (0/1/2/3+) improves on the
# base 0/1/2+. The table lets fit and precision be read side by side, so a term
# is only kept if the fit gain is worth the sample-size it spends.

library(balancer)
library(dplyr)
library(ggplot2)
library(flextable)
library(officer)

source("R/01_config.R")

df <- readRDS(file.path(out_dir, "analysis_data.rds"))
cv <- code_covariates(df); df <- cv$data
cci <- cv$bin                                  # cci_1, cci_2p

# centred age powers: centring before squaring keeps the quadratic and cubic
# from lining up too closely with linear age. The column space is unchanged, so
# the likelihood-ratio test against linear age is exactly the same as with raw
# powers, but the balancing is better behaved.
age_c        <- df$agediag - mean(df$agediag)
df$age_sq    <- age_c^2
df$age_cube  <- age_c^3

# the candidate sets, in the order they appear in the table. cont / bin are the
# covariates handed to the weighting step; model_terms are the outcome-model
# predictors (the same covariates) used for AIC and the likelihood-ratio test.
base_cont <- cv$cont                           # agediag
base_bin  <- cci

# finer comorbidity coding (0/1/2/3+): split the base 2+ into separate 2 and 3+
# categories. The base 0/1/2+ is this model with the 2 and 3+ effects held equal,
# so the two are nested and the likelihood-ratio test against the base is valid
# (1 df); balancing on the extra category also lets its cost in ESS be seen.
df$cci_2  <- as.integer(df$cci_n_conditions == 2)
df$cci_3p <- as.integer(df$cci_n_conditions >= 3)
cci_fine  <- c(base_bin[1], "cci_2", "cci_3p")   # reuse the "1" dummy, split 2+

sets <- list(
  list(label = "Age + comorbidity (base)",
       cont = base_cont,                          bin = base_bin),
  list(label = "Base + age squared",
       cont = c(base_cont, "age_sq"),             bin = base_bin),
  list(label = "Base + age squared + age cubed",
       cont = c(base_cont, "age_sq", "age_cube"), bin = base_bin),
  list(label = "Base + season",
       cont = base_cont,                          bin = c(base_bin, season_terms)),
  list(label = "Base + calendar year",
       cont = base_cont,                          bin = c(base_bin, year_term)),
  list(label = "Base + season + calendar year",
       cont = base_cont,                          bin = c(base_bin, season_terms, year_term)),
  list(label = "Comorbidity 0/1/2/3+ (vs base 0/1/2+)",
       cont = base_cont,                          bin = cci_fine),
  list(label = "Comorbidity 0/1/2/3+ + season + calendar year",
       cont = base_cont,                          bin = c(cci_fine, season_terms, year_term)))

# fit the base outcome model once; every candidate is tested against it
base_model <- lm(reformulate(c(base_cont, base_bin), "wait"), df)
ll_base    <- as.numeric(logLik(base_model))
k_base     <- length(coef(base_model))
aic_base   <- AIC(base_model)

ess_summary <- function(ess_df) {
  data.frame(median_ess         = median(ess_df$ess),
             min_ess            = min(ess_df$ess),
             mean_pct_reduction = mean(ess_df$pct_reduction),
             n_below_threshold  = sum(ess_df$ess < ess_threshold))
}

# work through the candidates: outcome-model fit for the AIC / LR test, then a
# fresh standardisation for the effective-sample-size cost.
results <- data.frame()
for (s in sets) {
  model <- lm(reformulate(c(s$cont, s$bin), "wait"), df)
  k     <- length(coef(model))
  is_base <- identical(s$cont, base_cont) && identical(s$bin, base_bin)
  if (is_base) {
    lr_chisq <- NA_real_; lr_df <- NA_integer_; p_value <- NA_real_; d_aic <- 0
  } else {
    lr_chisq <- 2 * (as.numeric(logLik(model)) - ll_base)
    lr_df    <- k - k_base
    p_value  <- pchisq(lr_chisq, lr_df, lower.tail = FALSE)
    d_aic    <- AIC(model) - aic_base
  }
  
  ess <- tryCatch({
    fit <- run_standardise(patient_data          = df,
                           continuous_covariates = s$cont,
                           binary_covariates     = s$bin,
                           lambda                = lambda_main)
    ess_summary(hospital_ess(fit, s$label))
  }, error = function(e) {
    data.frame(median_ess = NA_real_, min_ess = NA_real_,
               mean_pct_reduction = NA_real_, n_below_threshold = NA_integer_)
  })
  
  results <- rbind(results, data.frame(
    label = s$label, params = k, AIC = AIC(model), dAIC = d_aic,
    LR_chisq = lr_chisq, LR_df = lr_df, p_value = p_value, ess,
    stringsAsFactors = FALSE))
}

write.csv(results, file.path(out_dir, "adjustment_justification.csv"), row.names = FALSE)
cat("Adjustment-set fit and effective-sample-size cost:\n")
print(results %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)

# --- Word table -------------------------------------------------------------
# formatting: the base row shows dashes where a test against itself is undefined.
fmt_num <- function(x, digits = 1) ifelse(is.na(x), "-", formatC(x, format = "f", digits = digits))
fmt_p   <- function(p) ifelse(is.na(p), "-", formatC(p, format = "g", digits = 2))
fmt_lr  <- function(chisq, df) ifelse(is.na(chisq), "-",
                                      sprintf("%s (%d)", formatC(chisq, format = "f", digits = 1), df))

disp <- data.frame(
  set       = results$label,
  params    = as.character(results$params),
  AIC       = fmt_num(results$AIC, 1),
  dAIC      = fmt_num(results$dAIC, 1),
  LR        = fmt_lr(results$LR_chisq, results$LR_df),
  p         = fmt_p(results$p_value),
  med_ess   = fmt_num(results$median_ess, 1),
  min_ess   = fmt_num(results$min_ess, 1),
  reduction = fmt_num(results$mean_pct_reduction, 1),
  stringsAsFactors = FALSE, check.names = FALSE)

# insert the two family heading rows (age non-linearity, then calendar terms)
heading <- function(text) {
  r <- as.data.frame(as.list(c(text, rep("", ncol(disp) - 1))), stringsAsFactors = FALSE)
  names(r) <- names(disp); r
}
combined <- rbind(disp[1, ],                       # base
                  heading("Assess relationship with age"),
                  disp[2:3, ],                      # age^2, age^3
                  heading("Assess impact of calendar/year adjustment"),
                  disp[4:6, ],                      # season, year, season + year
                  heading("Assess comorbidity coding"),
                  disp[7:8, ])                      # comorbidity 0/1/2/3+, and with season + year

row_base    <- 1
row_head_ag <- 2
rows_age    <- 3:4
row_head_cal<- 5
rows_cal    <- 6:8
row_head_cci<- 9
rows_cci    <- 10:11
head_rows   <- c(row_head_ag, row_head_cal, row_head_cci)

FONT <- 8
ft <- flextable(combined)
ft <- set_header_labels(ft,
                        set = "Adjustment set", params = "Parameters", AIC = "AIC",
                        dAIC = "Difference in AIC", LR = "Likelihood ratio (LR), chi-square (df)",
                        p = "p-value", med_ess = "Median ESS", min_ess = "Minimum ESS",
                        reduction = "Mean reduction in ESS (%)")

# the two family headings span the full width
for (r in head_rows) ft <- merge_at(ft, i = r, j = seq_len(ncol(combined)), part = "body")
ft <- bold(ft, i = head_rows, part = "body")
ft <- bold(ft, i = row_base, j = "set", part = "body")
ft <- italic(ft, i = rows_age, j = "set", part = "body")   # italic + indent the family members
ft <- italic(ft, i = rows_cal, j = "set", part = "body")
ft <- italic(ft, i = rows_cci, j = "set", part = "body")
ft <- bold(ft, part = "header")

ft <- align(ft, part = "all", align = "center")
ft <- align(ft, j = "set", align = "left", part = "all")
ft <- align(ft, i = head_rows, align = "left", part = "body")

box  <- fp_border(color = "black",  width = 1)
soft <- fp_border(color = "grey65", width = 0.5)
ft <- border_outer(ft, border = box, part = "all")     # box around the whole table
ft <- hline(ft, i = row_base,        border = soft, part = "body")
ft <- hline(ft, i = max(rows_age),   border = soft, part = "body")
ft <- hline(ft, i = max(rows_cal),   border = soft, part = "body")
ft <- hline_bottom(ft, border = soft, part = "header")

ft <- fontsize(ft, size = FONT, part = "all")
ft <- padding(ft, padding.top = 1.5, padding.bottom = 1.5, part = "all")
ft <- padding(ft, i = c(rows_age, rows_cal, rows_cci), j = "set", padding.left = 16, part = "body")
ft <- autofit(ft)
ft <- width(ft, j = "set",       width = 2.15)   # +30%
ft <- width(ft, j = "params",    width = 0.85)
ft <- width(ft, j = "AIC",       width = 0.64)   # +20%
ft <- width(ft, j = "dAIC",      width = 0.72)   # +20%
ft <- width(ft, j = "LR",        width = 1.15)
ft <- width(ft, j = "p",         width = 0.60)   # -20%
ft <- width(ft, j = "med_ess",   width = 0.67)
ft <- width(ft, j = "min_ess",   width = 0.67)
ft <- width(ft, j = "reduction", width = 0.74)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Supplementary table. Justification of the primary adjustment set. Starting",
  "from age (linear) and comorbidity (0/1/2+), each candidate is judged by",
  "outcome-model fit (AIC and a likelihood-ratio test against the base) and by",
  "the effective sample size the balancing costs (Kish ESS per hospital). The",
  "adopted primary model is the base plus season and calendar year; the final rows",
  "additionally check a finer comorbidity coding (0/1/2/3+), alone and with season",
  "and calendar year.")
footnote_text <- paste(
  "Age is centred before squaring, so the likelihood-ratio test is identical to",
  "raw polynomials. ESS = (sum of weights)^2 / sum of squared weights per",
  "hospital; the reduction is the mean percent fall from the raw count. A term is",
  "retained only where the fit gain justifies the effective-sample-size cost. The",
  "comorbidity 0/1/2/3+ row splits the base 2+ category; the base is nested within",
  "it, so its test (1 df) asks whether separating 3+ from 2 improves fit.")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "landscape")))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = fp_text(font.size = 9))))
print(doc, target = file.path(out_dir, "adjustment_justification.docx"))

# --- effective sample size scatter for the primary model --------------------
# one point per hospital: effective (y) against true (x) sample size. The dashed
# line is y = x (no loss); a point below it has spent sample size on balancing.
# Aesthetic follows the summary figure (theme_classic, darkblue points).
axis_title_size <- 9
axis_text_size  <- 8
col_base        <- "darkblue"
pt_size         <- 1.6

fit_primary <- readRDS(file.path(out_dir, "fit_primary.rds"))
ess_primary <- hospital_ess(fit_primary, "primary")
lim <- c(0, max(ess_primary$n))

p_ess <- ggplot(ess_primary, aes(n, ess)) +
  geom_abline(slope = 1, intercept = 0, colour = "gray30",
              linewidth = 0.7, linetype = "dashed") +
  geom_hline(yintercept = ess_threshold, colour = "gray70",
             linewidth = 0.4, linetype = "dotted") +
  geom_point(shape = 16, colour = col_base, size = pt_size) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = axis_title_size),
        axis.text  = element_text(size = axis_text_size),
        legend.position = "none") +
  scale_x_continuous("True sample size (patients per hospital)", limits = lim) +
  scale_y_continuous("Effective sample size", limits = lim) +
  coord_equal(xlim = lim, ylim = lim)

ggsave(file.path(out_dir, "ess_scatter_primary.png"), p_ess,
       width = 90, height = 90, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "ess_scatter_primary.pdf"), p_ess,
       width = 90, height = 90, units = "mm")

cat("\nadjustment-set justification table and ESS scatter written.\n")