# heatmap_diagnostics  how sparse is the stability grid, and why
# -----------------------------------------------------------------------------
# Numerical companion to 12_heatmap_stability.R. It answers, in numbers, the
# question the grey cells raise: how much of the hospital x period grid can be
# shown, which hospitals and periods are worst, how the levers (wider periods,
# higher-volume subset) change that, and what range of standardised values the
# colour scale has to span. Sparsity is about counts, so most of this needs only
# patient counts; the value distribution uses the same case-mix model as the plot.
#
# All console output is also written to heatmap_diag_report.txt, and the final
# block is a compact plot-design summary meant to be copied out whole.

library(dplyr)
library(tidyr)
library(lubridate)

source("R/01_config.R")
min_cell_n <- cell_suppress_n

df <- readRDS(file.path(out_dir, "analysis_data.rds"))
start_date <- min(df$diagmdy)

# tee everything below to a text file as well as the console
report_file <- file.path(out_dir, "heatmap_diag_report.txt")
while (sink.number() > 0) sink()          # clear any sink left open by a failed run
sink(report_file, split = TRUE)

# the cell grid for a given period width: one row per hospital per period, with
# the patient count and a three-way status. Absent hospital-periods become n = 0.
cell_grid <- function(months) {
  d <- df %>%
    mutate(bin = floor(as.numeric(interval(start_date, diagmdy), "months") / months) + 1)
  labels <- sort(unique(d$bin))
  d$bin <- factor(d$bin, levels = labels)
  counts <- d %>% count(hosp, diag_hosp_canon, bin, name = "n")
  hospitals <- d %>% distinct(hosp, diag_hosp_canon)
  crossing(hospitals, bin = factor(labels, levels = labels)) %>%
    left_join(counts, by = c("hosp", "diag_hosp_canon", "bin")) %>%
    mutate(n = replace_na(n, 0),
           status = case_when(n == 0        ~ "empty",
                              n < min_cell_n ~ "suppressed",
                              TRUE           ~ "shown"))
}

# --- 1. period-width sweep: does a wider period buy back coverage? ------------
sweep_widths <- c(4, 6, 12)
period_sweep <- data.frame()
for (w in sweep_widths) {
  cell <- cell_grid(w)
  n_periods <- length(levels(cell$bin))
  per_hosp  <- cell %>% group_by(hosp) %>% summarise(shown = sum(status == "shown"), .groups = "drop")
  period_sweep <- rbind(period_sweep, data.frame(
    period_months   = w,
    n_periods       = n_periods,
    cells           = nrow(cell),
    pct_shown       = round(100 * mean(cell$status == "shown"), 1),
    pct_suppressed  = round(100 * mean(cell$status == "suppressed"), 1),
    pct_empty       = round(100 * mean(cell$status == "empty"), 1),
    pct_grey        = round(100 * mean(cell$status != "shown"), 1),
    hosp_fully_shown= sum(per_hosp$shown == n_periods),
    hosp_none_shown = sum(per_hosp$shown == 0),
    median_shown    = median(per_hosp$shown)))
}
write.csv(period_sweep, file.path(out_dir, "heatmap_diag_period_sweep.csv"), row.names = FALSE)
cat("Coverage by period width (grey = any suppressed or empty cell):\n")
print(period_sweep, row.names = FALSE)

# --- everything below uses the six-month grid the plot defaults to -----------
cell <- cell_grid(6)
n_periods <- length(levels(cell$bin))

# --- 2. how the plot looks: overall greyness and greyness by period ----------
cat(sprintf("\nSix-month grid: %d hospitals x %d periods = %d cells.\n",
            n_distinct(cell$hosp), n_periods, nrow(cell)))
cat(sprintf("Shown %.1f%%, suppressed %.1f%%, empty %.1f%%  ->  %.1f%% of the plot is grey.\n",
            100 * mean(cell$status == "shown"), 100 * mean(cell$status == "suppressed"),
            100 * mean(cell$status == "empty"), 100 * mean(cell$status != "shown")))

by_period <- cell %>% group_by(bin) %>%
  summarise(shown = sum(status == "shown"),
            pct_grey = round(100 * mean(status != "shown"), 1), .groups = "drop")
cat("\nGreyness by period (window edges are often worst):\n")
print(as.data.frame(by_period), row.names = FALSE)

# --- 3. distribution of shown periods per hospital ---------------------------
per_hosp <- cell %>%
  group_by(hosp, diag_hosp_canon) %>%
  summarise(total_n           = sum(n),
            periods_shown     = sum(status == "shown"),
            periods_suppressed= sum(status == "suppressed"),
            periods_empty     = sum(status == "empty"),
            .groups = "drop")
write.csv(per_hosp %>% arrange(periods_shown, total_n),
          file.path(out_dir, "heatmap_diag_hospital.csv"), row.names = FALSE)

shown_dist <- per_hosp %>% count(periods_shown, name = "hospitals") %>%
  mutate(pct_of_hospitals = round(100 * hospitals / sum(hospitals), 1))
cat("\nHospitals by number of periods shown (0 = never appears on the plot):\n")
print(as.data.frame(shown_dist), row.names = FALSE)

# --- 4. it is a volume story: shown periods against hospital volume ----------
vol_band <- per_hosp %>%
  mutate(volume_band = cut(total_n,
                           breaks = c(0, 40, 60, 100, 150, Inf),
                           labels = c("<40", "40-59", "60-99", "100-149", "150+"))) %>%
  group_by(volume_band) %>%
  summarise(hospitals          = n(),
            median_total_n     = median(total_n),
            mean_periods_shown = round(mean(periods_shown), 2),
            pct_fully_shown    = round(100 * mean(periods_shown == n_periods), 0),
            .groups = "drop")
cat("\nPeriods shown by hospital volume band (the sparseness is concentrated in small sites):\n")
print(as.data.frame(vol_band), row.names = FALSE)

# --- 5. raw cell-count distribution ------------------------------------------
count_dist <- cell %>%
  mutate(band = cut(n, breaks = c(-1, 0, min_cell_n - 1, Inf),
                    labels = c("0 (empty)", sprintf("1-%d (suppressed)", min_cell_n - 1),
                               sprintf(">=%d (shown)", min_cell_n)))) %>%
  count(band, name = "cells") %>%
  mutate(pct = round(100 * cells / sum(cells), 1))
cat("\nCell patient-count distribution:\n")
print(as.data.frame(count_dist), row.names = FALSE)

# --- 6. standardised-value distribution: what the colour scale must span ------
# same case-mix model as the plot (age + cci), so the numbers match the fill.
cv <- code_covariates(df); dm <- cv$data
pm <- lm(as.formula(paste("wait ~", paste(c(cv$cont, cv$bin), collapse = " + "))), data = dm)
dm$pred <- predict(pm)
grand   <- mean(dm$wait)
start2  <- min(dm$diagmdy)
dm <- dm %>% mutate(bin = floor(as.numeric(interval(start2, diagmdy), "months") / 6) + 1)
std_shown <- dm %>%
  group_by(hosp, bin) %>%
  summarise(n = n(), std = mean(wait) - mean(pred) + grand, .groups = "drop") %>%
  filter(n >= min_cell_n) %>% pull(std)
qs <- quantile(std_shown, c(0, .01, .05, .25, .5, .75, .95, .99, 1))
cat(sprintf("\nStandardised days across the %d shown cells (grand mean %.1f):\n", length(std_shown), grand))
cat(sprintf("  min %.1f  p1 %.1f  p5 %.1f  p25 %.1f  median %.1f  p75 %.1f  p95 %.1f  p99 %.1f  max %.1f\n",
            qs[1], qs[2], qs[3], qs[4], qs[5], qs[6], qs[7], qs[8], qs[9]))
cat(sprintf("  cells beyond p1-p99 (candidates for a capped scale): %d (%.1f%%)\n",
            sum(std_shown < qs[2] | std_shown > qs[8]), 100 * mean(std_shown < qs[2] | std_shown > qs[8])))

# --- 7. volume-floor sweep: cost of a higher-volume subset -------------------
floors <- c(40, 60, 80, 100, 120)
volume_sweep <- data.frame()
for (f in floors) {
  keep <- per_hosp$total_n >= f
  volume_sweep <- rbind(volume_sweep, data.frame(
    min_total_volume = f,
    hospitals_kept   = sum(keep),
    pct_hospitals    = round(100 * mean(keep), 0),
    pct_patients     = round(100 * sum(per_hosp$total_n[keep]) / sum(per_hosp$total_n), 0),
    mean_periods_shown_kept = round(mean(per_hosp$periods_shown[keep]), 2)))
}
write.csv(volume_sweep, file.path(out_dir, "heatmap_diag_volume_sweep.csv"), row.names = FALSE)
cat("\nRestricting to higher-volume hospitals (trade coverage against how many hospitals remain):\n")
print(volume_sweep, row.names = FALSE)

# --- consolidated block, meant to be copied whole ----------------------------
g6  <- period_sweep$pct_grey[period_sweep$period_months == 6]
g12 <- period_sweep$pct_grey[period_sweep$period_months == 12]
cat("\n=================  PLOT DESIGN SUMMARY (copy this block)  =================\n")
cat(sprintf("hospitals (rows)                : %d\n", nrow(per_hosp)))
cat(sprintf("periods                         : %d at 6 months, %d at 12 months\n",
            period_sweep$n_periods[period_sweep$period_months == 6],
            period_sweep$n_periods[period_sweep$period_months == 12]))
cat(sprintf("grey share of grid              : %.1f%% at 6 months, %.1f%% at 12 months\n", g6, g12))
cat(sprintf("hospitals fully shown (6 month) : %d of %d\n",
            period_sweep$hosp_fully_shown[period_sweep$period_months == 6], nrow(per_hosp)))
cat(sprintf("hospitals never shown (6 month) : %d\n",
            period_sweep$hosp_none_shown[period_sweep$period_months == 6]))
cat(sprintf("std days median (shown cells)   : %.1f\n", qs[5]))
cat(sprintf("std days p5-p95 (suggested fill): %.1f to %.1f\n", qs[3], qs[7]))
cat(sprintf("std days full range             : %.1f to %.1f\n", qs[1], qs[9]))
cat(sprintf("hospitals kept at >=60 patients : %d (%.0f%% of hospitals, %.0f%% of patients)\n",
            volume_sweep$hospitals_kept[volume_sweep$min_total_volume == 60],
            volume_sweep$pct_hospitals[volume_sweep$min_total_volume == 60],
            volume_sweep$pct_patients[volume_sweep$min_total_volume == 60]))
cat("==========================================================================\n")

cat("\nfull report written to heatmap_diag_report.txt; CSVs: period_sweep, hospital, volume_sweep\n")
sink()