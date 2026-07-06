# 12  stability heat map: standardised waiting time by hospital and six-month period
# -----------------------------------------------------------------------------
# One value per hospital per six-month period: the age- and comorbidity-
# standardised mean days to decision-to-treat. Hospitals are ordered fastest (top)
# to slowest (bottom) by their overall mean across the whole window.
#
# Inclusion. Only diagnosing hospitals that pass the analysis QC are shown: at
# least min_per_year patients in every year of the window. This is already applied
# in analysis_data.rds; it is repeated here so the figure is self-contained.
#
# Blank cells. A cell is blank only when the hospital had fewer than cell_suppress_n
# patients in that six-month window - too few to show reliably or publish under
# small-number rules. It is not the standardisation, which only shifts a cell's
# value and never removes it. That threshold lives in config.

# ----------------------------  plot appearance  -----------------------------
# every visual choice lives here; nothing further down needs editing to retune it.

fill_low      <- "#2166ac"  # colour for the shortest waits
fill_mid      <- "#f7f7f7"  # colour at the average wait
fill_high     <- "#b2182b"  # colour for the longest waits
fill_cap_lo   <- 0.05       # lower quantile the tile colours are anchored to
fill_cap_hi   <- 0.95       # upper quantile the tile colours are anchored to
legend_lo     <- 10         # lower end of the legend scale (display only)
legend_hi     <- 50         # upper end of the legend scale (display only)
legend_breaks <- seq(legend_lo, legend_hi, by = 10)

blank_fill    <- "grey96"   # base tile colour for a blanked (small-count) cell
hatch_colour  <- "grey75"   # cross-hatch line colour: lighter = softer
hatch_alpha   <- 0.4        # cross-hatch opacity: 1 = solid, lower = fainter
hatch_density <- 0.3       # share of the cell the hatch lines cover
hatch_spacing <- 0.008      # gap between hatch lines: larger = sparser lines
hatch_size    <- 0.1       # hatch line thickness

tile_border   <- "white"    # line drawn between tiles
tile_border_w <- 0.15

axis_text_size    <- 12     # period labels along the x-axis
axis_title_size   <- 14     # axis titles (the y-axis label)
legend_text_size  <- 13     # numbers along the legend scale
legend_title_size <- 14     # legend title
legend_title_gap  <- 10      # space between the legend title and the scale below it
legend_box_gap    <- -2     # space between the plot and the whole legend: smaller = closer

fig_width  <- 7             # inches
fig_height <- 9

# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggpattern)

source("R/01_config.R")

df <- readRDS(file.path(out_dir, "analysis_data.rds"))
cv <- code_covariates(df); df <- cv$data

# diagnosing-site QC (>= min_per_year patients in every window year), matching the
# analysis. Already applied when analysis_data.rds is built, so this is a no-op -
# but it makes the inclusion rule explicit and guarantees only QC-passing
# hospitals are plotted.
wyears <- sort(unique(df$ydiag))
qc_ok <- df %>%
  count(diag_hosp_canon, ydiag) %>%
  complete(diag_hosp_canon, ydiag = wyears, fill = list(n = 0L)) %>%
  group_by(diag_hosp_canon) %>%
  summarise(min_year = min(n), .groups = "drop") %>%
  filter(min_year >= min_per_year) %>%
  pull(diag_hosp_canon)
df <- df %>% filter(diag_hosp_canon %in% qc_ok)

# case-mix expected wait on the analysis cohort: a pooled model on age +
# comorbidity gives each patient an expected wait, and each cell is shifted by how
# far its patients sit from that expectation (indirect standardisation).
pm      <- lm(as.formula(paste("wait ~", paste(c(cv$cont, cv$bin), collapse = " + "))), data = df)
df$pred <- predict(pm)
grand   <- mean(df$wait)

# six-month periods aligned to calendar halves (Jan-Jun, Jul-Dec), not to
# whatever day the data happens to start on. anchor_date is the 1 Jan or 1 Jul
# on or before the earliest diagnosis, so every window is a clean calendar half.
first_date  <- min(df$diagmdy)
anchor_date <- if (month(first_date) <= 6) {
  as.Date(paste0(year(first_date), "-01-01"))
} else {
  as.Date(paste0(year(first_date), "-07-01"))
}
df <- df %>%
  mutate(pnum = floor(as.numeric(interval(anchor_date, diagmdy), "months") / 6) + 1)
pn <- sort(unique(df$pnum))
period_labels <- character(length(pn))
for (i in seq_along(pn)) {
  p <- pn[i]
  s <- anchor_date %m+% months(6 * (p - 1))
  e <- (anchor_date %m+% months(6 * p)) %m-% days(1)
  period_labels[i] <- paste0(format(s, "%b"), "-", format(e, "%b %Y"))
}
df <- df %>% mutate(period = factor(pnum, levels = pn, labels = period_labels))

# one cell per hospital and period: patient count and standardised mean, with the
# cells below the small-count threshold blanked
cell <- df %>%
  group_by(diag_hosp_canon, period) %>%
  summarise(n = n(), std = mean(wait) - mean(pred) + grand, .groups = "drop") %>%
  complete(diag_hosp_canon, period, fill = list(n = 0L)) %>%
  mutate(std = ifelse(n < cell_suppress_n, NA_real_, std))

# order hospitals by their OVERALL standardised mean across all their patients in
# the window (not the mean of the shown cells). Using every patient keeps the rank
# stable and stops a hospital being lifted to the top by one small, fast six-month
# cell while the rest of its record is hatched. slowest first -> fastest at top.
hosp_order <- df %>%
  group_by(diag_hosp_canon) %>%
  summarise(overall = mean(wait) - mean(pred) + grand, .groups = "drop") %>%
  arrange(desc(overall)) %>%
  pull(diag_hosp_canon)
cell <- cell %>% mutate(diag_hosp_canon = factor(diag_hosp_canon, levels = hosp_order))

# diverging scale anchored at the fill_cap quantiles (same three colours as
# before, at the same data positions), but the legend itself is drawn out to
# legend_lo/legend_hi so it reads 0 to 50 - a logical, easy-to-read scale - even
# though the data never reaches those extremes. Beyond the anchored range the
# tile colour is flat (fill_low below lo, fill_high above hi), same as before.
lo <- quantile(cell$std, fill_cap_lo, na.rm = TRUE)
hi <- quantile(cell$std, fill_cap_hi, na.rm = TRUE)

stops     <- c(legend_lo, lo, grand, hi, legend_hi)
stop_cols <- c(fill_low, fill_low, fill_mid, fill_high, fill_high)

data_cells  <- cell %>% filter(!is.na(std))
blank_cells <- cell %>% filter(is.na(std))

p <- ggplot() +
  geom_tile(data = data_cells, aes(period, diag_hosp_canon, fill = std),
            colour = tile_border, linewidth = tile_border_w) +
  scale_fill_gradientn(colours = stop_cols,
                       values = scales::rescale(stops, from = c(legend_lo, legend_hi)),
                       limits = c(legend_lo, legend_hi), breaks = legend_breaks,
                       oob = scales::oob_squish, name = "Mean waiting time\n(days)") +
  geom_tile_pattern(data = blank_cells, aes(period, diag_hosp_canon),
                    fill = blank_fill, colour = tile_border, linewidth = tile_border_w,
                    pattern = "crosshatch", pattern_colour = hatch_colour,
                    pattern_fill = hatch_colour, pattern_alpha = hatch_alpha,
                    pattern_angle = 45, pattern_density = hatch_density,
                    pattern_spacing = hatch_spacing, pattern_size = hatch_size) +
  labs(x = NULL, y = "Hospitals (ordered by mean waiting time)") +
  theme_minimal(base_size = axis_text_size) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size,
                                    margin = margin(b = legend_title_gap)),
        legend.box.spacing = unit(legend_box_gap, "pt"))

# ggpattern draws the cross-hatch via clipped pattern fills, which the default
# pdf() device does not render reliably. cairo_pdf supports it; a PNG is saved
# alongside as a robust fallback.
ggsave(file.path(out_dir, "heatmap_stability.pdf"), p, width = fig_width, height = fig_height,
       device = cairo_pdf)
ggsave(file.path(out_dir, "heatmap_stability.png"), p, width = fig_width, height = fig_height,
       dpi = 300, bg = "white")
# write the cells in plot order (plot_row 1 = top of the figure = fastest) with a
# status flag, so the csv reads top-to-bottom the same way the figure does
plot_row <- data.frame(diag_hosp_canon = rev(hosp_order),
                       plot_row = seq_along(hosp_order))
cell_out <- cell %>%
  left_join(plot_row, by = "diag_hosp_canon") %>%
  mutate(status = case_when(n == 0    ~ "empty",
                            is.na(std) ~ sprintf("blanked (<%d)", cell_suppress_n),
                            TRUE       ~ "shown")) %>%
  arrange(plot_row, period) %>%
  select(plot_row, diag_hosp_canon, period, n, std, status)
write.csv(cell_out, file.path(out_dir, "heatmap_cells.csv"), row.names = FALSE)
cat(sprintf("heat map saved: %d hospitals x %d periods, %d of %d cells blank (<%d patients).\n",
            n_distinct(cell$diag_hosp_canon), nlevels(cell$period),
            sum(is.na(cell$std)), nrow(cell), cell_suppress_n))