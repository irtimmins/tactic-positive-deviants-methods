# 13  summary figure: four panels combined with patchwork (a-d)
# -----------------------------------------------------------------------------
#   (a) caterpillar of sustained performance (y from 0 days)
#   (b) caterpillar of improvement (change over the period)
#   (c) scatter: CCI 0-1 stratum mean (x) vs the overall mean (y)
#   (d) scatter: CCI 2+  stratum mean (x) vs the overall mean (y)
# The two scatters carry credible intervals in BOTH directions and a dashed y=x
# line, so a hospital off the line performs differently in that stratum than
# overall. Aesthetic (theme_classic, grey reference line, segment CIs, palette)
# follows the house style; top performers are not highlighted yet.

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

source("R/01_config.R")

# --- tunable style (adjust here) --------------------------------------------
axis_title_size <- 10       # axis titles
axis_text_size  <- 9       # axis tick labels
axis_lineheight <- 0.9     # line spacing within the two-line c/d axis titles

cat_ci_alpha <- 0.30       # caterpillar CI opacity (a, b) - lighter
cat_ci_lwd   <- 0.25
cat_pt_size  <- 0.8
cat_xtitle_gap   <- 2.5    # pt gap from the x-axis line to the "Hospitals" title
a_ytitle_vjust   <- 0.4      # panel a y-title position within its strip: raise to move
b_ytitle_vjust   <- 0      # panel b y-title position (tuned separately from a, since
# the two panels' tick labels differ in width). Raise to
# move the title TOWARD the axis, lower to move it away.

sc_ci_alpha  <- 0.15       # scatter CI opacity (c, d) - lighter
sc_ci_lwd    <- 0.30
sc_pt_size   <- 1.6        # scatter dot size (bigger than the caterpillars)
sc_pt_alpha <- 0.7

# panel labels a-d (patchwork tags): position is npc within each panel including
# its top margin, so lab_y near 1 places the tag in the top margin, above the plot
lab_size       <- 14
lab_x          <- 0.02     # 0 = flush left; larger = further right
lab_y          <- 1.00     # 1 = top; lower to bring the tag down into the panel
lab_top_margin <- 16       # pt of space above each panel for the tag to sit in

# house palette (kept for when top performers are highlighted later)
col_base    <- "darkblue"
col_high    <- "darkorange3"   # overall shortest waits
col_improve <- "#009E73"       # greatest improvement

read_ranks <- function(file) {
  if (!file.exists(file)) return(NULL)
  read.csv(file, stringsAsFactors = FALSE) %>%
    transmute(diag_hosp = as.character(diag_hosp), post_mean, ci_lo, ci_hi, exp_rank)
}

# shared theme so all four panels match; top margin leaves room for the a-d label
panel_theme <- theme_classic(base_size = 11) +
  theme(axis.title  = element_text(size = axis_title_size, lineheight = axis_lineheight),
        axis.text   = element_text(size = axis_text_size),
        legend.position = "none",
        plot.margin = margin(t = lab_top_margin, r = 4, b = 2, l = 2))

# --- caterpillar panel -------------------------------------------------------
# hospitals ordered by expected posterior rank; segment = credible interval; a
# grey reference line at the given level. from_zero floors the y-axis at 0.
caterpillar_panel <- function(d, ylab, ref_line, from_zero = FALSE, ytitle_vjust = 1) {
  d <- d %>% arrange(exp_rank) %>% mutate(rank = row_number())
  p <- ggplot(d, aes(rank, post_mean)) +
    geom_hline(yintercept = ref_line, linewidth = 1, alpha = 0.6, colour = "gray30") +
    geom_segment(aes(xend = rank, y = ci_lo, yend = ci_hi),
                 colour = col_base, alpha = cat_ci_alpha, linewidth = cat_ci_lwd) +
    geom_point(shape = 16, colour = col_base, size = cat_pt_size) +
    panel_theme +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          axis.title.x = element_text(size = axis_title_size, margin = margin(t = cat_xtitle_gap)),
          axis.title.y = element_text(size = axis_title_size, lineheight = axis_lineheight,
                                      vjust = ytitle_vjust)) +
    scale_x_continuous("Hospitals") +
    scale_y_continuous(ylab, breaks = breaks_width(10))
  if (from_zero) p <- p + coord_cartesian(ylim = c(0, max(d$ci_hi) * 1.02))
  p
}

# --- scatter panel -----------------------------------------------------------
# stratum estimate (x) vs overall estimate (y), credible intervals both ways,
# equal aspect and a dashed y=x line so agreement is a 45-degree diagonal.
scatter_panel <- function(overall, stratum, xlab, ylab) {
  d <- overall %>%
    select(diag_hosp, y = post_mean, y_lo = ci_lo, y_hi = ci_hi) %>%
    inner_join(stratum %>% select(diag_hosp, x = post_mean, x_lo = ci_lo, x_hi = ci_hi),
               by = "diag_hosp")
  lim <- range(c(d$x_lo, d$x_hi, d$y_lo, d$y_hi), na.rm = TRUE)
  ggplot(d, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, colour = "gray30",
                linewidth = 0.7, linetype = "solid") +
    geom_segment(aes(x = x_lo, xend = x_hi, y = y, yend = y),      # horizontal CI (stratum)
                 colour = col_base, alpha = sc_ci_alpha, linewidth = sc_ci_lwd) +
    geom_segment(aes(x = x, xend = x, y = y_lo, yend = y_hi),      # vertical CI (overall)
                 colour = col_base, alpha = sc_ci_alpha, linewidth = sc_ci_lwd) +
    geom_point(shape = 16, colour = col_base, size = sc_pt_size, alpha = sc_pt_alpha) +
    panel_theme +
    scale_x_continuous(xlab, breaks = breaks_width(10)) +
    scale_y_continuous(ylab, breaks = breaks_width(10)) +
    coord_equal(xlim = lim, ylim = lim)
}

# --- data --------------------------------------------------------------------
sus <- read_ranks(file.path(out_dir, "ranks_sustained.csv"))
imp <- read_ranks(file.path(out_dir, "ranks_improve.csv"))
s01 <- read_ranks(file.path(out_dir, "ranks_strata_01.csv"))
s2  <- read_ranks(file.path(out_dir, "ranks_strata_2.csv"))

# --- panels ------------------------------------------------------------------
panel_a <- caterpillar_panel(sus, "Mean waiting time (days)",
                             ref_line = mean(sus$post_mean), from_zero = TRUE,
                             ytitle_vjust = a_ytitle_vjust)
panel_b <- caterpillar_panel(imp, "Change in waiting time (days)",
                             ref_line = 0, from_zero = FALSE,
                             ytitle_vjust = b_ytitle_vjust)
panel_c <- scatter_panel(sus, s01,
                         "Patients with 0 or 1 comorbidities,\nMean waiting time (days), per hospital",
                         "All patients,\nMean waiting time (days), per hospital")
panel_d <- scatter_panel(sus, s2,
                         "Patients with 2+ comorbidities,\nMean waiting time (days), per hospital",
                         "All patients,\nMean waiting time (days), per hospital")

# patchwork aligns the panel axes across the grid. The shared y-title strip is
# sized to the widest title (the two-line c/d titles), so a/b's title cannot be
# pulled in with its own margin - use a_ytitle_vjust / b_ytitle_vjust to
# reposition each title within that strip instead.
figure <- panel_a + panel_b + panel_c + panel_d +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = lab_size, face = "bold"),
        plot.tag.position = c(lab_x, lab_y))

ggsave(file.path(out_dir, "figure_summary_panels.png"), figure,
       width = 180, height = 178, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "figure_summary_panels.pdf"), figure,
       width = 180, height = 178, units = "mm", bg = "white")

cat("summary figure saved (panels a-d).\n")