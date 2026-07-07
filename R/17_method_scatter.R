# 17  method-comparison scatter panels
# -----------------------------------------------------------------------------
# The main model's shrunk posterior mean (y-axis) against the other generators,
# one figure per estimand (sustained, improvement), six panels each:
#   (a) regression-based direct standardisation, shrunk
#   (b) no adjustment (raw mean), shrunk
#   (c) funnel plot, indirectly standardised mean (sustained) / change (improvement)
#   (d) classic funnel z-score
#   (e) funnel z-score with over-dispersion correction
#   (f) funnel z-score with Winsorised over-dispersion correction
# Panels a-c are on the day scale, so the solid diagonal is equality (the indirect
# estimate in c is a different standardisation, shown on the same scale for
# reference). Panels d-f are funnel z-scores, a different metric, so those show
# association rather than agreement. Each panel carries Spearman's rho. Equal-size
# panels on a 2-column grid (a b / c d / e f), axes shared down each column.

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

source("R/01_config.R")

# --- tunable style -----------------------------------------------------------
axis_title_size <- 8
axis_text_size  <- 8
axis_lineheight <- 0.9
pt_size         <- 1.6
pt_alpha        <- 0.7
rho_size        <- 3.0
col_base        <- "darkblue"

lab_size        <- 14      # a-f tag size
lab_x           <- 0.02
lab_y           <- 1.00
lab_top_margin  <- 16

panel_theme <- theme_classic(base_size = 11) +
  theme(axis.title      = element_text(size = axis_title_size, lineheight = axis_lineheight),
        axis.text       = element_text(size = axis_text_size),
        legend.position = "none",
        plot.margin     = margin(t = lab_top_margin, r = 4, b = 2, l = 2))

rho_label <- function(x, y) {
  r <- cor(x, y, method = "spearman", use = "complete.obs")
  sprintf("rho = %.2f", r)
}

# same-units panel: shrunk mean (y) vs another day-scale estimate (x), a solid
# y = x line so agreement is the diagonal. Fixed square limits shared by the figure.
same_units_panel <- function(d, xvar, xlab, ylab, lim) {
  dd <- data.frame(x = d[[xvar]], y = d$main_shrunk_mean)
  ggplot(dd, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, colour = "gray30", linewidth = 0.7) +
    geom_point(shape = 16, colour = col_base, size = pt_size, alpha = pt_alpha) +
    annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1,
             label = rho_label(dd$x, dd$y), size = rho_size, colour = "gray30") +
    panel_theme +
    scale_x_continuous(xlab, breaks = breaks_width(10)) +
    scale_y_continuous(ylab, breaks = breaks_width(10)) +
    coord_cartesian(xlim = lim, ylim = lim)
}

# different-units panel: shrunk mean (y, days) vs a funnel z-score (x), free x, y
# fixed to the shared range, a dashed line at z = 0 (the funnel target). No diagonal.
diff_units_panel <- function(d, xvar, xlab, ylab, lim) {
  dd <- data.frame(x = d[[xvar]], y = d$main_shrunk_mean)
  xr <- range(dd$x, na.rm = TRUE)
  ggplot(dd, aes(x, y)) +
    geom_vline(xintercept = 0, colour = "gray70", linewidth = 0.5, linetype = "dashed") +
    geom_point(shape = 16, colour = col_base, size = pt_size, alpha = pt_alpha) +
    annotate("text", x = xr[1], y = lim[2], hjust = 0, vjust = 1,
             label = rho_label(dd$x, dd$y), size = rho_size, colour = "gray30") +
    panel_theme +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab, breaks = breaks_width(10)) +
    coord_cartesian(ylim = lim)
}

# six equal panels: a b / c d / e f.
build_figure <- function(d, ylab, x_reg, x_raw, x_ind, lim = NULL) {
  if (is.null(lim)) {
    rng <- range(c(d$main_shrunk_mean, d$reg_shrunk_mean, d$raw_shrunk_mean, d$indirect),
                 na.rm = TRUE)
    pad <- 0.04 * diff(rng)
    lim <- c(rng[1] - pad, rng[2] + pad)
  }
  pa <- same_units_panel(d, "reg_shrunk_mean", x_reg, ylab, lim)
  pb <- same_units_panel(d, "raw_shrunk_mean", x_raw, ylab, lim)
  pc <- same_units_panel(d, "indirect",        x_ind, ylab, lim)
  pd <- diff_units_panel(d, "funnel_z",        "Funnel z-score (classic)", ylab, lim)
  pe <- diff_units_panel(d, "funnel_z_re",     "Funnel z-score (over-dispersion)", ylab, lim)
  pf <- diff_units_panel(d, "funnel_z_re_rob",
                         "Funnel z-score\n(Winsorised over-dispersion)", ylab, lim)
  
  design <- "AB\nCD\nEF"
  pa + pb + pc + pd + pe + pf +
    plot_layout(design = design) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag          = element_text(size = lab_size, face = "bold"),
          plot.tag.position = c(lab_x, lab_y))
}

# --- data --------------------------------------------------------------------
sus <- read.csv(file.path(out_dir, "method_scatter_sustained.csv"), stringsAsFactors = FALSE)
imp <- read.csv(file.path(out_dir, "method_scatter_improve.csv"),   stringsAsFactors = FALSE)

fig_sus <- build_figure(
  sus,
  ylab  = "mean waiting time\n(balancer weighted direct standardisation\nwith shrinkage)",
  x_reg = "mean waiting time\n(regression-based direct standardisation\nwith shrinkage)",
  x_raw = "mean waiting time\n(no adjustment but still with shrinkage)",
  x_ind = "mean waiting time\n(funnel plot-based indirect standardisation)",
  lim   = c(10, 50))

fig_imp <- build_figure(
  imp,
  ylab  = "change in mean waiting time\n(balancer weighted direct standardisation\nwith shrinkage)",
  x_reg = "change in mean waiting time\n(regression-based direct standardisation\nwith shrinkage)",
  x_raw = "change in mean waiting time\n(no adjustment but still with shrinkage)",
  x_ind = "change in mean waiting time\n(funnel plot-based indirect standardisation)",
  lim   = NULL)

fig_w <- 180; fig_h <- 230    # physical size, same for both figures and both formats
png_dpi <- 300                # slightly lower-resolution png; the figure size is unchanged
ggsave(file.path(out_dir, "figure_method_scatter_sustained.png"), fig_sus,
       width = fig_w, height = fig_h, units = "mm", dpi = png_dpi, bg = "white")
ggsave(file.path(out_dir, "figure_method_scatter_sustained.pdf"), fig_sus,
       width = fig_w, height = fig_h, units = "mm", bg = "white")
ggsave(file.path(out_dir, "figure_method_scatter_improve.png"), fig_imp,
       width = fig_w, height = fig_h, units = "mm", dpi = png_dpi, bg = "white")
ggsave(file.path(out_dir, "figure_method_scatter_improve.pdf"), fig_imp,
       width = fig_w, height = fig_h, units = "mm", bg = "white")

cat("method-comparison scatter figures saved (sustained, improvement; panels a-f).\n")