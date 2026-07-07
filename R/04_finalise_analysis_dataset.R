# 04  finalise the analysis dataset (provider identity + volume)
# -----------------------------------------------------------------------------
# Applies the crosswalks from 03 to the clinically eligible set from 02:
# canonicalises diagnosing and treating codes to 5-digit sites, keeps patients
# whose diagnosing site is valid and whose treating site does bowel surgery,
# applies the hospital volume threshold, and writes the final analysis dataset,
# the hospital lookup, and the complete flowchart (clinical steps + provider
# steps). If the crosswalks are absent, codes pass through unchanged and only the
# volume cut is applied, so the pipeline still runs before the QC is curated.

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/01_config.R")

df   <- readRDS(clinical_rds)
flow <- read.csv(flow_clinical_csv, stringsAsFactors = FALSE)
if (!all(c("box", "kind", "label", "n", "n_hosp") %in% names(flow)))
  stop(flow_clinical_csv, " is from an older run of 02 - re-run 02 first so it ",
       "writes the boxed flowchart, then re-run 04.")
.prev <- nrow(df)   # box 5 milestone: the clinically eligible cohort
fc_box <- function(box, label, n, n_hosp = NA_integer_) {
  flow  <<- rbind(flow, data.frame(box = box, kind = "box", label = label,
                                   n = as.integer(n), n_hosp = as.integer(n_hosp)))
  .prev <<- as.integer(n)
}
fc_excl <- function(box, reason, n_after) {
  flow  <<- rbind(flow, data.frame(box = box, kind = "exclusion", label = reason,
                                   n = .prev - as.integer(n_after), n_hosp = NA_integer_))
  .prev <<- as.integer(n_after)
}

read_codes <- function(path, col) if (file.exists(path)) trimws(read.csv(path, colClasses = "character")[[col]]) else NULL

diag_x  <- if (file.exists(diag_xwalk_csv))  read.csv(diag_xwalk_csv,  colClasses = "character") else NULL
treat_x <- if (file.exists(treat_xwalk_csv)) read.csv(treat_xwalk_csv, colClasses = "character") else NULL
diag_keep  <- read_codes(diag_include_csv,  "canonical_code")
treat_keep <- read_codes(treat_include_csv, "canonical_code")

# canonicalise to 5-digit sites (identity where no crosswalk row matches) ------
df <- df %>%
  mutate(diag_hosp = trimws(as.character(diag_hosp)),
         tx_site   = trimws(as.character(SITETRET)),
         diag_hosp_canon = canonicalise_hosp(diag_hosp, diag_x),
         tx_hosp_canon   = canonicalise_hosp(tx_site,   treat_x))

# box 6: hospital-level QC (provider identity + volume) ----------------------
# keep valid diagnosing sites (the analysis unit)
if (!is.null(diag_keep)) {
  df <- df %>% filter(diag_hosp_canon %in% diag_keep)
  fc_excl(6, "Patients with invalid hospital site code", nrow(df))
}
# keep patients treated at a bowel-cancer-surgery site
if (!is.null(treat_keep)) {
  df <- df %>% filter(tx_hosp_canon %in% treat_keep)
  fc_excl(6, "Not treated at bowel surgery site", nrow(df))
}

# per-year volume: a site must have at least min_per_year patients in EVERY year
# of the window - diagnosing patients for the diagnosing unit, resections for the
# treating site. The treating cut is applied first (it drops patients), then the
# diagnosing cut is recomputed on what remains.
wyears <- sort(unique(df$ydiag))
keep_by_year <- function(d, unit) {
  d %>% count(unit = .data[[unit]], ydiag) %>%
    tidyr::complete(unit, ydiag = wyears, fill = list(n = 0L)) %>%
    group_by(unit) %>% summarise(min_year = min(n), .groups = "drop") %>%
    filter(min_year >= min_per_year) %>% pull(unit)
}

tx_ok <- keep_by_year(df, "tx_hosp_canon")
df <- df %>% filter(tx_hosp_canon %in% tx_ok)
fc_excl(6, sprintf("Treating site < %d resections per year", min_per_year), nrow(df))

diag_ok <- keep_by_year(df, "diag_hosp_canon")
df <- df %>% filter(diag_hosp_canon %in% diag_ok)
fc_excl(6, sprintf("Diagnosing site < %d patients per year", min_per_year), nrow(df))

fc_box(6, "Analysis cohort", nrow(df), n_distinct(df$diag_hosp_canon))

# total window volume per diagnosing site, for the lookup
df <- df %>% left_join(count(df, diag_hosp_canon, name = "volume"), by = "diag_hosp_canon")

# numeric hospital id for the estimation step (contiguous after the filters) ---
df <- df %>% mutate(hospid = as.integer(factor(diag_hosp_canon)), hosp = hospid)

hosp_lookup <- df %>%
  distinct(hospid, diag_hosp_canon, diag_hosp, diag_hosp_name, canalliance, volume) %>%
  arrange(hospid)

# flowchart for the paper: boxes down the page, exclusions in a side arrow ----
cm <- function(x) formatC(x, format = "d", big.mark = ",")
boxes <- flow[flow$kind == "box", ]
boxes <- boxes[order(boxes$box), ]
cat("\nflowchart:\n")
for (b in boxes$box) {
  ex <- flow[flow$kind == "exclusion" & flow$box == b, ]
  if (nrow(ex)) {
    cat("        |\n        |--- excluded:\n")
    for (i in seq_len(nrow(ex)))
      cat(sprintf("        |       %-46s %7s\n", ex$label[i], cm(ex$n[i])))
    cat("        v\n")
  }
  bx   <- boxes[boxes$box == b, ]
  hosp <- if (!is.na(bx$n_hosp)) sprintf("  (%d diagnosing hospitals)", bx$n_hosp) else ""
  cat(sprintf("  [ Box %d ]  %s\n", b, bx$label))
  cat(sprintf("             N = %s%s\n", cm(bx$n), hosp))
}
cat(sprintf("\nmedian wait %.0f days (IQR %.0f-%.0f)\n",
            median(df$wait), quantile(df$wait, .25), quantile(df$wait, .75)))

saveRDS(df,          file.path(out_dir, "analysis_data.rds"))
saveRDS(hosp_lookup, file.path(out_dir, "hosp_lookup.rds"))

# the structured funnel, plus tidy box and exclusion tables for a drawn figure
write.csv(flow, file.path(out_dir, "analysis_flowchart.csv"), row.names = FALSE)
write.csv(boxes[, c("box", "label", "n", "n_hosp")],
          file.path(out_dir, "flowchart_boxes.csv"), row.names = FALSE)
ex_all <- flow[flow$kind == "exclusion", c("box", "label", "n")]
names(ex_all) <- c("into_box", "reason", "n_excluded")
write.csv(ex_all, file.path(out_dir, "flowchart_exclusions.csv"), row.names = FALSE)

# flowchart figure -----------------------------------------------------------
# The console funnel above, drawn as a figure: milestone boxes down the page,
# each transition's exclusions in a panel to the right. Written as an image
# (PDF and PNG) and, if devEMF is available, as an editable Word file (an EMF
# vector you can ungroup into shapes). Wrapped so a plotting/package problem
# cannot lose the data outputs already saved above.
try({
  bx  <- boxes[order(boxes$box), ]
  nbx <- nrow(bx)
  
  # exclusions feeding each transition are labelled with the box they lead into
  gap_ex <- lapply(seq_len(nbx - 1),
                   function(i) flow[flow$kind == "exclusion" & flow$box == bx$box[i + 1],
                                    c("label", "n")])
  
  # one readable size for every label; nothing bold
  font_sz  <- 3.8
  wrap_box <- 38      # characters per line inside a main box
  wrap_ex  <- 27      # characters per line for an exclusion reason
  line_h   <- 0.42
  
  # main boxes as single sentences with N folded in. Figure wording is kept here
  # so it can be reworded without a 02 re-run; boxes 1/3/6 use the 02 labels (box 1
  # carries the dates). Only the final box names the hospitals.
  disp <- bx$label
  disp[bx$box == 2] <- "Patients with stage 1-3 colon cancer"
  disp[bx$box == 4] <- "Patients that underwent elective colon resection"
  disp[bx$box == 5] <- "Patients with valid waiting time (decision to treat recorded, plausible wait)"
  suffix <- sprintf(" (N~=~%s)", cm(bx$n))
  li <- which.max(bx$box)
  if (!is.na(bx$n_hosp[li]))
    suffix[li] <- sprintf(" (N~=~%s,~%d~diagnosing~hospitals)", cm(bx$n[li]), bx$n_hosp[li])
  # protect the spaces inside the N clause with "~" so wrapping keeps it on one line
  wrap_keep <- function(s) gsub("~", " ", paste(strwrap(s, width = wrap_box), collapse = "\n"))
  titles    <- vapply(paste0(disp, suffix), wrap_keep, character(1), USE.NAMES = FALSE)
  n_tlines <- lengths(strsplit(titles, "\n", fixed = TRUE))
  box_h_i  <- 0.5 + n_tlines * 0.48
  
  # wrap each exclusion reason (its count sits in a right-hand column); the panel
  # needs a line for the "Excluded" heading plus every wrapped reason line
  gap_wrapped <- lapply(gap_ex, function(ex)
    if (nrow(ex)) lapply(ex$label, function(s) strwrap(s, width = wrap_ex)) else list())
  ex_h <- vapply(gap_wrapped, function(wl)
    if (length(wl)) (1 + sum(lengths(wl))) * line_h + 0.3 else 0, numeric(1))
  
  # vertical layout, top down: per-box height, each gap sized to fit its panel
  gap_min <- 1.1
  ytop <- numeric(nbx); ybot <- numeric(nbx); y <- 0
  for (i in seq_len(nbx)) {
    ytop[i] <- y; ybot[i] <- y - box_h_i[i]; y <- ybot[i]
    if (i < nbx) y <- y - max(gap_min, ex_h[i] + 0.5)
  }
  
  box_df <- data.frame(xmin = 1, xmax = 5.4, ymin = ybot, ymax = ytop, title = titles)
  arr_df <- data.frame(x = 3.2, y = ybot[-nbx], yend = ytop[-1])
  
  # exclusion panels: heading, then wrapped reasons on the left with counts right
  pan_df <- data.frame(); lab_df <- data.frame(); num_df <- data.frame(); conn_df <- data.frame()
  for (i in seq_len(nbx - 1)) {
    wl <- gap_wrapped[[i]]
    if (!length(wl)) next
    ex <- gap_ex[[i]]
    g_top <- ybot[i]; g_bot <- ytop[i + 1]; mid <- (g_top + g_bot) / 2
    nlines <- 1 + sum(lengths(wl))
    h <- nlines * line_h + 0.3; p_top <- mid + h / 2; p_bot <- mid - h / 2
    pan_df  <- rbind(pan_df,  data.frame(xmin = 6.6, xmax = 10.0, ymin = p_bot, ymax = p_top))
    conn_df <- rbind(conn_df, data.frame(x = 3.2, y = mid, xend = 6.6, yend = mid))
    lab_df  <- rbind(lab_df, data.frame(x = 6.7, y = p_top - 0.28, label = "Excluded"))
    row <- 1
    for (k in seq_along(wl)) {
      for (m in seq_along(wl[[k]])) {
        yy <- p_top - 0.28 - row * line_h
        lab_df <- rbind(lab_df, data.frame(x = 6.7, y = yy, label = wl[[k]][m]))
        if (m == 1) num_df <- rbind(num_df, data.frame(x = 9.9, y = yy, label = cm(ex$n[k])))
        row <- row + 1
      }
    }
  }
  
  fc <- ggplot() +
    geom_segment(data = conn_df, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey40", linewidth = 0.5,
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    geom_segment(data = arr_df, aes(x = x, y = y, xend = x, yend = yend),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 linewidth = 0.5, colour = "grey40") +
    geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "grey40", linewidth = 0.5) +
    geom_rect(data = pan_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "grey55", linewidth = 0.4) +
    geom_text(data = box_df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = title),
              size = font_sz, lineheight = 0.9) +
    geom_text(data = lab_df, aes(x = x, y = y, label = label),
              hjust = 0, size = font_sz, lineheight = 0.9) +
    geom_text(data = num_df, aes(x = x, y = y, label = label),
              hjust = 1, size = font_sz) +
    coord_cartesian(xlim = c(0.4, 10.1), clip = "off") +
    theme_void()
  
  fig_h <- max(6, (ytop[1] - ybot[nbx]) * 0.42)
  ggsave(file.path(out_dir, "flowchart.pdf"), fc, width = 8.5, height = fig_h)
  ggsave(file.path(out_dir, "flowchart.png"), fc, width = 8.5, height = fig_h,
         dpi = 200, bg = "white")
  cat("flowchart figure written to flowchart.pdf and flowchart.png\n")
  
  # editable Word flowchart. The rvg route is broken on some officer/rvg version
  # pairs, so instead draw the figure to an EMF vector (devEMF) and embed that
  # with officer::body_add_img. Word shows it as a vector picture; right-click >
  # Group > Ungroup turns the boxes, arrows and labels into editable shapes. The
  # standalone .emf is also left in place. Best-effort; never stops the pipeline.
  if (requireNamespace("devEMF", quietly = TRUE)) {
    emf_path <- file.path(out_dir, "flowchart.emf")
    drawn <- tryCatch({
      devEMF::emf(file = emf_path, width = 8.5, height = fig_h)
      print(fc); dev.off(); TRUE
    }, error = function(e) { try(dev.off(), silent = TRUE)
      message("EMF export skipped: ", conditionMessage(e)); FALSE })
    if (drawn) {
      cat("editable vector written to flowchart.emf\n")
      if (requireNamespace("officer", quietly = TRUE)) {
        tryCatch({
          library(officer)
          # scale the EMF to fit a portrait page, keeping its aspect ratio
          s   <- min(7.0 / 8.5, 9.2 / fig_h)
          doc <- read_docx()
          doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "portrait")))
          doc <- body_add_img(doc, src = emf_path, width = 8.5 * s, height = fig_h * s)
          print(doc, target = file.path(out_dir, "flowchart.docx"))
          cat("editable flowchart written to flowchart.docx",
              "(right-click the figure > Group > Ungroup to edit the shapes)\n")
        }, error = function(e) message("Word embed skipped: ", conditionMessage(e)))
      }
    }
  } else {
    cat("devEMF not installed - editable flowchart skipped (PDF and PNG written).",
        "Run install.packages('devEMF') to enable the editable Word/EMF output.\n")
  }
})