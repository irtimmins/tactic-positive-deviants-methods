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
  fc_excl(6, "Diagnosing hospital not resolvable to one site", nrow(df))
}
# keep patients treated at a bowel-cancer-surgery site
if (!is.null(treat_keep)) {
  df <- df %>% filter(tx_hosp_canon %in% treat_keep)
  fc_excl(6, "Treated at a non-bowel-surgery site", nrow(df))
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
fc_excl(6, sprintf("Treating site < %d resections/year", min_per_year), nrow(df))

diag_ok <- keep_by_year(df, "diag_hosp_canon")
df <- df %>% filter(diag_hosp_canon %in% diag_ok)
fc_excl(6, sprintf("Diagnosing site < %d patients/year", min_per_year), nrow(df))

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