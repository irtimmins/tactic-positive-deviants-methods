# 05  Table 1: cohort characteristics with waiting time, written to Word
# -----------------------------------------------------------------------------
# One row per characteristic level: number and percent of patients (or mean (SD)
# for a continuous row), plus the mean (SD) waiting time in that group. Written
# to Word with flextable (needs flextable and officer). Also a short
# representativeness check by cancer alliance.

library(dplyr)
library(flextable)
library(officer)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

tww <- "Urgent suspected cancer referral\n(two week wait)"
route_order <- c(tww, "GP referral", "Screening", "Other outpatient", "Inpatient elective")

d <- df %>% mutate(
  age_grp = cut(age, c(-Inf, 50, 60, 70, 80, Inf),
                labels = c("<50", "50-59", "60-69", "70-79", "80+")),
  sex     = factor(sexf, levels = c("Male", "Female")),
  stage_f = factor(as.character(stage), levels = c("1", "2", "3")),
  eth_f   = factor(as.character(eth),
                   levels = c("White","Asian","Black","Mixed","Other","Unknown")),
  imd_f   = factor(imd),
  cci_f   = factor(ifelse(cci_n_conditions >= 2, "2+", as.character(cci_n_conditions)),
                   levels = c("0", "1", "2+")),
  route_f  = factor(recode(as.character(route), TWW = tww)),
  season_f = factor(case_when(q2 == 1 ~ "Apr-Jun",
                              q3 == 1 ~ "Jul-Sep",
                              q4 == 1 ~ "Oct-Dec",
                              TRUE    ~ "Jan-Mar"),
                    levels = c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")),
  year_f   = factor(ydiag)
)
d$route_f <- factor(d$route_f, levels = c(route_order, setdiff(levels(d$route_f), route_order)))

# one block per characteristic: a heading row (the label, no data alongside),
# then one indented row per level with n (%) and mean (SD) wait.
cat_block <- function(d, var, label) {
  levels_tab <- d %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(Level = .data[[var]]) %>%
    summarise(np = n(), wm = mean(wait), wsd = sd(wait), .groups = "drop") %>%
    transmute(item = as.character(Level),
              patients = sprintf("%s (%.1f)", formatC(np, format = "d", big.mark = ","),
                                 100 * np / nrow(d)),
              wait = sprintf("%.1f (%.1f)", wm, wsd),
              is_head = FALSE, indent = TRUE)
  bind_rows(tibble(item = label, patients = "", wait = "", is_head = TRUE, indent = FALSE),
            levels_tab)
}


cat(sprintf("Age at diagnosis: mean %.1f years (SD %.1f)\n", mean(d$age), sd(d$age)))

tab1 <- bind_rows(
  tibble(item = "Patients, total",
         patients = formatC(nrow(d), format = "d", big.mark = ","),
         wait = sprintf("%.1f (%.1f)", mean(d$wait), sd(d$wait)),
         is_head = FALSE, indent = FALSE),
  cat_block(d, "age_grp",  "Age group"),
  cat_block(d, "sex",      "Sex"),
  cat_block(d, "stage_f",  "Stage at diagnosis"),
  cat_block(d, "eth_f",    "Ethnicity"),
  cat_block(d, "imd_f",    "Deprivation quintile"),
  cat_block(d, "cci_f",    "Charlson comorbidity index"),
  cat_block(d, "route_f",  "Route to diagnosis"),
  cat_block(d, "year_f",   "Calendar year"),
  cat_block(d, "season_f", "Season of diagnosis")
)

head_rows   <- which(tab1$is_head)                       # characteristic headings
indent_rows <- which(tab1$indent)                        # the levels under each heading
total_row   <- which(tab1$item == "Patients, total")

disp <- tab1[, c("item", "patients", "wait")]
ft <- flextable(disp)
ft <- set_header_labels(ft, item = "Characteristic",
                        patients = "Patients, n (%)", wait = "Waiting time, days\nmean (SD)")
ft <- bold(ft, part = "header")
ft <- bold(ft, j = "item", i = c(head_rows, total_row), part = "body")   # headings and the total
ft <- align(ft, j = "item", align = "left", part = "all")
ft <- align(ft, j = c("patients", "wait"), align = "center", part = "all")
ft <- valign(ft, valign = "top", part = "body")
ft <- add_footer_lines(ft, paste("Patients as n (%); age at diagnosis as mean (SD).",
                                 "Waiting time is mean (SD) days from diagnosis to decision-to-treat."))
ft <- fontsize(ft, size = 8, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
ft <- padding(ft, j = "item", padding.left = 3, part = "body")                     # headings flush
ft <- padding(ft, i = indent_rows, j = "item", padding.left = 18, part = "body")   # levels indented
ft <- line_spacing(ft, space = 1.3, part = "all")

# bordering to match the other tables: a black box round the whole table, a soft
# header underline, and very faint grey rules between the rows.
ft <- border_remove(ft)
box   <- fp_border(color = "black",  width = 1)
soft  <- fp_border(color = "grey65", width = 0.5)
faint <- fp_border(color = "grey85", width = 0.5)
ft <- border_outer(ft, border = box, part = "all")
ft <- hline_bottom(ft, border = soft, part = "header")
ft <- hline(ft, i = seq_len(nrow(disp) - 1), border = faint, part = "body")

ft <- autofit(ft)
ft <- width(ft, j = "item",     width = 1.92)   # -20%
ft <- width(ft, j = "patients", width = 1.1)
ft <- width(ft, j = "wait",     width = 1.2)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

save_as_docx(ft, path = file.path(out_dir, "table1_characteristics.docx"))
write.csv(disp, file.path(out_dir, "table1_characteristics.csv"), row.names = FALSE)
cat("Table 1 written to Output/table1_characteristics.docx\n")

# representativeness: patient share by cancer alliance -------------------------
alliance <- df %>%
  count(canalliance, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))
print(as.data.frame(alliance), row.names = FALSE)
write.csv(alliance, file.path(out_dir, "table1_alliance.csv"), row.names = FALSE)