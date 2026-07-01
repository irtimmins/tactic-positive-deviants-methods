# 05  Table 1: cohort characteristics with waiting time, written to Word
# -----------------------------------------------------------------------------
# One row per characteristic level: number and percent of patients (or mean (SD)
# for a continuous row), plus the mean (SD) waiting time in that group. Written
# to Word with flextable (needs flextable and officer). Also a short
# representativeness check by cancer alliance.

library(dplyr)
library(flextable)

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
  cci_f   = factor(ifelse(cci_n_conditions >= 3, "3+", as.character(cci_n_conditions)),
                   levels = c("0", "1", "2", "3+")),
  route_f = factor(recode(as.character(route), TWW = tww)),
  year_f  = factor(ydiag)
)
d$route_f <- factor(d$route_f, levels = c(route_order, setdiff(levels(d$route_f), route_order)))

# a categorical block: one row per level with n, %, and mean (SD) wait
cat_block <- function(d, var, label) {
  d %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(Level = .data[[var]]) %>%
    summarise(np = n(), wm = mean(wait), wsd = sd(wait), .groups = "drop") %>%
    transmute(Characteristic = label, Level = as.character(Level),
              patients = sprintf("%s (%.1f)", formatC(np, format = "d", big.mark = ","),
                                 100 * np / nrow(d)),
              wait = sprintf("%.1f (%.1f)", wm, wsd))
}


cat(sprintf("Age at diagnosis: mean %.1f years (SD %.1f)\n", mean(d$age), sd(d$age)))

tab1 <- bind_rows(
  tibble(Characteristic = "Patients, total", Level = "",
         patients = formatC(nrow(d), format = "d", big.mark = ","),
         wait = sprintf("%.1f (%.1f)", mean(d$wait), sd(d$wait))),
  cat_block(d, "age_grp", "Age group"),
  cat_block(d, "sex",     "Sex"),
  cat_block(d, "stage_f", "Stage at diagnosis"),
  cat_block(d, "eth_f",   "Ethnicity"),
  cat_block(d, "imd_f",   "Deprivation quintile"),
  cat_block(d, "cci_f",   "Charlson comorbidity index"),
  cat_block(d, "route_f", "Route to diagnosis"),
  cat_block(d, "year_f",  "Calendar year")
)

# only label the first row of each characteristic (blank the repeats)
tab1$Characteristic[c(FALSE, tab1$Characteristic[-1] == tab1$Characteristic[-nrow(tab1)])] <- ""

ft <- flextable(tab1)
ft <- set_header_labels(ft, Characteristic = "Characteristic", Level = "",
                        patients = "Patients, n (%)", wait = "Waiting time, days\nmean (SD)")
ft <- bold(ft, part = "header")
ft <- bold(ft, j = "Characteristic", i = which(tab1$Characteristic != ""), part = "body")
ft <- align(ft, j = c("patients", "wait"), align = "right", part = "all")
ft <- valign(ft, valign = "top", part = "body")
ft <- fontsize(ft, size = 9, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
ft <- line_spacing(ft, space = 1.3, part = "all")
ft <- add_footer_lines(ft, paste("Patients as n (%); age at diagnosis as mean (SD).",
                                 "Waiting time is mean (SD) days from diagnosis to decision-to-treat."))
ft <- autofit(ft)

save_as_docx(ft, path = file.path(out_dir, "table1_characteristics.docx"))
write.csv(tab1, file.path(out_dir, "table1_characteristics.csv"), row.names = FALSE)
cat("Table 1 written to Output/table1_characteristics.docx\n")

# representativeness: patient share by cancer alliance -------------------------
alliance <- df %>%
  count(canalliance, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))
print(as.data.frame(alliance), row.names = FALSE)
write.csv(alliance, file.path(out_dir, "table1_alliance.csv"), row.names = FALSE)