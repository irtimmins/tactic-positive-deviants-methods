# Table 1: patient and system characteristics, with median (IQR) waiting time
# in each group. Also a short representativeness check by cancer alliance.

library(dplyr)
library(tidyr)

source("R/00_config_funcs.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# summarise the outcome within levels of a categorical variable
by_cat <- function(d, var, label) {
  d %>%
    mutate(.lev = factor(.data[[var]])) %>%
    group_by(.lev) %>%
    summarise(
      n = n(),
      pct = 100 * n() / nrow(d),
      wait_median = median(wait),
      wait_q1 = quantile(wait, .25),
      wait_q3 = quantile(wait, .75),
      .groups = "drop"
    ) %>%
    transmute(variable = label, level = as.character(.lev),
              n, pct, wait_median, wait_q1, wait_q3)
}

# a continuous row (age) reported as mean (SD) rather than counts
age_row <- df %>% summarise(
  variable = "Age at diagnosis", level = sprintf("mean %.1f (SD %.1f)", mean(age), sd(age)),
  n = n(), pct = 100, wait_median = median(wait),
  wait_q1 = quantile(wait, .25), wait_q3 = quantile(wait, .75)
)

age_grp <- df %>%
  mutate(age_grp = cut(age, c(-Inf, 50, 60, 70, 80, Inf),
                       labels = c("<50","50-59","60-69","70-79","80+"))) %>%
  by_cat("age_grp", "Age group")

tab1 <- bind_rows(
  tibble(variable = "All patients", level = "",
         n = nrow(df), pct = 100, wait_median = median(df$wait),
         wait_q1 = quantile(df$wait, .25), wait_q3 = quantile(df$wait, .75)),
  age_row,
  age_grp,
  by_cat(df, "sexf",       "Sex"),
  by_cat(df, "stage",      "Stage"),
  by_cat(df, "cci_group",  "Charlson group"),
  by_cat(df, "cci_strata", "Comorbidity stratum"),
  by_cat(df, "eth",        "Ethnicity"),
  by_cat(df, "imd",        "Deprivation quintile"),
  by_cat(df, "route",      "Route to diagnosis"),
  by_cat(df, "period",     "Window half")
) %>%
  mutate(across(c(pct, wait_median, wait_q1, wait_q3), ~ round(., 1)))

print(as.data.frame(tab1), row.names = FALSE)
write.csv(tab1, file.path(out_dir, "table1_characteristics.csv"), row.names = FALSE)

# representativeness: patient share by cancer alliance
alliance <- df %>%
  count(canalliance, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(desc(n))
print(as.data.frame(alliance), row.names = FALSE)
write.csv(alliance, file.path(out_dir, "table1_alliance.csv"), row.names = FALSE)
