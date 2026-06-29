#test
df <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/colon_cohort_cci_2015_2022.rds")
names(df)

summary(df)
library(dplyr)

cat("=== 1. variables, type, completeness ===\n")
tibble(
  variable = names(df),
  type     = vapply(df, function(x) class(x)[1], character(1)),
  pct_complete = round(100 * vapply(df, function(x) mean(!is.na(x)), numeric(1)), 1)
) %>% as.data.frame() %>% print()

cat("\n=== 2. levels for categorical columns (<=25 distinct) ===\n")
for (v in names(df)) {
  x <- df[[v]]
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    u <- unique(x[!is.na(x)])
    if (length(u) <= 25)
      cat(sprintf("  %-26s : %s\n", v, paste(sort(as.character(u)), collapse = " | ")))
    else
      cat(sprintf("  %-26s : <%d distinct - id/free-text>\n", v, length(u)))
  }
}
