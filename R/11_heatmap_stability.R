# Stability heat map: standardised mean days to decision-to-treat by hospital
# and six-month period across the window. Hospitals are ordered by their overall
# standardised mean. Indirect standardisation (a pooled prognostic model giving
# the expected wait in each cell) is used here because six-month cells are small;
# this is a descriptive view of stability rather than the headline estimand.

library(dplyr)
library(lubridate)
library(ggplot2)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# six-month period index within the window
start_date <- min(df$diagmdy)
df <- df %>%
  mutate(bin = floor(as.numeric(interval(start_date, diagmdy), "months") / 6) + 1,
         bin = factor(bin, labels = paste0("H", sort(unique(bin)))))

# pooled prognostic model for expected wait (main case-mix: age + cci)
cv <- code_covariates(df); df <- cv$data
pm <- lm(as.formula(paste("wait ~", paste(c(cv$cont, cv$bin), collapse = " + "))),
         data = df)
df$pred <- predict(pm)
grand <- mean(df$wait)

cell <- df %>%
  group_by(hosp, diag_hosp_canon, bin) %>%
  summarise(n = n(), std = mean(wait) - mean(pred) + grand, .groups = "drop") %>%
  filter(n >= 10)                       # suppress very small cells

order_by <- cell %>% group_by(diag_hosp_canon) %>%
  summarise(overall = mean(std), .groups = "drop") %>% arrange(overall)
cell <- cell %>% mutate(diag_hosp_canon = factor(diag_hosp_canon, levels = order_by$diag_hosp_canon))

p <- ggplot(cell, aes(bin, diag_hosp_canon, fill = std)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = grand, low = "#2c7bb6", mid = "#ffffbf",
                       high = "#d7191c", name = "Std. days") +
  labs(x = "Six-month period", y = "Hospital (ordered by overall mean)",
       title = "Stability of standardised waits over the window") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "heatmap_stability.pdf"), p, width = 7, height = 9)
write.csv(cell, file.path(out_dir, "heatmap_cells.csv"), row.names = FALSE)
cat("heat map saved.\n")