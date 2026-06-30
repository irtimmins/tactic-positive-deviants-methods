# 00  master - run the colon positive-deviance analysis end to end
# -----------------------------------------------------------------------------
# Run from the project root (the folder containing R/ and Output/). Each step
# writes to Output/ and the next reads it, so order matters. Settings live in
# R/01_config.R; edit base_dir there first.
#
# Step 03 builds the provider crosswalks from the ODS API and the curated site
# Excel. It is a run-once, then-curate step: run it, review the two printed
# tables and the crosswalk CSVs (and add manual overrides), and only then run 04
# onward. It is commented out of the automatic run below for that reason.
#
# Prerequisites: balancer, rstan, lme4, dplyr, tidyr, ggplot2, lubridate; for 03
# also readxl, httr2, tidyverse; and the Stan model file R/dp_normal_cont.stan.

r_dir <- "R"
step <- function(file) {
  message("\n========== ", file, " ==========")
  source(file.path(r_dir, file), local = new.env())
}

step("02_build_analysis_dataset.R")     # clinical funnel -> eligible set (raw codes)
step("03_provider_level_qc.R") # run once, then curate (needs API + Excel)
step("04_finalise_analysis_dataset.R")  # apply crosswalks + volume -> final dataset
step("05_table1_characteristics.R")     # Table 1 and alliance representativeness
step("06_estimation_weights.R")         # case-mix standardisation (balancing weights)
step("07_shrinkage.R")                  # Bayesian shrinkage of hospital estimates
step("08_ranks_caterpillars.R")         # posterior ranks, caterpillars, candidates
step("09_method_comparison.R")          # ranking-method comparison + funnel plots
step("10_strata_concordance.R")         # consistency across comorbidity strata
step("11_heatmap_stability.R")          # stability over six-month periods

message("\nAnalysis complete. Outputs in Output/.")
