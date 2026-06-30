# 00  master - run the colon positive-deviance analysis end to end
# -----------------------------------------------------------------------------
# Run from the project root (the folder containing R/ and Output/). Each step
# writes its outputs to Output/ and the next step reads them, so the order
# matters. Settings live in R/01_config.R; edit base_dir there first.
#
# Prerequisites: R packages balancer, rstan, lme4, dplyr, tidyr, ggplot2,
# lubridate, and the Stan model file R/dp_normal_cont.stan. The shrinkage and
# estimation steps are the slow ones.

r_dir <- "R"
step <- function(file) {
  message("\n========== ", file, " ==========")
  source(file.path(r_dir, file), local = new.env())
}

step("02_build_analysis_dataset.R")   # full cohort -> analysis sample + flowchart
step("03_table1_characteristics.R")   # Table 1 and alliance representativeness
step("04_estimation_weights.R")       # case-mix standardisation (balancing weights)
step("05_shrinkage.R")                # Bayesian shrinkage of hospital estimates
step("06_ranks_caterpillars.R")       # posterior ranks, caterpillars, candidates
step("07_method_comparison.R")        # ranking-method comparison + funnel plots
step("08_strata_concordance.R")       # consistency across comorbidity strata
step("09_heatmap_stability.R")        # stability over six-month periods

message("\nAnalysis complete. Outputs in Output/.")
