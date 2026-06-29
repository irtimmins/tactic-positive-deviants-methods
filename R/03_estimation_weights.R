# Directly standardised hospital waiting times via balancing weights.
# Mirrors the Estimation_Step / Post_Estimation steps of the template:
#  - build the standardised balance matrix and reweight each hospital to the
#    population covariate means using balancer::standardize()
#  - a pooled prognostic model gives residuals for the augmented (residual
#    balancing) estimate and the population-mean prediction
#  - per-hospital weighted and augmented means, effective n and pooled SEs
# Produces inputs for the shrinkage step for: the sustained estimand (whole
# window), the improvement estimand (second half minus first), and the two
# comorbidity strata.

library(balancer)
library(dplyr)

source("00_config_funcs.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# reweight one data frame to its own population means and summarise by hospital.
# cont/bin are the balance covariates; the prognostic model uses the same set.
run_standardise <- function(d, cont, bin, lambda = lambda_main, uplim = NULL) {
  d <- d %>% arrange(hosp)
  X    <- make_std_matrix(d, cont, bin)
  Z    <- d$hosp
  tgt0 <- rep(0, ncol(X))
  args <- list(X = X, target = tgt0, Z = Z, lambda = lambda, exact_global = FALSE)
  if (!is.null(uplim)) args$uplim <- uplim
  std_out <- do.call(standardize, args)
  d$w <- extract_weights(std_out)
  d$y <- d$wait                     # site_summary works on a column named y
  
  # pooled prognostic model on the raw covariates for residual balancing
  form <- as.formula(paste("wait ~", paste(c(cont, bin), collapse = " + ")))
  pm <- lm(form, data = d)
  d$resid     <- resid(pm)
  d$canonical <- mean(fitted(pm))   # constant for a pooled model (= grand mean)
  
  site <- site_summary(d) %>%
    left_join(distinct(d, hosp, diag_hosp), by = "hosp")
  list(d = d, site = site, std_out = std_out, model = pm)
}

# balance vs effective-sample-size trade-off across lambda --------------------
# prognostic weights for the imbalance metric come from the pooled model on the
# standardised covariates, so imbalance is summarised on the outcome scale.
balance_tradeoff <- function(d, cont, bin, grid = lambda_grid) {
  d  <- d %>% arrange(hosp)
  X  <- make_std_matrix(d, cont, bin)
  Z  <- d$hosp
  pm <- lm(as.formula(paste("y_std ~", paste(colnames(X), collapse = " + "))),
           data = data.frame(y_std = d$wait, X))
  beta <- coef(pm)[-1]
  beta[is.na(beta)] <- 0
  
  hosp_means <- rowsum(X, Z) / as.numeric(table(Z))    # hospital x covariate
  unw <- as.numeric(abs((hosp_means %*% beta)))
  
  res <- t(sapply(grid, function(l) {
    so <- standardize(X, rep(0, ncol(X)), Z, lambda = l, exact_global = FALSE)
    w  <- extract_weights(so)
    wm <- (t(so$weights) %*% X)                        # hospital x covariate
    wt <- as.numeric(abs(wm %*% beta))
    ne <- tapply(w, Z, function(x) sum(x)^2 / sum(x^2))
    c(lambda = l,
      bias_removed = 1 - mean(wt) / mean(unw),
      mean_eff_n   = mean(ne),
      mean_deff    = mean(ne / as.numeric(table(Z))))
  }))
  as.data.frame(res)
}

# sustained estimand (whole window) ------------------------------------------
trade <- balance_tradeoff(df, cont_vars, bin_primary)
cat("balance vs effective n by lambda:\n"); print(round(trade, 3))
write.csv(trade, file.path(out_dir, "lambda_tradeoff.csv"), row.names = FALSE)

fit_primary <- run_standardise(df, cont_vars, bin_primary, lambda_main)
fit_full    <- run_standardise(df, cont_vars, bin_full,    lambda_main)

site_sustained <- fit_primary$site
saveRDS(fit_primary, file.path(out_dir, "fit_primary.rds"))
saveRDS(fit_full,    file.path(out_dir, "fit_full.rds"))
saveRDS(site_sustained, file.path(out_dir, "site_sustained.rds"))

# improvement estimand -------------------------------------------------------
# standardise each half of the window to the same (whole-window) target by
# reusing the population means; the difference is the change in standardised
# performance. hospitals need enough volume in both halves.
half_n <- df %>% count(hosp, period) %>%
  tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)
ok_hosp <- half_n %>% filter(first > min_volume / 2, second > min_volume / 2) %>% pull(hosp)

site_p1 <- run_standardise(filter(df, period == "first",  hosp %in% ok_hosp),
                           cont_vars, bin_primary, lambda_main)$site
site_p2 <- run_standardise(filter(df, period == "second", hosp %in% ok_hosp),
                           cont_vars, bin_primary, lambda_main)$site

site_improve <- site_p1 %>%
  select(hosp, diag_hosp, stand1 = stand_adj, se1 = se_adj_pool, n1 = n) %>%
  inner_join(site_p2 %>% select(hosp, stand2 = stand_adj, se2 = se_adj_pool, n2 = n),
             by = "hosp") %>%
  mutate(delta = stand2 - stand1,            # negative = faster over time
         se_delta = sqrt(se1^2 + se2^2))
saveRDS(site_improve, file.path(out_dir, "site_improve.rds"))

# comorbidity strata ---------------------------------------------------------
# within a stratum comorbidity is near-constant, so balance on age and calendar
# terms only and standardise to the stratum population.
strata_cont <- "agediag"
for (st in levels(df$cci_strata)) {
  fit_st <- run_standardise(filter(df, cci_strata == st),
                            strata_cont, bin_primary, lambda_main)
  saveRDS(fit_st$site,
          file.path(out_dir, sprintf("site_strata_%s.rds", gsub("[^0-9a-zA-Z]", "", st))))
}

cat("\nstandardisation done. sustained hospitals:", nrow(site_sustained),
    "| improvement hospitals:", nrow(site_improve), "\n")