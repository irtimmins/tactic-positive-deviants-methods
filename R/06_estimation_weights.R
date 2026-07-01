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
library(ggplot2)

source("R/01_config.R")
df <- readRDS(file.path(out_dir, "analysis_data.rds"))

# reweight one data frame to its own population means and summarise by hospital.
# cont/bin are the balance covariates; the prognostic model uses the same set.
run_standardise <- function(d, cont, bin, lambda = lambda_main, uplim = NULL,
                            ref = NULL, target_data = NULL) {
  d <- d %>% arrange(hosp)
  if (is.null(ref)) ref <- ref_moments(d, cont, bin)
  X <- make_std_matrix_ref(d, cont, bin, ref)
  
  # drop covariates that are constant within d: they carry no balancing
  # information and would divide by zero (e.g. a within-period indicator after
  # the data have been split by period).
  keep <- apply(X, 2, function(col) { s <- sd(col); is.finite(s) && s > 0 })
  X <- X[, keep, drop = FALSE]
  
  Z <- d$hosp
  args <- list(X = X, target = rep(0, ncol(X)), Z = Z,
               lambda = lambda, exact_global = FALSE)
  if (!is.null(uplim)) args$uplim <- uplim
  std_out <- do.call(standardize, args)
  d$w <- extract_weights(std_out)
  d$y <- d$wait                     # site_summary works on a column named y
  
  # pooled prognostic model on the raw covariates for residual balancing
  form <- as.formula(paste("wait ~", paste(c(cont, bin), collapse = " + ")))
  pm <- lm(form, data = d)
  d$resid <- resid(pm)
  # model prediction averaged over the target population (defaults to d itself);
  # for the later improvement period this is the baseline population so the
  # augmented estimate also targets the baseline case-mix.
  d$canonical <- if (is.null(target_data)) mean(fitted(pm))
  else mean(predict(pm, newdata = target_data))
  
  site <- site_summary(d) %>%
    left_join(distinct(d, hosp, diag_hosp = diag_hosp_canon), by = "hosp")
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

# main analysis: case-mix (age + cci) standardisation ------------------------
# NB the main model deliberately excludes season and calendar year (they are a
# sensitivity analysis, script 12) and all other patient factors.
cv    <- code_covariates(df)
trade <- balance_tradeoff(cv$data, cv$cont, cv$bin)
cat("balance vs effective n by lambda:\n"); print(round(trade, 3))
write.csv(trade, file.path(out_dir, "lambda_tradeoff.csv"), row.names = FALSE)

# bias-variance trade-off curve (after Keele et al). Each point is one lambda:
# how much case-mix imbalance the weights remove (percent bias reduced) against
# the average effective sample size they keep. Small lambda sits top-left (most
# bias removed, fewest effective patients); larger lambda moves down and right.
# The elbow is the practical sweet spot; the working value lambda_main is marked.
trade_curve <- trade %>%
  mutate(pct_bias_removed = 100 * bias_removed,
         is_main = abs(lambda - lambda_main) < 1e-9)

p_trade <- ggplot(trade_curve, aes(mean_eff_n, pct_bias_removed)) +
  geom_path(colour = "grey60") +
  geom_point(aes(colour = is_main), size = 2) +
  geom_text(aes(label = lambda), vjust = -0.8, size = 2.8) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick"),
                      guide = "none") +
  labs(x = "Average effective sample size per hospital",
       y = "Case-mix bias removed (%)",
       title = "Bias reduction against effective sample size across lambda",
       subtitle = "each point is a lambda value; the working value is highlighted") +
  theme_bw()

ggsave(file.path(out_dir, "lambda_tradeoff.pdf"), p_trade, width = 7, height = 5)

fit_main <- run_standardise(cv$data, cv$cont, cv$bin, lambda_main)
site_sustained <- fit_main$site
saveRDS(fit_main,       file.path(out_dir, "fit_primary.rds"))
saveRDS(site_sustained, file.path(out_dir, "site_sustained.rds"))

# improvement estimand -------------------------------------------------------
# baseline period = first half of the window; later period = second half. Both
# are standardised to the SAME reference: the case-mix of the baseline period.
# The period split is the design here, so calendar terms are not covariates.
half_n <- df %>% count(hosp, period) %>%
  tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)
ok_hosp <- half_n %>% filter(first >= min_per_year, second >= min_per_year) %>% pull(hosp)

base_dat <- filter(df, period == "first",  hosp %in% ok_hosp)
late_dat <- filter(df, period == "second", hosp %in% ok_hosp)
cvb <- code_covariates(base_dat); cvl <- code_covariates(late_dat)
ref_base <- ref_moments(cvb$data, cvb$cont, cvb$bin)    # the fixed target

site_p1 <- run_standardise(cvb$data, cvb$cont, cvb$bin, lambda_main,
                           ref = ref_base, target_data = cvb$data)$site
site_p2 <- run_standardise(cvl$data, cvl$cont, cvl$bin, lambda_main,
                           ref = ref_base, target_data = cvb$data)$site

site_improve <- site_p1 %>%
  select(hosp, diag_hosp, stand1 = stand_adj, se1 = se_adj_pool, n1 = n) %>%
  inner_join(site_p2 %>% select(hosp, stand2 = stand_adj, se2 = se_adj_pool, n2 = n),
             by = "hosp") %>%
  mutate(delta = stand2 - stand1,            # negative = faster over time
         se_delta = sqrt(se1^2 + se2^2))
saveRDS(site_improve, file.path(out_dir, "site_improve.rds"))

# comorbidity strata ---------------------------------------------------------
# within a stratum comorbidity is near-constant, so adjust for age only and
# standardise to the stratum population.
for (st in levels(df$cci_strata)) {
  cvs <- code_covariates(filter(df, cci_strata == st), cci = "none")
  fit_st <- run_standardise(cvs$data, cvs$cont, cvs$bin, lambda_main)
  saveRDS(fit_st$site,
          file.path(out_dir, sprintf("site_strata_%s.rds", gsub("[^0-9a-zA-Z]", "", st))))
}