// Normal-normal shrinkage with a switchable prior on the between-hospital sd,
// for the prior-sensitivity analysis (script 16). Same non-centred model as
// dp_normal_cont.stan; only the sigma_true prior changes, selected by tau_prior:
//   1  half-normal(0, tau_scale)
//   2  uniform(0, tau_upper)          (flat)
//   3  half-student_t(tau_df, 0, tau_scale)
//   4  half-cauchy(0, tau_scale)      (the main model's prior)
// The uniform needs a finite upper bound; sigma_true is bounded above by
// tau_upper, set large for the other priors so it never binds.

data {
  int<lower=0> J;                    // number of hospitals
  real y_site_obs[J];                // observed standardised mean per hospital
  real<lower=0> sigma_site_obs[J];   // its standard error
  real prior_mu_mean;                // prior mean for the grand mean
  real<lower=0> prior_mu_sd;         // prior sd for the grand mean
  int<lower=1, upper=4> tau_prior;   // which prior on the between-hospital sd
  real<lower=0> tau_scale;           // scale for the normal / student_t / cauchy priors
  real<lower=0> tau_upper;           // upper bound (defines the uniform prior)
  real<lower=0> tau_df;              // degrees of freedom for the student_t prior
}

parameters {
  real mu_true;                          // grand mean
  real<lower=0, upper=tau_upper> sigma_true;  // between-hospital sd
  vector[J] z_site;                      // standardised hospital deviations
}

transformed parameters {
  real y_site_true[J];                   // latent hospital means
  for (j in 1:J)
    y_site_true[j] = mu_true + sigma_true * z_site[j];
}

model {
  mu_true ~ normal(prior_mu_mean, prior_mu_sd);

  if (tau_prior == 1)
    sigma_true ~ normal(0, tau_scale);                 // half-normal (lower = 0)
  else if (tau_prior == 2)
    sigma_true ~ uniform(0, tau_upper);                // flat
  else if (tau_prior == 3)
    sigma_true ~ student_t(tau_df, 0, tau_scale);      // half-student_t
  else
    sigma_true ~ cauchy(0, tau_scale);                 // half-cauchy

  z_site     ~ normal(0, 1);
  y_site_obs ~ normal(y_site_true, sigma_site_obs);
}
