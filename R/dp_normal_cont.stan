// Normal-normal shrinkage model for hospital-level summaries, continuous
// outcome version of dp_normal.stan. The template file is written for a
// proportion in [0,1] (beta prior, cauchy(0,0.5)); a waiting time in days is
// unbounded and on a much larger scale, so the [0,1] bounds are dropped and the
// prior scales are passed in from R. The structure (a single normal random
// effect with a half-cauchy scale) is unchanged.

data {
  int<lower=0> J;                    // number of hospitals
  real y_site_obs[J];                // observed standardised mean per hospital
  real<lower=0> sigma_site_obs[J];   // its standard error
  real prior_mu_mean;                // prior mean for the grand mean
  real<lower=0> prior_mu_sd;         // prior sd for the grand mean
  real<lower=0> prior_tau_scale;     // half-cauchy scale for between-hospital sd
}

parameters {
  real mu_true;                      // grand mean
  real<lower=0> sigma_true;          // between-hospital sd
  real y_site_true[J];               // latent hospital means
}

model {
  mu_true    ~ normal(prior_mu_mean, prior_mu_sd);
  sigma_true ~ cauchy(0, prior_tau_scale);

  y_site_true ~ normal(mu_true, sigma_true);
  y_site_obs  ~ normal(y_site_true, sigma_site_obs);
}
