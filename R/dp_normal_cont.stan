// Normal-normal shrinkage model for hospital-level summaries, continuous
// outcome version of dp_normal.stan. The template file is written for a
// proportion in [0,1] (beta prior, cauchy(0,0.5)); a waiting time in days is
// unbounded and on a much larger scale, so the [0,1] bounds are dropped and the
// prior scales are passed in from R.
//
// Non-centred parameterisation: hospital effects are modelled as standardised
// deviations z_site ~ normal(0,1) and rescaled to y_site_true. This is the same
// model as the centred form but avoids the funnel geometry that causes divergent
// transitions when the between-hospital sd is small relative to the per-hospital
// standard errors (e.g. low-signal strata such as high comorbidity).

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
  vector[J] z_site;                  // standardised hospital deviations
}

transformed parameters {
  real y_site_true[J];               // latent hospital means
  for (j in 1:J)
    y_site_true[j] = mu_true + sigma_true * z_site[j];
}

model {
  mu_true    ~ normal(prior_mu_mean, prior_mu_sd);
  sigma_true ~ cauchy(0, prior_tau_scale);

  z_site     ~ normal(0, 1);
  y_site_obs ~ normal(y_site_true, sigma_site_obs);
}
