// Save this file as src/stan_files/
data {
  int Z; // number of all entries
  int S_N; // number of all SNPs
  int S[Z]; // index of all SNPs
  int N[Z]; // number tries
  int Y[Z]; // number of hits response
  int X[Z]; // index of all individuals
}

parameters {
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_beta;
  real <lower = 1> nu_alpha;
  real <lower = 1> nu_beta;
  real <lower = 0> tau_alpha;
  real <lower = 0> tau_beta;
  real z_beta[S_N];
  real z_alpha[S_N];
}


transformed parameters {
  real alpha[S_N];
  real beta[S_N];
  
  for(s in 1:S_N) {
    alpha[s] = mu_alpha + z_alpha[s]*sigma_alpha/sqrt(tau_alpha);
    beta[s] = mu_beta + z_beta[s]*sigma_beta/sqrt(tau_beta);
  }
}

model {
  for(i in 1:Z) {
    Y[i] ~ binomial_logit(N[i], alpha[S[i]] + beta[S[i]]*X[i]);
  }

  mu_alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);
  
  nu_alpha ~ gamma(2, 0.1);
  nu_beta ~ gamma(2, 0.1);

  sigma_alpha ~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
  
  z_alpha ~ normal(0, 1);
  z_beta ~ normal(0, 1);
  
  tau_alpha ~ gamma(nu_alpha/2, nu_alpha/2);
  tau_beta ~ gamma(nu_beta/2, nu_beta/2);
}
