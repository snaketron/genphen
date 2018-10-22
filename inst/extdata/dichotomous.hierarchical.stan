data {
  int Z; // number of all entries
  int S_N; // number of all SNPs
  int S[Z]; // index of all SNPs
  int N[Z]; // number tries
  int Y[Z]; // number of hits response
  int X[Z]; // index of all individuals
}

parameters {
  real alpha[S_N];
  real beta[S_N];
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sd_alpha;
  real <lower = 0> sd_beta;
  real <lower = 1> nu_alpha;
  real <lower = 1> nu_beta;
  // real <lower = 1> sigma;
  // real p;
}


model {
  for(i in 1:Z) {
    Y[i] ~ binomial_logit(N[i], alpha[S[i]] + beta[S[i]]*X[i]);
  }
  // vector[Z] mu;
  // for(i in 1:Z) {
  //   mu[i] = alpha[S[i]] + beta[S[i]]*X[i];
  // }
  // Y ~ binomial_logit(N, p);
  // p ~ normal(mu, sigma);
  // for(j in 1:S_N) {
  //   alpha[j] ~ student_t(nu_alpha, mu_alpha, sd_alpha);
  //   beta[j] ~ student_t(nu_beta, mu_beta, sd_beta);
  // }

  mu_alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);

  nu_alpha ~ gamma(2, 0.1);
  nu_beta ~ gamma(2, 0.1);

  sd_alpha ~ cauchy(0, 1);
  sd_beta ~ cauchy(0, 1);
  // sigma ~ cauchy(0, 1);
}
