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
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_beta;
  real <lower = 1> nu_alpha;
  real <lower = 1> nu_beta;
}


model {
  for(i in 1:Z) {
    Y[i] ~ binomial_logit(N[i], alpha[S[i]] + beta[S[i]]*X[i]);
  }

  mu_alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);

  nu_alpha ~ gamma(2, 0.1);
  nu_beta ~ gamma(2, 0.1);

  sigma_alpha ~ cauchy(0, 1);
  sigma_beta ~ cauchy(0, 1);
}
