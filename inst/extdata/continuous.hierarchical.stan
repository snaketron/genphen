data {
  int Z; // number of all entries
  int S_N; // number of all SNPs
  int S[Z]; // index of all SNPs
  real Y[Z]; // number of hits response
  int X[Z]; // index of all individuals
}

parameters {
  real alpha[S_N];
  real beta[S_N];
  real <lower = 1> nu;
  real <lower = 0> sigma;
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_beta;
  real <lower = 1> nu_alpha;
  real <lower = 1> nu_beta;
}

model {
  for(i in 1:Z) {
    Y[i] ~ student_t(nu, alpha[S[i]] + beta[S[i]]*X[i], sigma);
  }
  
  for(j in 1:S_N) {
    alpha[j] ~ student_t(nu_alpha, mu_alpha, sigma_alpha);
    beta[j] ~ student_t(nu_beta, mu_beta, sigma_beta);
  }
  
  mu_alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);
  
  nu ~ gamma(2, 0.1);
  nu_alpha ~ gamma(2, 0.1);
  nu_beta ~ gamma(2, 0.1);

  sigma ~ cauchy(0, 1);
  sigma_alpha ~ cauchy(0, 1);
  sigma_beta ~ cauchy(0, 1);
}
