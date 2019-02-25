// Save this file as src/stan_files/
data {
  int Z; // number of all entries
  vector [Z] Y; // number of hits response
  vector [Z] X; // index of all individuals
}

parameters {
  real alpha;
  real beta;
  real <lower = 0> sigma;
  real <lower = 1> nu;
}

model {
  Y ~ student_t(nu, alpha + beta*X, sigma);
  alpha ~ normal(0, 100);
  beta ~ normal(0, 10);
  nu ~ gamma(2.0, 0.1);
  sigma ~ cauchy(0, 5);
}
