data {
  int Z; // number of all entries
  real Y[Z]; // number of hits response
  int X[Z]; // index of all individuals
}

parameters {
  real alpha;
  real beta;
  real <lower = 0> sigma;
  real <lower = 1> nu;
}

model {
  for(i in 1:Z) {
    Y[i] ~ student_t(nu, alpha + beta*X[i], sigma);
    
  }
  alpha ~ normal(0, 100);
  beta ~ normal(0, 10);
  nu ~ gamma(2.0, 0.1);
  sigma ~ cauchy(0, 1);
}
