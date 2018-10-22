data {
  int Z; // number of all entries
  int N[Z]; // number tries
  int Y[Z]; // number of hits response
  int X[Z]; // index of all individuals
}

parameters {
  real alpha;
  real beta;
}


model {
  for(i in 1:Z) {
    Y[i] ~ binomial_logit(N[i], alpha+beta*X[i]);
  }
  alpha ~ normal(0, 100);
  beta ~ normal(0, 10);
}
