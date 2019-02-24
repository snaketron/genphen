// Save this file as src/stan_files/
data {
  int Z; // number of all entries
  int  N [Z]; // number tries
  int Y [Z]; // number of hits response
  vector [Z] X; // index of all individuals
}

parameters {
  real alpha;
  real beta;
}


model {
  Y ~ binomial_logit(N, alpha+beta*X);
  alpha ~ student_t(1, 0, 100);
  beta ~ student_t(1, 0, 10);
}
