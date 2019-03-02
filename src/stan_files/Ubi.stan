data {
  int N; // number of all entries
  int Ntc; // number of continuous traits
  int Ntd; // number of dichotomous traits
  int Ns; // number of all SNPs
  real Yc[N, Ntc] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Ns] X [N]; // index of all individuals
}

parameters {
  vector [Ns] alpha [Ntc+Ntd];
  vector [Ns] beta [Ntc+Ntd];
  vector <lower = 1> [Ntc] nu;
  vector <lower = 0> [Ntc] sigma;
}



model {
  if(Ntc > 0) {
    for(t in 1:Ntc) {
      for(i in 1:N) {
        Yc[i,t] ~ student_t(nu[t], alpha[t] + rows_dot_product(X[i], to_vector(beta[t])), sigma[t]);
      }
    }
  }
  
  if(Ntd > 0) {
    for(d in 1:Ntd) {
      for(i in 1:N) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntc] + rows_dot_product(X[i], to_vector(beta[d+Ntc])));
      }
    }
  }
  
  for(t in 1:(Ntc+Ntd)) {
    alpha[t] ~ student_t(1, 0, 100);
    beta[t] ~ student_t(1, 0, 10);
  }
  nu ~ gamma(2, 0.1);
  sigma ~ cauchy(0, 5);
}
