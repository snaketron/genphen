data {
  int N; // number of all entries
  int Ntq; // number of continuous traits
  int Ntd; // number of dichotomous traits
  int Ns; // number of all SNPs
  int Nsk; // number of SNPs x substitutions
  real Yq[N, Ntq] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Nsk] X [N]; // index mapping function
}


parameters {
  vector [Nsk] alpha [Ntq+Ntd];
  vector [Nsk] beta [Ntq+Ntd];
  vector <lower = 0> [Ntq] sigma;
}


model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + X[i] .* beta[t], sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + X[i] .* beta[d+Ntq]);
      }
    }
  }
  
  sigma ~ cauchy(0, 5);
  for(t in 1:(Ntq+Ntd)) {
    beta[t] ~ student_t(1, 0, 10);
    alpha[t] ~ student_t(1, 0, 100);
  }
}

generated quantities {
  vector [Nsk] Y_hat [Ntq+Ntd, 2]; 
  
  
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Y_hat[t, 1] = to_vector(normal_rng(alpha[t] + beta[t], sigma[t]));
        Y_hat[t, 2] = to_vector(normal_rng(alpha[t] + -beta[t], sigma[t]));
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Y_hat[d+Ntq, 1] = to_vector(bernoulli_rng(inv_logit(alpha[d+Ntq] + beta[d+Ntq])));
        Y_hat[d+Ntq, 2] = to_vector(bernoulli_rng(inv_logit(alpha[d+Ntq] + -beta[d+Ntq])));
      }
    }
  }
}
