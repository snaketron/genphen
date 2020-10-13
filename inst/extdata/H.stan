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
  vector <lower = 0> [Ntq] sigma;
  vector [Ntq+Ntd] mu_alpha;
  vector <lower = 0> [Ntq+Ntd] sigma_alpha;
  vector <lower = 0> [Ntq+Ntd] tau_alpha;
  vector <lower = 1> [Ntq+Ntd] nu_alpha;
  vector [Ntq+Ntd] mu_beta;
  vector <lower = 0> [Ntq+Ntd] sigma_beta;
  vector <lower = 0> [Ntq+Ntd] tau_beta;
  vector <lower = 1> [Ntq+Ntd] nu_beta;
  vector [Nsk] za [Ntq+Ntd];
  vector [Nsk] zb [Ntq+Ntd];
}


transformed parameters {
  vector [Nsk] alpha [Ntq+Ntd];
  vector [Nsk] beta [Ntq+Ntd];

  for(t in 1:(Ntq+Ntd)) {
    alpha[t] = mu_alpha[t] + za[t]*sigma_alpha[t]/sqrt(tau_alpha[t]);
    beta[t] = mu_beta[t] + zb[t]*sigma_beta[t]/sqrt(tau_beta[t]);
  }
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
  for(t in 1:(Ntq+Ntd)) {
    za[t] ~ std_normal();
    zb[t] ~ std_normal();
  }
  
  sigma ~ cauchy(0, 5);
  mu_beta ~ student_t(1, 0, 10);
  sigma_beta ~ cauchy(0, 5);
  nu_beta ~ gamma(2, 0.1);
  tau_beta ~ gamma(nu_beta/2, nu_beta/2);
  
  mu_alpha ~ student_t(1, 0, 100);
  sigma_alpha ~ cauchy(0, 5);
  nu_alpha ~ gamma(2, 0.1);
  tau_alpha ~ gamma(nu_alpha/2, nu_alpha/2);
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
