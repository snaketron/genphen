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
  vector [Ntq+Ntd] alpha;
  vector <lower = 0> [Ntq] sigma;
  vector <lower = 1> [Ntq] nu;
  vector [Ntq+Ntd] mu_beta;
  vector <lower = 0> [Ntq+Ntd] sigma_beta;
  vector <lower = 0> [Ntq+Ntd] tau_beta;
  vector <lower = 1> [Ntq+Ntd] nu_beta;
  vector [Nsk] z [Ntq+Ntd];
}


transformed parameters {
  vector [Nsk] beta [Ntq+Ntd];

  for(t in 1:(Ntq+Ntd)) {
    beta[t] = mu_beta[t] + z[t]*sigma_beta[t]/sqrt(tau_beta[t]);
  }
}


model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ student_t(nu[t], alpha[t] + rows_dot_product(X[i], beta[t]), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + rows_dot_product(X[i], beta[d+Ntq]));
      }
    }
  }
  
  alpha ~ student_t(1, 0, 100);
  for(t in 1:(Ntq+Ntd)) {
    z[t] ~ normal(0, 1);
  }
  nu ~ gamma(2, 0.1);
  sigma ~ cauchy(0, 5);
  mu_beta ~ student_t(1, 0, 10);
  sigma_beta ~ cauchy(0, 5);
  nu_beta ~ gamma(2, 0.1);
  tau_beta ~ gamma(nu_beta/2, nu_beta/2);
}
