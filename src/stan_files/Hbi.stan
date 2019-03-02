data {
  int N; // number of all entries
  int Ntc; // number of continuous traits
  int Ntd; // number of dichotomous traits
  int Ns; // number of all SNPs
  int Nsk; // number of SNPs x substitutions
  real Yc[N, Ntc] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  int X [N, Ns]; // index mapping function 
}


parameters {
  vector [Ntc+Ntd] alpha;
  vector <lower = 0> [Ntc] sigma;
  vector <lower = 1> [Ntc] nu;
  vector [Ntc+Ntd] mu_beta;
  vector <lower = 0> [Ntc+Ntd] sigma_beta;
  vector <lower = 0> [Ntc+Ntd] tau_beta;
  vector <lower = 1> [Ntc+Ntd] nu_beta;
  vector [Ns] z_beta [Ntc+Ntd];
}


transformed parameters {
  vector [Ns] beta [Ntc+Ntd];
  
  for(t in 1:(Ntc+Ntd)) {
    for(s in 1:Ns) {
      beta[t][s] = mu_beta[t] + z_beta[t][s]*sigma_beta[t]/sqrt(tau_beta[t]);
    }
  }
}


model {
  if(Ntc > 0) {
    for(t in 1:Ntc) {
      for(i in 1:N) {
        Yc[i,t] ~ student_t(nu[t], alpha[t] + rows_dot_product(X[i], beta[t]), sigma[t]);
      }
    }
  }
  
  if(Ntd > 0) {
    for(d in 1:Ntd) {
      for(i in 1:N) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntc] + rows_dot_product(X[i], beta[t]));
      }
    }
  }
  
  for(t in 1:(Ntc+Ntd)) {
    alpha[t] ~ student_t(1, 0, 100);
    z_beta[t] ~ normal(0, 1);
  }
  nu ~ gamma(2, 0.1);
  sigma ~ cauchy(0, 5);
  
  mu_beta ~ student_t(1, 0, 10);
  sigma_beta ~ cauchy(0, 5);
  nu_beta ~ gamma(2, 0.1);
  tau_beta ~ gamma(nu_beta/2, nu_beta/2);
}
