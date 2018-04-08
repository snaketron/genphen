data {
  int Ny; // number of entries
  int Nx; // number of groups
  int Ns; // number of sites
  int S[Ny]; // number of sites
  real Y[Ny]; // continuous response
  int X[Ny]; // group membership indices
  real E_sigma; // empirical SD for each group
  real E_mu; // empirical mean for each group
}

parameters {
  real mu[Nx];
  real<lower=0> sigma[Nx];
  real<lower=1> nu[Ns];
}

model {
  for(i in 1:Ny) {
    Y[i] ~ student_t(nu[S[i]], mu[X[i]], sigma[X[i]]);
  }
  for(j in 1:Nx) {
    mu[j] ~ normal(E_mu, E_sigma*100.0);
    sigma[j] ~ uniform(E_sigma/100.0, E_sigma*100.0);
  }
  for(k in 1:Ns) {
    nu[k] ~ gamma(2.0, 0.1);
  }
}
