data {
  int Ny; // number of entries
  int Nx; // number of groups
  real Y[Ny]; // continuous response
  int X[Ny]; // group membership indices
  real E_sigma[Nx]; // empirical SD for each group
  real E_mu[Nx]; // empirical mean for each group
}

parameters {
  real mu[Nx];
  real<lower=0> sigma[Nx];
  real<lower=1> nu[Nx];
}

model {
  for(i in 1:Ny) {
    Y[i] ~ student_t(nu[X[i]], mu[X[i]], sigma[X[i]]);
  }
  for(j in 1:Nx) {
    mu[j] ~ normal(E_mu[j], E_sigma[j]);
    sigma[j] ~ uniform(E_sigma[j]/100, E_sigma[j]*100);
    nu[j] ~ gamma(2, 0.1);
  }
}
