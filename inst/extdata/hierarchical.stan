data {
  int Ny; // number of entries
  int Nx; // number of groups
  real Y[Ny]; // continuous response
  int X[Ny]; // group membership indices
  real E_S_low; // empirical SD low 99% HDI boundary
  real E_S_high; // empirical SD high 99% HDI boundary
  real E_M; // empirical mean for each group
  real E_M_low; // empirical mean low 99% HDI boundary
  real E_M_high; // empirical mean high 99% HDI boundary
}

parameters {
  real mu[Nx];
  real<lower=0> sigma[Nx];
  real<lower=1> nu;
  real M;
  real <lower=0>S;
}

model {
  for(i in 1:Ny) {
    Y[i] ~ student_t(nu, mu[X[i]], sigma[X[i]]);
  }
  for(j in 1:Nx) {
    mu[j] ~ normal(M, S);
    // sigma[j] ~ uniform(E_S_low/10.0, E_S_high*10.0);
  }
  nu ~ gamma(2.0, 0.1);
  M ~ normal(E_M, E_S_high*10.0);
  S ~ uniform(E_M_low/10, E_M_high*10);
}
