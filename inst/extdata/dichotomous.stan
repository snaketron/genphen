data {
  int Ny; // number of entries
  int Nx; // number of groups
  int Y[Ny]; // continuous response
  int X[Ny]; // group membership indices
}

parameters {
  real <lower = 0, upper = 1> mu[Nx];
}

model {
  for(i in 1:Ny) {
    Y[i] ~ bernoulli(mu[X[i]]);
  }
  for(j in 1:Nx) {
    mu[j] ~ beta(1.0/2.0, 1.0/2.0);
  }
}
