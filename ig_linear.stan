functions {
  real IG_log (real x, real mu, real lambda){
    // vector [num_elements (x)] prob;
    real lprob;
    lprob = log((lambda/(2*pi()*(x^3)))^0.5 * exp(-lambda*(x - mu)^2/(2*mu^2*x)));
    return lprob;
  }
}


data {
  int<lower=0> I;
  int<lower=0> J;
  matrix[I,J] x;
  matrix[I,J] t;
}

parameters {
  real mu;
  real<lower=0> w;
}


model {
  w ~ gamma(1,1); //scale
  mu ~ normal(0, 1); //shape
  for (i in 1:I){
    for (j in 1:J) {
      x[i,j] ~ IG_log (mu * t[i,j], w * mu^2 * t[i,j]^2);
    }
  }
}
