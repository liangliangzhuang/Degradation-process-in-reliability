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
  w ~ gamma(1,1);
  mu ~ normal(0, 100/w);
  for (i in 1:I){
    for (j in 1:J) {
      x[i,j] ~ normal(mu * t[i,j], 1/w * t[i,j]);
    }
  }
}
