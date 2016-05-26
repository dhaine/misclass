// 3-level hierarchical logistic regression model
data {
  // Define variables in data
  int<lower=1> N;  // nbr of observations
  int<lower=0,upper=1> y[N];  // outcome
  int<lower=1> K;  // number of fixed effects 
  matrix[N, K] X;  // FE design matrix 
  vector[K] X_means;  // column means of X
  // Random effects for herd
  int<lower=1> J_1[N];  // herd RE levels
  int<lower=1> N_1;  // nbr of herds
  int<lower=1> K_1;  // nbr of REs
  vector[N] Z_1;  // RE design matrix
  // Random effects for cows within herds
  int<lower=1> J_2[N];  // cow RE levels
  int<lower=1> N_2;  // nbr of cows
  int<lower=1> K_2;  // nbr of REs
  vector[N] Z_2;  // RE design matrix
}
parameters {
  vector[K] b;  // fixed effects
  real temp_Intercept;  // temporary intercept
  real<lower=0> sd_1;  // RE std dev (sigma)
  vector[N_1] z_1;  // unscaled REs (alpha)
  real<lower=0> sd_2;  // RE std dev (sigma)
  vector[N_2] z_2;  // unscaled REs (alpha)
}
transformed parameters{
  vector[N] eta;  // linear predictor
  vector[N_1] r_1;  // REs
  vector[N_2] r_2;  // REs
  // compute linear predictor
  eta <- X * b + temp_Intercept;
  r_1 <- sd_1 * (z_1);  // scale REs
  r_2 <- sd_2 * (z_2);  // scale REs
  for (n in 1:N) {
    eta[n] <- eta[n] + r_1[J_1[n]] * Z_1[n] + r_2[J_2[n]] * Z_2[n];
  }
}
model {
  // priors
  sd_1 ~ student_t(3, 0, 10);
  z_1 ~ normal(0, 1);
  sd_2 ~ student_t(3, 0, 10);
  z_2 ~ normal(0, 1);

  // likelihood contribution
     y ~ bernoulli_logit(eta);
}
generated quantities {
  real b_Intercept;  // herd-level intercept
  b_Intercept <- temp_Intercept - dot_product(X_means, b);
}
