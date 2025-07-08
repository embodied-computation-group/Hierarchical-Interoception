data {
  
  int param;

  int<lower=0> N;
  vector[N] x;
  array[N] int y;
  array[N] int ns;
  
  matrix[N,param] design_matrix;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {


  vector [param] expo_alpha;
  vector [param] expo_beta;
  
  // real <lower = 0> int_alpha;
  // real <lower = 0> int_beta;

  // vector [param1] asym_alpha;
  // vector [param1] asym_beta;
  
  
  real intercept_alpha;
  real intercept_beta;
  
}


transformed parameters {
  
  vector[N] X_norm;
  vector[N] alpha;

  vector[N] beta;
  
  real mu_intercept_alpha = exp(intercept_alpha);
  real mu_intercept_beta = exp(intercept_beta);
  
  
  // vector[N] mu_expo_alpha = design_matrix * expo_alpha;
  // vector[N] mu_expo_beta = design_matrix * expo_beta;

  // vector[N] mu_asym_alpha = exp(design_matrix1 * asym_alpha);
  
  // vector[N] mu_asym_beta = exp(design_matrix1 * asym_beta);

  alpha =  mu_intercept_alpha * (design_matrix[,1] ^ expo_alpha[1]) .* (design_matrix[,2] ^ expo_alpha[2]);

  beta =  mu_intercept_beta * (design_matrix[,1] ^ expo_beta[1]) .* (design_matrix[,2] ^ expo_beta[2]);
  
  X_norm = 1/(1+exp(- (1/beta) .* (x - alpha)));

  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  target += normal_lpdf(expo_alpha | -1, 3);
  target += normal_lpdf(expo_beta | -1, 3);
  
  // target += normal_lpdf(int_alpha | 0, 3);
  // target += normal_lpdf(int_beta | 0, 3);
  // 
  target += normal_lpdf(intercept_alpha | 2, 2);
  target += normal_lpdf(intercept_beta | 2, 2);
  
  // target += normal_lpdf(asym_alpha | 0, 3);
  // target += normal_lpdf(asym_beta | 0, 3);
  

  
  
   y ~ binomial(ns, X_norm);
}


generated quantities{
  
  vector[N] log_lik;
  
  for(i in 1:N){
     log_lik[i] = binomial_lpmf(y[i] | ns[i], X_norm[i]);
  }
  
}

