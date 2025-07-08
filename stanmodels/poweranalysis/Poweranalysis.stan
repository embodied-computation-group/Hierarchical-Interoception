data{

  //Constants
  int<lower=1> T; // Total number of trials in the data
  int<lower=1> S; // Total number of subjects in the data
  array[T] int S_id; //n vector of integeres that signify participant numbers


  vector[T] condition;
  
  array[T] int Y; // Vector of binary responses

  matrix[T, 2] X;  // design matrix (first column being intercept i.e. 1)

  

}
transformed data{
 int N = 5;
}
parameters{
  
  vector[N] gm;  // Group means 

  vector<lower = 0>[N]  tau_u;   // Between participant scales

  matrix[N, S] z_expo;    // Participant deviation from the group means

}
transformed parameters{
  vector[S] alpha_int =(gm[1]+(tau_u[1] * z_expo[1,]))';
  vector[S] beta_int = (gm[2]+(tau_u[2] * z_expo[2,]))';
  vector[S] lapse = (inv_logit(gm[3]+(tau_u[3] * z_expo[3,])) / 2)';
  vector[S] alpha_dif =(gm[4]+(tau_u[4] * z_expo[4,]))';
  vector[S] beta_dif = (gm[5]+(tau_u[5] * z_expo[5,]))';
  
  
    
}

model{
  // Defining priors.


  target += normal_lpdf(gm[1] | 0,50);
  target += normal_lpdf(gm[2] |  -1.3,1.3);
  //target += normal_lpdf(gm_lapse | -4, 2);
  target += normal_lpdf(gm[3] | -4, 2);
  target += normal_lpdf(gm[4] | 0,10);
  target += normal_lpdf(gm[5] | 0,2);
  
  // 
  // // target += normal_lpdf(gm[4] | -3, 2);
  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u[1] | 0, 50);
  target += normal_lpdf(tau_u[2] | 0, 1.3);
  target += normal_lpdf(tau_u[3] | 0, 2);
  
  target += normal_lpdf(tau_u[4] | 0,10);
  target += normal_lpdf(tau_u[5] | 0,2);
  
  // Computing the likelihood. The cummulative normal is used here:
    
  profile("likelihood") {
  for(n in 1:T) {
    Y[n] ~ bernoulli(lapse[S_id[n]] + (1 - 2 * lapse[S_id[n]]) * (0.5+0.5*erf(exp(beta_int[S_id[n]] + beta_dif[S_id[n]] * condition[n])*(X[n,2]-(alpha_int[S_id[n]]+ alpha_dif[S_id[n]] * condition[n])) / sqrt(2))));
  }
  }

}

