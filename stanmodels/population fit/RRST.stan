functions{
  
  real psycho(real X, real alpha, real beta, real lapse){
   if(X == 999){
     return(0.5);
   }else{
    return((0.5) + (1 - (0.5) - (inv_logit(lapse) / 2)) .* (1-exp(-10^(exp(beta)*(X-alpha)))));
   }
  
}
}

data{

  //Constants
  int<lower=1> T; // Total number of trials in the data
  int<lower=1> S; // Total number of subjects in the data
  array[T] int S_id; //n vector of integeres that signify participant numbers
  
  array[T] int Y; // Vector of binary responses

  vector[T] X; // Vector of deltaBPM values that match the binary response



}
transformed data{
  int<lower=1> N=3; //Number of free parameters

}
parameters{
  
  vector[N] gm;  // Group means 

  vector<lower = 0>[N]  tau_u;   // Between participant scales

  cholesky_factor_corr[N] L_u;    // Between participant cholesky decomposition

  matrix[N, S] z_expo;    // Participant deviation from the group means
  
}

transformed parameters{
  
  
  // Extracting individual deviations for each subject for each parameter
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  matrix[S, N] param;
  
  // adding the participant deviation to the group means
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
  
  vector[S] beta = param[,1];
  vector[S] lapse = param[,2];
  vector[S] alpha = param[,3];
  
  
  
}

model{
  // Defining priors.

    gm[1] ~ normal(1.2,1.9); //global mean of beta
    
    gm[2] ~ normal(-4, 2); //global mean of lapse
    
    gm[3] ~ normal(-2, 1); //global mean of lapse
    
    // gm_alpha ~ normal(-2,1); //global mean of beta
    
  
    L_u ~ lkj_corr_cholesky(2);
    
    to_vector(z_expo) ~ std_normal();
    
    tau_u[1] ~ normal(0 , 1.9);
    tau_u[2] ~ normal(0 , 2);
    tau_u[3] ~ normal(0 , 1);
    
    // tau_u_alpha ~ normal(0 , 1);
    

    // alpha ~ normal(gm_alpha , tau_u_alpha);
  

  for (n in 1:T){
    
      Y[n] ~ bernoulli(psycho(X[n], alpha[S_id[n]], beta[S_id[n]], lapse[S_id[n]]));
      
  }
  
}


generated quantities{
  
  vector[T] log_lik;
  
  for (n in 1:T){
    
      log_lik[n] = bernoulli_lpmf(Y[n] | psycho(X[n], alpha[S_id[n]], beta[S_id[n]], lapse[S_id[n]]));
      
  }
}

