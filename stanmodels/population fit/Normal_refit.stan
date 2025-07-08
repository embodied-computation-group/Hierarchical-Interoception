data{
  //Constants
  int<lower=1> T; //n Trials
  int<lower=1> S; //n participants
  array[T] int S_id; //n participants
  
  int<lower=1> N_alpha;
  int<lower=1> N_beta;
  int<lower=1> N_lapse;
  

  matrix[N_alpha,T] X_alpha;
  matrix[N_beta,T] X_beta;
  matrix[N_lapse,T] X_lapse;
  
  
  array[T] int Y;
  vector[T] X;
  
  array[T] int npx;  // total number of observations per X

  

}
transformed data{
  int<lower=1> N=N_beta+N_lapse+N_alpha; //n free parameters
  int<lower=1> N_centered=N_beta+N_alpha; //n free parameters

}
parameters{
  // Group means 
  vector[N] gm;
  // Between participant scales
  vector<lower = 0>[N]  tau_u;
  // Between participant cholesky decomposition
  cholesky_factor_corr[N_centered] L_u;
  // Participant deviation 
  matrix[N_lapse, S] z_expo;  
  
  array[S] vector[N_centered] param;        // Centered individual parameters


}
transformed parameters{
  matrix[N_lapse,S] lapse_p =  gm[N] + z_expo * tau_u[N];
  

    
}

model{


  gm[1] ~ normal(-1.3,1.3); //global mean of beta
  
  gm[2] ~ normal(0, 3); //global difference beta
  
  gm[3] ~ normal(0,50); //global threshold

  gm[4] ~ normal(0,50); //global difference threshold

  gm[5] ~ normal(-4,2); //global mean of lapse


  to_vector(z_expo) ~ std_normal();
  
  tau_u[1] ~ normal(0 , 1.3);
  
  tau_u[2] ~ normal(0 , 1.3);
  
  tau_u[3] ~ normal(0 , 50);
  
  tau_u[4] ~ normal(0 , 50);
  
  tau_u[5] ~ normal(0 , 2);
    
  
  
  param ~ multi_normal_cholesky(gm[1:N_centered], diag_pre_multiply(tau_u[1:N_centered], L_u)); 

  L_u ~ lkj_corr_cholesky(2);


  ///Recomposition

  vector[T] alpha;
  vector[T] beta;
  vector[T] lapse;
  
  
  
  for(n in 1:T){

    beta[n] = exp((X_beta[1,n] * param[S_id[n],1] + X_beta[2,n] * param[S_id[n],2]));
    
    alpha[n] = X_alpha[1,n] * param[S_id[n],3] + X_alpha[2,n] * param[S_id[n],4];
    
    lapse[n] = inv_logit(dot_product(X_lapse[,n], lapse_p[,S_id[n]])) / 2;
    
    }

  Y ~ binomial(npx, lapse + (1 - 2 * lapse) .* (0.5+0.5*erf(beta .* (X-alpha) / sqrt(2))));;
 
}

generated quantities{
  
  matrix[N_centered,N_centered] correlation_matrix;
  correlation_matrix = L_u*L_u';

  vector[T] log_lik;
  

  for (n in 1:T){
    log_lik[n] = binomial_lpmf(Y[n] | npx[n], (inv_logit(dot_product(X_lapse[,n], lapse_p[,S_id[n]])) / 2) + (1 - 2 * (inv_logit(dot_product(X_lapse[,n], lapse_p[,S_id[n]])) / 2)) .* (0.5+0.5*erf((exp((X_beta[1,n] * param[S_id[n],1] + X_beta[2,n] * param[S_id[n],2]))) * (X[n]-(X_alpha[1,n] * param[S_id[n],3] + X_alpha[2,n] * param[S_id[n],4])) / sqrt(2))));  
  }
}
