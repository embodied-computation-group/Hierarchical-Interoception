
Model_recovery_big = function(parameters){
  
  
  data <- read_csv(here::here("model_parameter_recovery","matlabdata","big",paste0("R_",parameters$model,"_50_80_all_data.csv"))) %>% 
    filter(idx == parameters$id)
  
  
  # getting mixed and non centered parameterization for the 4 models:
  model = "gaussian"
  mod_gaus_mixed = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("mixed_param.stan")))
  mod_gaus_noncent = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("full_noncentered.stan")))
  
  
  model = "hyper"
  mod_hyper_mixed = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("mixed_param.stan")))
  mod_hyper_noncent = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("full_noncentered.stan")))
  
  model = "gumbel"
  mod_gumbel_mixed = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("mixed_param.stan")))
  mod_gumbel_noncent = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("full_noncentered.stan")))
  
  model = "logit"
  mod_logit_mixed = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("mixed_param.stan")))
  mod_logit_noncent = cmdstanr::cmdstan_model(here::here("model_parameter_recovery","stanmodels",model,paste0("full_noncentered.stan")))
  
  ############ fitting normal
  normal_fit1 = fit_model(data,mod_gaus_mixed, model = "gaussian")
  
  if(normal_fit1[[2]]$mean_div != 0){
    normal_fit2 = fit_model(data,mod_gaus_noncent, model = "gaussian")
    if(normal_fit1[[2]]$mean_div < normal_fit2[[2]]$mean_div){
      normal_fit = normal_fit1
    }else{
      normal_fit = normal_fit2
    }
  }else{
    normal_fit = normal_fit1
  }
  
  normal_param = extract_parameters(normal_fit[[3]], model = "gaussian", id = parameters$id)
  
  
  ############ fitting hyperbolicf
  hyper_fit1 = fit_model(data,mod_hyper_mixed, model = "hyper")
  
  if(hyper_fit1[[2]]$mean_div != 0){
    hyper_fit2 = fit_model(data,mod_hyper_noncent, model = "hyper")
    
    if(hyper_fit1[[2]]$mean_div < hyper_fit2[[2]]$mean_div){
      hyper_fit = hyper_fit1
    }else{
      hyper_fit = hyper_fit2
    }
  }else{
    hyper_fit = hyper_fit1
  }
  
  hyper_param = extract_parameters(hyper_fit[[3]], model = "hyper", id = parameters$id)
  
  
  
  ############ fitting gumbel
  gumbel_fit1 = fit_model(data,mod_gumbel_mixed, model = "gumbel")
  
  
  if(gumbel_fit1[[2]]$mean_div != 0){
    gumbel_fit2 = fit_model(data,mod_gumbel_noncent, model = "gumbel")
    if(gumbel_fit1[[2]]$mean_div < gumbel_fit2[[2]]$mean_div){
      gumbel_fit = gumbel_fit1
    }else{
      gumbel_fit = gumbel_fit2
    }
  }else{
    gumbel_fit = gumbel_fit1
  }
  
  gumbel_param = extract_parameters(gumbel_fit[[3]], model = "gumbel", id = parameters$id)
  
  
  ############ fitting logistic
  logistic_fit1 = fit_model(data,mod_logit_mixed, model = "logit")
  
  if(logistic_fit1[[2]]$mean_div != 0){
    logistic_fit2 = fit_model(data,mod_logit_noncent, model = "logit")
    
    if(logistic_fit1[[2]]$mean_div < logistic_fit2[[2]]$mean_div){
      logistic_fit = logistic_fit1
    }else{
      logistic_fit = logistic_fit2
    }
  }else{
    logistic_fit = logistic_fit1
  }
  
  logistic_param = extract_parameters(logistic_fit[[3]], model = "logit", id = parameters$id)
  
  
  
  #compare models
  
  loo = loo::loo_compare(list(gaussian = normal_fit[[1]],
                              hyperbolic = hyper_fit[[1]],
                              gumbel = gumbel_fit[[1]],
                              logit = logistic_fit[[1]]))
  
  #combine diagnositics
  
  fit_diagnostics = rbind(normal_fit[[2]],
                          hyper_fit[[2]],
                          gumbel_fit[[2]],
                          logistic_fit[[2]]) %>% mutate(replicate = parameters$replicate,
                                                        subjects = parameters$subjects,
                                                        trials = parameters$trials,
                                                        dataset = parameters$model,
                                                        idx = parameters$id)
  
  
  #combine model comparison
  result_modelcompar = inner_join(data.frame(loo) %>% rownames_to_column(var = "models") %>% 
                                    dplyr::select(models, elpd_diff,se_diff) %>% mutate(elpd_ratio = elpd_diff/se_diff), fit_diagnostics) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model,
           idx = parameters$id)
  
  subj_level_param = rbind(normal_param[[1]],hyper_param[[1]],gumbel_param[[1]],logistic_param[[1]]) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model,
           idx = parameters$id)
  
  group_level_param = rbind(normal_param[[2]],hyper_param[[2]],gumbel_param[[2]],logistic_param[[2]]) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model,
           idx = parameters$id)
  
  real_subj_level_param = rbind(normal_param[[3]],hyper_param[[3]],gumbel_param[[3]],logistic_param[[3]]) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model,
           idx = parameters$id)
  
  
  
  return(list(fit_diagnostics,result_modelcompar,subj_level_param,group_level_param,real_subj_level_param))
  
}




fit_model = function(data,mod, model){
  
  
  
  datastan = list(Y = data$response,
                  N = nrow(data),
                  S = length(unique(data$participant)),
                  S_id = data$participant,
                  X = matrix(c(rep(1,nrow(data)), data$stimuli), ncol = 2, nrow = nrow(data)))
  
  
  #fitting
  fit <- mod$sample(
    data = datastan,
    iter_sampling = 1000,
    iter_warmup = 1000,
    chains = 4,
    # refresh = 0,
    init = 0,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  
  diags_norm = data.frame(fit$diagnostic_summary()) %>% 
    mutate(models = model) %>% group_by(models) %>% summarize(mean_div = mean(num_divergent),
                                                              mean_treedepth = mean(num_max_treedepth))
  
  rhat_norm = data.frame(fit$summary(c("gm[1]","gm[2]","gm[3]",
                                       "tau_u[1]","tau_u[2]","tau_u[3]"))) %>% 
    mutate(models = model) %>% group_by(models) %>% summarize(meanrhat = mean(rhat))
  
  
  
  options(mc.cores = 4)
  
  fit_loo = loo(fit$draws("log_lik"),moment_match=T)
  
  # fit_loomm = fit$loo(moment_match = T)
  # 
  # fit_loo = fit$loo()
  
  normal_loo_diag = sum(fit_loo$diagnostics$pareto_k > 0.7)
  
  diags = inner_join(diags_norm,rhat_norm) %>% mutate(pareto_k_over_0.7 = normal_loo_diag)
  
  return(list(fit_loo, diags,fit))
}


extract_parameters = function(fit, model, id){
  
  parameters = c("alpha","beta","lapse")
  
  subj_levelparam = fit$summary(parameters)
  
  subj_levelparam = subj_levelparam %>% mutate(id = gsub(".*\\[(\\d+)\\].*", "\\1", variable)) %>% 
    mutate(variable = ifelse(grepl("alpha",variable),"threshold",ifelse(grepl("beta",variable),"slope","lapse")))%>% 
    mutate(model = model, id = id)
  
  group_levelparam = fit$summary(c("gm[1]","gm[2]","gm[3]","tau_u[1]","tau_u[2]","tau_u[3]"))%>% 
    mutate(model = model, id = id)
  
  group_levelparam$variable = c("mu_threshold","mu_slope","mu_lapse","sd_threshold","sd_slope","sd_lapse") 
  
  real_subj_param = read.csv(here::here("model_parameter_recovery","matlabdata",paste0("R_",model,"_50_80_all_participant_information.csv"))) %>% 
    filter(idx == id)
  
  real_group_param = read.csv(here::here("model_parameter_recovery","matlabdata",paste0("R_",model,"_50_80_all_dataset_information.csv")))%>% 
    filter(idx == id) %>% mutate(X = NULL) %>% 
    pivot_longer(cols = c("mu_alpha","mu_beta","mu_lambda","tau_alpha","tau_beta","tau_lambda"), values_to = "simulated_values", names_to = "variable")
  
  
  group_levelparam$simulated_values = real_group_param$simulated_values
  
  return(list(subj_esti = subj_levelparam,group = group_levelparam, real = real_subj_param))
  
}





