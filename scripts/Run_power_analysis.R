
pacman::p_load(cmdstanr, tidyverse,posterior, tidybayes, furrr)

print("Running")

subjects= 30

parameter = "Slope"
effectsize = "0.5"

trials = 30

matlab_ana = function(parameters){
  
  effectsize = parameters$effectsize
  
  df <- read_csv(here::here("simulated_data","poweranalysis",paste0(parameters$parameter,"_",effectsize,"_120_150_data.csv")))%>% 
    group_by(dataset_idx, participant_idx, condition_is_treatment) %>% 
    mutate(trail = 1:150)%>% 
    filter(dataset_idx == parameters$dataset, participant_idx %in% 1:parameters$subjects, trail %in% 1:parameters$trials)
  
  
  dataset_info <- read_csv(here::here("simulated_data","poweranalysis",paste0(parameters$parameter,"_",effectsize,"_120_150_dataset_info.csv"))) %>% 
    filter(dataset_idx == parameters$dataset)
  
  pp_info <- read_csv(here::here("simulated_data","poweranalysis",paste0(parameters$parameter,"_",effectsize,"_120_150_ppt_info.csv")))%>% 
    filter(dataset_idx == parameters$dataset, participant_idx %in% 1:parameters$subjects)
  
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","poweranalysis","Poweranalysis.stan"))
  
  data = df %>% filter(dataset_idx == parameters$dataset, participant_idx %in% 1:parameters$subjects, trail %in% 1:parameters$trials)
  
  datastan = list(Y = data$response,
                  T = nrow(data),
                  S = length(unique(data$participant_idx)),
                  S_id = data$participant_idx,
                  condition = data$condition_is_treatment,
                  X = matrix(c(rep(1,nrow(data)), data$stimulus), ncol = 2, nrow = nrow(data)))
  
  
  #fitting
  fit <- mod$sample(
    data = datastan,
    iter_sampling = 1000,
    iter_warmup = 1000,
    chains = 4,
    parallel_chains = 4,
    refresh = 0,
    init = 0,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  sim_means = dataset_info %>% select(mu_alpha_int,mu_beta_int, mu_lambda_int, mu_alpha_effect, mu_beta_effect)
  
  sim_taus = dataset_info %>% select(tau_b_alpha_int,tau_b_beta_int,tau_b_lambda_int,tau_b_alpha_effect,tau_b_beta_effect,
                                     tau_w_alpha_int,tau_w_beta_int,tau_w_lambda_int) %>% 
    mutate(tau_alpha = sqrt((tau_b_alpha_int^2 + tau_w_alpha_int^2)),
           tau_beta = sqrt((tau_b_beta_int^2 + tau_w_beta_int^2 )),
           tau_lapse = sqrt((tau_b_lambda_int^2 + tau_w_lambda_int^2 )),
           tau_alpha_dif = tau_b_alpha_effect,
           tau_beta_dif = tau_b_beta_effect,
           
    ) %>% select(tau_alpha, tau_beta, tau_lapse,tau_alpha_dif,tau_beta_dif)
  
  groupmeans = data.frame(fit$summary("gm")) %>% mutate(variable = c("alpha_int","beta_int","lapse","alpha_dif","beta_dif"),
                                                        sum_div = sum(fit$diagnostic_summary()$num_divergent),
                                                        tree_div = sum(fit$diagnostic_summary()$num_max_treedepth),
                                                        sim_n = parameters$id,
                                                        dataset = parameters$dataset,
                                                        estimation_time = fit$time()$total,
                                                        parameter = "mu",
                                                        simulated = as.numeric(sim_means))
  
  groupvar = data.frame(fit$summary("tau_u")) %>% mutate(variable = c("alpha_int","beta_int","lapse","alpha_dif","beta_dif"),
                                                         sum_div = sum(fit$diagnostic_summary()$num_divergent),
                                                         tree_div = sum(fit$diagnostic_summary()$num_max_treedepth),
                                                         sim_n = parameters$id,
                                                         dataset = parameters$dataset,
                                                         estimation_time = fit$time()$total,
                                                         parameter = "tau",
                                                         simulated = as.numeric(sim_taus))
  
  grouplevel = rbind(groupmeans,groupvar)  
  
  
  pp_info = pp_info %>% select(target_parameter,effect_size,dataset_idx,participant_idx,             
                               alpha_control,beta_control,lambda_control,alpha_treatment,             
                               beta_treatment,lambda_treatment) %>% rename(alpha_int = alpha_control,
                                                                           beta_int = beta_control,
                                                                           lapse = lambda_control,
                                                                           lapse2 = lambda_treatment,
                                                                           alpha_dif = alpha_treatment,
                                                                           beta_dif = beta_treatment) %>% 
    pivot_longer(cols = c(alpha_int,beta_int,lapse,alpha_dif,beta_dif,lapse2),names_to = "variable", values_to= "simulated")
  
  subj_level_param = data.frame(fit$summary(c("alpha_int","beta_int","lapse","alpha_dif","beta_dif"))) %>% 
    mutate(participant_idx = as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", variable)),
           variable = gsub("\\[\\d+\\]", "", variable),
           sum_div = sum(fit$diagnostic_summary()$num_divergent),
           tree_div = sum(fit$diagnostic_summary()$num_max_treedepth),
           sim_n = parameters$id,
           dataset = parameters$dataset,
           estimation_time = fit$time()$total)
  
  subj = inner_join(subj_level_param,pp_info)
  
  alpha_draws = as_draws_df(fit$draws("gm[4]")) %>% .$`gm[4]`
  p_zero_alpha = sum(alpha_draws<0)/length(alpha_draws)
  
  return(list(grouplevel,subj,p_zero_alpha))
  
}


print("Running_v2")


parameters = expand.grid(subjects = subjects,
                         parameter = parameter,
                         effectsize = effectsize,
                         trials = trials,
                         dataset = 1:100) %>% 
  mutate(id = 1:nrow(.))

data_list <- split(parameters, parameters$id)

cores = parallelly::availableCores()

plan(multisession, workers = cores / 4)

print("Running_v3")

# results = matlab_ana(data_list[[1]])

possfit_model = possibly(.f = matlab_ana, otherwise = "Error")

results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE))

print("saving")

here::here()


saveRDS(results, here::here("power analysis",parameter,paste0("results_",parameter,"_subject_",subjects,"_trials_",trials,"_effectsize_",effectsize,".rds")))



