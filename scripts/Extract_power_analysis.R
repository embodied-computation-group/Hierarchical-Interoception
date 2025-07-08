pacman::p_load(cmdstanr, tidyverse,posterior, tidybayes, furrr)



hier_power = function(){
  
  # files = list.files(here::here("..","..","..","mnt","slow_scratch","H-intero","analyses","poweranalysis","results","Threshold"), full.names = T)
  files = list.files(here::here("results", "power analysis","Threshold"), full.names = T)
  
  extract_data <- function(file_name) {
    pattern <- "subject_(\\d+)_trials_(\\d+)_effectsize_(\\d+\\.\\d+)"
    matches <- regmatches(file_name, regexec(pattern, file_name))[[1]]
    data.frame(
      subject = as.numeric(matches[2]),
      trials = as.numeric(matches[3]),
      effect_size = as.numeric(matches[4])
    )
  }
  
  
  grouplevel = data.frame()
  sublevel = data.frame()
  p_val = data.frame()
  
  for (file in files){
    qq = readRDS(file)  
    sim_data = extract_data(file)
    
    l_grouplevel = map_dfr(qq,1) %>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    l_sublevel = map_dfr(qq,2)%>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    l_p_val = data.frame(p_val = as.numeric(map_dfr(qq,3))) %>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    
    
    grouplevel = rbind(grouplevel,l_grouplevel)
    sublevel = rbind(sublevel,l_sublevel)
    p_val = rbind(p_val,l_p_val)
  }
  
  # for the threshold:
  
  p_val_threshold = grouplevel %>% filter(variable == "alpha_dif" & parameter == "mu") %>% 
    mutate(signi = ifelse(mean - 1.96 * sd > 0 | mean + 1.96 * sd < 0, 1, 0))
  
  p_val_threshold %>% group_by(subjects,trials,effect_size) %>% summarize(mean = sum(signi)) %>% 
    ggplot(aes(x = effect_size, y = mean, ymin = mean, ymax = mean))+
    geom_pointrange()+
    facet_grid(trials~subjects, labeller = label_both)
  
  p_val %>% group_by(subjects,trials,effect_size) %>% summarize(mean = sum(p_val<0.05)) %>% 
    ggplot(aes(x = effect_size, y = mean, ymin = mean, ymax = mean))+
    geom_pointrange()+
    facet_grid(trials~subjects, labeller = label_both)
  ## error bars:
  
  data_alpha = p_val_threshold %>% group_by(subjects,trials,effect_size) %>% 
    summarize(
      sum_signi = sum(signi),   # Number of significant cases
      n = n(),                  # Total observations
      alpha_post = 1 + sum_signi,   # Updated Beta posterior alpha
      beta_post = 1 + (n - sum_signi), # Updated Beta posterior beta
      mean = sum_signi / n,   # Proportion of significant cases
      lower = qbeta(0.025, alpha_post, beta_post),  # 95% HDI lower bound
      upper = qbeta(0.975, alpha_post, beta_post)   # 95% HDI upper bound
    ) 
  
  
  data_alpha%>% 
    ggplot(aes(x = effect_size, y = mean, ymin = lower, ymax = upper))+
    geom_pointrange()+
    facet_grid(trials~subjects, labeller = label_both)+
    geom_hline(yintercept = 0.05)
  
  
  # modeling:
  
  p_val_threshold_data = p_val_threshold %>% group_by(subjects,trials,effect_size) %>% 
    summarize(
      sum_signi = sum(signi),   # Number of significant cases
      n = n())
  
  
  
  # mod = cmdstanr::cmdstan_model(here::here("Power analysis", "analysis","stanmodels","baseline.stan"),stanc_options = list("O1"))
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","poweranalysis","baseline.stan"),stanc_options = list("O1"))
  
  standata = list(x = p_val_threshold_data$effect_size,
                  N = nrow(p_val_threshold_data),
                  design_matrix = data.frame(subs = p_val_threshold_data$subjects,
                                             trials = p_val_threshold_data$trials),
                  param = 2,
                  ns = p_val_threshold_data$n,
                  y = p_val_threshold_data$sum_signi)
  
  
  
  fit_baseline <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    init = 0,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )

  
  # for the slope:
  
  files = list.files(here::here("results", "power analysis","Slope"), full.names = T)
  
  grouplevel_slope = data.frame()
  sublevel_slope = data.frame()
  p_val_slope = data.frame()
  
  for (file in files){
    qq = readRDS(file)  
    sim_data = extract_data(file)
    
    l_grouplevel = map_dfr(qq,1) %>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    l_sublevel = map_dfr(qq,2)%>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    l_p_val = data.frame(p_val = as.numeric(map_dfr(qq,3))) %>% mutate(subjects = sim_data$subject, trials = sim_data$trials, effect_size = sim_data$effect_size)
    
    
    grouplevel_slope = rbind(grouplevel_slope,l_grouplevel)
    sublevel_slope = rbind(sublevel_slope,l_sublevel)
    p_val_slope = rbind(p_val_slope,l_p_val)
  }
  
  
  p_val_slope = grouplevel_slope %>% filter(variable == "beta_dif" & parameter == "mu") %>% 
    mutate(signi = ifelse(mean - 1.96 * sd > 0 | mean + 1.96 * sd < 0, 1, 0))
  
  
  
  ## error bars:
  
  data_beta = p_val_slope %>% group_by(subjects,trials,effect_size) %>% 
    summarize(
      sum_signi = sum(signi),   # Number of significant cases
      n = n(),                  # Total observations
      alpha_post = 1 + sum_signi,   # Updated Beta posterior alpha
      beta_post = 1 + (n - sum_signi), # Updated Beta posterior beta
      mean = sum_signi / n,   # Proportion of significant cases
      lower = qbeta(0.025, alpha_post, beta_post),  # 95% HDI lower bound
      upper = qbeta(0.975, alpha_post, beta_post)   # 95% HDI upper bound
    ) 
  
  
  # data_beta%>% 
  #   ggplot(aes(x = effect_size, y = mean, ymin = lower, ymax = upper))+
  #   geom_pointrange()+
  #   facet_grid(trials~subjects, labeller = label_both)+
  #   geom_hline(yintercept = 0.05)
  
  
  # modeling:
  
  p_val_slope_data = p_val_slope %>% group_by(subjects,trials,effect_size) %>% 
    summarize(
      sum_signi = sum(signi),   # Number of significant cases
      n = n())
  
  
  standata = list(x = p_val_slope_data$effect_size,
                  N = nrow(p_val_slope_data),
                  design_matrix = data.frame(subs = p_val_slope_data$subjects,
                                             trials = p_val_slope_data$trials),
                  param = 2,
                  ns = p_val_slope_data$n,
                  y = p_val_slope_data$sum_signi)
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","baseline.stan"),stanc_options = list("O1"))

  
  fit_baseline_slope <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    init = 0,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  data = rbind(data_beta%>% mutate(parameter = "Slope"), data_alpha %>% mutate(parameter = "Threshold"))
  
  
  
  data = data %>% ungroup() %>% group_by(parameter) %>% mutate(s_subject = scale(subjects)[,1],
                                                               s_trial = scale(trials)[,1],
                                                               s_effectsize = scale(effect_size)[,1]) %>% 
    ungroup()
  
  
  
  
  
  return(list(data,fit_baseline,fit_baseline_slope))
  
  
}




psi_power = function(){
  
  
  files = list.files(here::here("simulated_data","poweranalysis"), full.names = T)
  
  
  files = files[grepl("ppt_info", files)]
  
  extract_info <- function(filepath) {
    match <- stringr::str_match(basename(filepath), "^(Slope|Threshold)_([0-9\\.]+)_")[,2:3]
    return(data.frame(Parameter = match[1], EffectSize = as.numeric(match[2]), stringsAsFactors = FALSE))
  }
  
  data_psi = data.frame()
  
  for (file in files){
    qq = read.csv(file)
    parameter = extract_info(file)[1,1]
    effect_size = extract_info(file)[1,2]
    for(sub in c(15,30,60,120)){
      for(t in c(30,60,90)){
        # Filter participants
        df <- qq %>% filter(participant_idx %in% 1:sub)
        
        if(parameter == "Threshold"){
          # Dynamically select the correct column names
          control_col <- paste0("alpha_estimate_control_", t)
          treatment_col <- paste0("alpha_estimate_treatment_", t)
        }else if(parameter == "Slope"){
          control_col <- paste0("beta_estimate_control_", t)
          treatment_col <- paste0("beta_estimate_treatment_", t)
        }
        # Ensure columns exist before proceeding
        if (!(control_col %in% colnames(df) && treatment_col %in% colnames(df))) {
          next  # Skip iteration if the expected columns are not present
        }
        
        
        # Summarize p-value per dataset
        df_summary <- df %>% 
          group_by(dataset_idx) %>% 
          summarize(ttest = t.test(!!sym(control_col), !!sym(treatment_col))$p.value, .groups = "drop")
        
        # Add metadata
        df_summary$trials <- t
        df_summary$subjects <- sub
        df_summary$effect_size = effect_size
        df_summary$parameter = parameter
        
        data_psi = rbind(data_psi,df_summary)
        
      }
    }
  }
  
  
  data_psi = data_psi  %>% 
    group_by(trials,subjects,effect_size,parameter) %>% 
    summarize(
      sum_signi = sum(ttest < 0.05),   # Number of significant cases
      n = n(),                  # Total observations
      alpha_post = 1 + sum_signi,   # Updated Beta posterior alpha
      beta_post = 1 + (n - sum_signi), # Updated Beta posterior beta
      mean = sum_signi / n,   # Proportion of significant cases
      lower = qbeta(0.025, alpha_post, beta_post),  # 95% HDI lower bound
      upper = qbeta(0.975, alpha_post, beta_post)   # 95% HDI upper bound
    ) 
  
  
  
  
  
  # modeling:
  
  p_val_threshold_data = data_psi %>% filter(parameter == "Threshold") 
  
  
  standata = list(x = p_val_threshold_data$effect_size,
                  N = nrow(p_val_threshold_data),
                  design_matrix = data.frame(subs = p_val_threshold_data$subjects,
                                             trials = p_val_threshold_data$trials),
                  param = 2,
                  ns = p_val_threshold_data$n,
                  y = p_val_threshold_data$sum_signi)
  
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","baseline.stan"),stanc_options = list("O1"))
  
  
  fit_psi <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  
  
  
  
  # modeling:
  
  p_val_slope_data = data_psi %>% filter(parameter == "Slope") 
  
  
  
  
  standata = list(x = p_val_slope_data$effect_size,
                  N = nrow(p_val_slope_data),
                  design_matrix = data.frame(subs = p_val_slope_data$subjects,
                                             trials = p_val_slope_data$trials),
                  param = 2,
                  ns = p_val_slope_data$n,
                  y = p_val_slope_data$sum_signi)
  
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","baseline.stan"),stanc_options = list("O1"))

  
  fit_psi_slope <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  
  data_psi = data_psi %>% ungroup() %>% group_by(parameter) %>% mutate(s_subject = scale(subjects)[,1],
                                                                       s_trial = scale(trials)[,1],
                                                                       s_effectsize = scale(effect_size)[,1]) %>% 
    ungroup()
  
  
  
  
  return(list(data_psi,fit_psi,fit_psi_slope))
  
  
}



psi_power_unc = function(){
  
  # files = list.files(here::here("..","..","..","mnt","slow_scratch","H-intero","analyses","poweranalysis","new_matlabdata"), full.names = T)
  files = list.files(here::here("simulated_data","poweranalysis"), full.names = T)
  
  files = files[grepl("ppt_info", files)]
  
  extract_info <- function(filepath) {
    match <- stringr::str_match(basename(filepath), "^(Slope|Threshold)_([0-9\\.]+)_")[,2:3]
    return(data.frame(Parameter = match[1], EffectSize = as.numeric(match[2]), stringsAsFactors = FALSE))
  }
  
  
  data_psi_unc = data.frame()
  
  for (file in files){
    qq = read.csv(file)
    parameter = extract_info(file)[1,1]
    effect_size = extract_info(file)[1,2]
    for(sub in c(15,30,60,120)){
      for(t in c(30,60,90)){
        # Filter participants
        df <- qq %>% filter(participant_idx %in% 1:sub)
        
        if(parameter == "Threshold"){
          # Dynamically select the correct column names
          control_col <- paste0("alpha_estimate_control_", t)
          control_se = paste0("alpha_se_control_",t)
          treatment_col <- paste0("alpha_estimate_treatment_", t)
          treatment_se = paste0("alpha_se_treatment_",t)
        }else if(parameter == "Slope"){
          control_col <- paste0("beta_estimate_control_", t)
          control_se <- paste0("beta_se_control_", t)
          treatment_col <- paste0("beta_estimate_treatment_", t)
          treatment_se <- paste0("beta_se_treatment_", t)
        }
        # Ensure columns exist before proceeding
        if (!(control_col %in% colnames(df) && treatment_col %in% colnames(df))) {
          next  # Skip iteration if the expected columns are not present
        }
        
        
        # Summarize p-value per dataset
        df_summary <- df %>% 
          group_by(dataset_idx) %>% 
          summarize(ttest = weights::wtd.t.test(!!sym(control_col), !!sym(treatment_col), weight = 1/(!!sym(control_se))^2,weighty = 1/(!!sym(treatment_se))^2)$coefficients[3])
        
        # Add metadata
        df_summary$trials <- t
        df_summary$subjects <- sub
        df_summary$effect_size = effect_size
        df_summary$parameter = parameter
        
        data_psi_unc = rbind(data_psi_unc,df_summary)
        
      }
    }
  }
  
  
  
  data_psi_unc = data_psi_unc %>% 
    group_by(trials,subjects,effect_size,parameter) %>% 
    summarize(
      sum_signi = sum(ttest < 0.05),   # Number of significant cases
      n = n(),                  # Total observations
      alpha_post = 1 + sum_signi,   # Updated Beta posterior alpha
      beta_post = 1 + (n - sum_signi), # Updated Beta posterior beta
      mean = sum_signi / n,   # Proportion of significant cases
      lower = qbeta(0.025, alpha_post, beta_post),  # 95% HDI lower bound
      upper = qbeta(0.975, alpha_post, beta_post)   # 95% HDI upper bound
    ) 
  
  
  
  
  
  
  
  
  p_val_threshold_data = data_psi_unc %>% filter(parameter == "Threshold") 
  
  
  standata = list(x = p_val_threshold_data$effect_size,
                  N = nrow(p_val_threshold_data),
                  design_matrix = data.frame(subs = p_val_threshold_data$subjects,
                                             trials = p_val_threshold_data$trials),
                  param = 2,
                  ns = p_val_threshold_data$n,
                  y = p_val_threshold_data$sum_signi)
  
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","baseline.stan"),stanc_options = list("O1"))

  
  fit_psi <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  # modeling:
  
  p_val_slope_data = data_psi_unc %>% filter(parameter == "Slope") 
  
  
  standata = list(x = p_val_slope_data$effect_size,
                  N = nrow(p_val_slope_data),
                  design_matrix = data.frame(subs = p_val_slope_data$subjects,
                                             trials = p_val_slope_data$trials),
                  param = 2,
                  ns = p_val_slope_data$n,
                  y = p_val_slope_data$sum_signi)
  
  mod = cmdstanr::cmdstan_model(here::here("stanmodels","baseline.stan"),stanc_options = list("O1"))
  
  
  fit_psi_slope <- mod$sample(
    data = standata,
    chains = 4,
    refresh = 50,
    iter_warmup = 1000,
    iter_sampling = 1000,
    parallel_chains = 4,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  
  data_psi_unc = data_psi_unc %>% ungroup() %>% group_by(parameter) %>% mutate(s_subject = scale(subjects)[,1],
                                                                               s_trial = scale(trials)[,1],
                                                                               s_effectsize = scale(effect_size)[,1]) %>% 
    ungroup()
  
  
  return(list(data_psi_unc,fit_psi,fit_psi_slope))
  
  
}



hier = hier_power()

write.csv(hier[[1]], here::here("Extracted","Hierarchical","df.csv"))
hier[[2]]$save_object(here::here("Extracted","Hierarchical","thresholdmodel.rds"))
hier[[3]]$save_object(here::here("Extracted","Hierarchical","slopemodel.rds"))

psi = psi_power()

write.csv(psi[[1]], here::here("Extracted","Psi","df.csv"))
psi[[2]]$save_object(here::here("Extracted","Psi","thresholdmodel.rds"))
psi[[3]]$save_object(here::here("Extracted","Psi","slopemodel.rds"))



psi_un = psi_power_unc()
write.csv(psi_un[[1]], here::here("Extracted","Psi_un","df.csv"))
psi_un[[2]]$save_object(here::here("Extracted","Psi_un","thresholdmodel.rds"))
psi_un[[3]]$save_object(here::here("Extracted","Psi_un","slopemodel.rds"))



