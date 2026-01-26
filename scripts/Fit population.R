
### HRDT

run_model = function(){
  
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr)
  
  source(here::here("Fitting population","R scripts", "Recover_data.R"))
  
  raw_hrd = read.csv(here::here("Fitting population","Data","raw_hrd.csv"))
  
  hrd_data = prep_data(raw_hrd, ses = 1) 
  
  
  hrd_data = hrd_data %>% mutate(nTrials = nTrials+1,
                                 Modality = ifelse(Modality == "Extero",1,ifelse(Modality == "Intero",0,NA)))
  
  df = hrd_data
  
  # model = "Gumbel"
  model = "Normal"
  
  mod_norm_prior = cmdstanr::cmdstan_model(here::here("stanmodels","population fit",paste0(model,"_refit.stan")),stanc_options = list("O1"))
  
  
  print(paste0("Running! ", model,": Extero_priors_noncentered"))
  df = df %>% arrange(Modality,s)
  
  datastan = list(T = nrow(df),
                  S = length(unique(df$s)),
                  S_id = as.numeric(df$s),
                  X = df %>% .$x,
                  X_lapse = as.matrix(t(data.frame(int = rep(1,nrow(df))))),
                  X_alpha = as.matrix(t(data.frame(int = rep(1,nrow(df)),
                                                   session = df %>% .$Modality))),
                  X_beta = as.matrix(t(data.frame(int = rep(1,nrow(df)),
                                                  session = df %>% .$Modality))),
                  N_alpha = 2,
                  N_beta = 2,
                  N_lapse = 1,
                  npx = rep(1,nrow(df)),
                  Y = df %>% .$y
  )
  
  
  fit_norm <- mod_norm_prior$sample(
    data = datastan,
    iter_sampling = 2000,
    iter_warmup = 2000,
    chains = 4,
    init = 0,
    parallel_chains = 4,
    refresh = 10,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
  fit_norm$save_object(file.path(here::here("Fitting population", "Saved models","Extero","new_priors",paste0("Extero_",model,".RDS"))))
  
  
}

run_model()


### RRST

raw_rrst = read.csv(here::here("..","..","..","mnt","slow_scratch","H-intero","analyses","fit population","RRST","data","raw_rrst.csv"))

ACC = raw_rrst %>% group_by(id) %>% summarize(Acc = sum(ResponseCorrect) / n()) 

grouplevel  = raw_rrst %>% 
  group_by(id) %>% summarize(Acc = sum(ResponseCorrect) / n()) %>% 
  summarize(mean = summary(Acc)[[4]],IQR = IQR(Acc), median = summary(Acc)[[3]], q1 = summary(Acc)[[2]], q3 = summary(Acc)[[5]])


badids = ACC %>% 
  mutate(q1 = grouplevel$q1, q3 = grouplevel$q3, IQR = grouplevel$IQR) %>% 
  filter(!(Acc > q1 - 1.5*IQR & Acc < q3 + 1.5*IQR))



raw_rrst %>% filter(id == 263) %>% mutate(Confidence = as.numeric(Confidence)) %>% 
  pivot_longer(cols = c("ResponseCorrect","DecisionRT","Confidence")) %>% 
  ggplot(aes(x = StimulusLevel, y = value, col = as.factor(Decision)))+
  facet_wrap(~name, ncol = 1, scales = "free")+
  geom_jitter(width = 0.05, height = 0.1)+geom_smooth()

raw_rrst %>% mutate(Confidence = as.numeric(Confidence)) %>% 
  pivot_longer(cols = c("ResponseCorrect","DecisionRT","Confidence")) %>% group_by(StimulusLevel,name) %>% 
  summarize(mean_resp = mean(value, na.rm = T), sd_resp = sd(value, na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = StimulusLevel, y = mean_resp))+
  facet_wrap(~name, ncol = 1, scales = "free")+geom_smooth()+
  geom_pointrange(aes(x = StimulusLevel, y = mean_resp, ymin = mean_resp - 2 * sd_resp, ymax = mean_resp + 2 * sd_resp))


data = raw_rrst %>% mutate(Confidence = as.numeric(Confidence),
                           DecisionRT = as.numeric(DecisionRT)) %>% 
  mutate(Confidence  = Confidence/100) %>% 
  filter(DecisionRT > 0.1) %>% drop_na()%>% 
  arrange(id) %>% select(id,participant_id, StimulusLevel,DecisionRT, Confidence, ResponseCorrect)

#data = data %>% filter(!id %in% bad_ids$id)%>%mutate(id = as.numeric(as.factor(id)))%>%arrange(id)
data = data %>% filter(!id %in% badids$id)%>%mutate(id = as.numeric(as.factor(id)))%>%arrange(id)


Xx = ifelse(log(data %>% .$StimulusLevel) == -Inf, 999, log(data %>% .$StimulusLevel))


datastan = list(T = nrow(data),
                S = length(unique(data$id)),
                S_id = as.numeric(data$id),
                X = Xx,
                Y = data %>% .$ResponseCorrect
)

mod = cmdstanr::cmdstan_model(here::here("stanmodels","population fit","RRST.stan"))

fit_norm <- mod$sample(
  data = datastan,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  # init = 0,
  parallel_chains = 4,
  refresh = 10,
  adapt_delta = 0.90,
  max_treedepth = 12
)

fit_norm$save_object(here::here("results","Fit population","RRST_gumbelpriors_n_centered_3000_0.9.RDS"))
