#Script used to check the proportion of datasets rejected from the power analysis due to divergent transitions
#Author: Arthur S. Courtin
#Licence: MIT
library(tidyverse)

sample_sizes<-c(15,30,60,120)
trial_numbers<-c(30,60,90)
effect_sizes<-c('0.0','0.2','0.5','0.8')
target_parameters<-c('Threshold','Slope')

n_rows <- 100 *
  length(sample_sizes) *
  length(trial_numbers) *
  length(effect_sizes) *
  length(target_parameters)

power_analysis_divergence <- tibble(
  dataset = integer(n_rows),
  sum_div = numeric(n_rows),
  target_parameter = character(n_rows),
  effect_size = character(n_rows),
  trial_number = integer(n_rows),
  sample_size = integer(n_rows)
)

for(sdx in 1:length(sample_sizes)){
  for(tdx in 1:length(trial_numbers)){
    for(edx in 1:length(effect_sizes)){
      for(pdx in 1:length(target_parameters)){
        fn<-
            paste0(
              'results/power analysis/',target_parameters[pdx],
              '/results_',target_parameters[pdx],
              '_subject_',sample_sizes[sdx],
              '_trials_',trial_numbers[tdx],
              '_effectsize_',effect_sizes[edx],
              '.rds'
              )
        if(file.exists(fn)){
          r<-readRDS(fn)
          for(idx in 1:100){
            d<-r[[idx]][[1]]
            cdx<-
              idx+
              (pdx-1)*100+
              (edx-1)*100*length(target_parameters)+
              (tdx-1)*100*length(target_parameters)*length(effect_sizes)+
              (sdx-1)*100*length(target_parameters)*length(effect_sizes)*length(trial_numbers)
            
            power_analysis_divergence$dataset[cdx]<-idx
            power_analysis_divergence$sum_div[cdx]<-d$sum_div[1]
            power_analysis_divergence$target_parameter[cdx]<-target_parameters[pdx]
            power_analysis_divergence$effect_size[cdx]<-effect_sizes[edx]
            power_analysis_divergence$trial_number[cdx]<-trial_numbers[tdx]
            power_analysis_divergence$sample_size[cdx]<-sample_sizes[sdx]
          } 
        }
      }
    }
  }
}
power_analysis_divergence_summary<-
  power_analysis_divergence%>% 
  filter(target_parameter %in% target_parameters) %>% 
  group_by(target_parameter,effect_size,trial_number,sample_size) %>% 
  summarise(n_with_div=sum(sum_div>0),n_valid=sum(!is.na(sum_div))) %>% 
  rename(n_trials=trial_number,n_participants=sample_size)

power_analysis_divergence_summary %>% 
  filter(target_parameter=='Threshold') %>% 
  ggplot()+
  geom_point(aes(x=as.numeric(effect_size),y=n_with_div/n_valid))+
  geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
  geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
  facet_grid(cols=vars(n_trials),rows=vars(n_participants),labeller = label_both)+
  theme_classic()+
  ylim(c(0,1))+
  labs(
    x='Standardized effect size',
    y='Proportion of datasets with divergence'
    )
  
ggsave('figures/FigureDivThreshold.png',units = 'cm',width = 16,height = 16)

power_analysis_divergence_summary %>% 
  filter(target_parameter=='Slope') %>% 
  ggplot()+
  geom_point(aes(x=as.numeric(effect_size),y=n_with_div/n_valid))+
  geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
  geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
  facet_grid(cols=vars(n_trials),rows=vars(n_participants),labeller = label_both)+
  theme_classic()+
  ylim(c(0,1))+
  labs(
    x='Standardized effect size',
    y='Proportion of datasets with divergence'
  )

ggsave('figures/FigureDivSlope.png',units = 'cm',width = 16,height = 16)





