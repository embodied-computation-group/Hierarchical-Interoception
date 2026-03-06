#Script used to check the proportion of datasets rejected from the recovery analysis due to divergent transitions
#Author: Arthur S. Courtin
#Licence: MIT
library(tidyverse)

gu<-readRDS(here::here("results","recovery","Gumbel_datasets.RData"))
ga<-readRDS(here::here("results","recovery","Gaussian_datasets.RData"))

n_rows <- 200

recovery_analysis_divergence <- tibble(
  generative = c(rep('gaussian',100),rep('gumbell',100)),
  sum_div = numeric(n_rows)
)

for(idx in 1:100){
  mr<-ga[idx][[1]][[2]] %>% 
    filter(variable=='gm[1]')
  
  recovery_analysis_divergence$sum_div[idx]<-sum(mr$mean_div)>0
  
  mr<-gu[idx][[1]][[2]] %>% 
    filter(variable=='gm[1]')
  
  recovery_analysis_divergence$sum_div[idx+100]<-sum(mr$mean_div)>0
}


recovery_analysis_divergence_summary<-
  recovery_analysis_divergence%>% 
  group_by(generative) %>% 
  summarise(n_with_div=sum(sum_div>0),n_valid=sum(!is.na(sum_div)))

recovery_analysis_divergence_summary %>% 
  ggplot()+
  geom_point(aes(x=as.numeric(effect_size),y=n_with_div/n_valid))+
  geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
  geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
  facet_grid(cols=vars(n_trials),rows=vars(n_participants),labeller = label_both)+
  theme_classic()+
  ylim(c(0,1))+
  labs(
    title='Threshold',
    x='Standardized effect size',
    y='Proportion of datasets with divergence'
    )
  
ggsave('figures/FigureRec.png',units = 'cm',width = 16,height = 16)





