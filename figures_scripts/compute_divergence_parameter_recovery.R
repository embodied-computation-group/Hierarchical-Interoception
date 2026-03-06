#Script used to check the proportion of datasets rejected from the recovery analysis due to divergent transitions
#Author: Arthur S. Courtin
#Licence: MIT
library(tidyverse)

gu<-readRDS(here::here("results","recovery","Gumbel_datasets.RData"))
ga<-readRDS(here::here("results","recovery","Gaussian_datasets.RData"))

n_rows <- 200

recovery_analysis_divergence <- tibble(
  generative = c(rep('gaussian',100),rep('gumbell',100)),
  sum_div = rep(0,200),
  init_error = rep(0,200)
)

for(idx in 1:100){
  if(is.list(ga[[idx]])){
      mr<-ga[idx][[1]][[2]] %>% 
        filter(variable=='gm[1]')
    
    recovery_analysis_divergence$sum_div[idx]<-sum(mr$mean_div)>0
  }else{
    recovery_analysis_divergence$init_error[idx]<-1
  }
  if(is.list(gu[[idx]])){
    mr<-gu[idx][[1]][[2]] %>% 
      filter(variable=='gm[1]')
    
    recovery_analysis_divergence$sum_div[idx+100]<-sum(mr$mean_div)>0
  }else{
    recovery_analysis_divergence$init_error[idx]<-1
  }
}


recovery_analysis_divergence_summary<-
  recovery_analysis_divergence%>% 
  group_by(generative) %>% 
  summarise(n_with_error=sum(init_error),n_with_div=sum(sum_div>0))


print(recovery_analysis_divergence_summary)



