library(shiny)
library(tidyverse)
library(flextable)
library(posterior)

get_p = function(df){
  alpha = df$intercept_alpha * (df$subjects ^ df$expo_alpha1) * (df$trials ^ df$expo_alpha2)
  beta = df$intercept_beta * (df$subjects ^ df$expo_beta1) * (df$trials ^ df$expo_beta2)
  
  return(1/(1+exp(- (1/beta) * (df$effect_size - alpha)))
  )
}

all_data <- read_csv(here::here("results","power analysis","Extracted","all_draws.csv"))


ui <- navbarPage("Bayesian Power Simulation Explorer",
                 
                 # --- Page 1: Power grid ---
                 tabPanel("Power Grid",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("parameter", "Parameter:", choices = c("Threshold", "Slope")),
                              numericInput("es", "Effect size (d):", value = 0.5, min = 0, max = 1, step = 0.01),
                              sliderInput("subjects", "N subjects:", min = 15, max = 120, value = c(10, 60), step = 1),
                              sliderInput("trials", "N trials:", min = 30, max = 90, value = c(10, 60), step = 1),
                              radioButtons(
                                inputId = "powerlevel",
                                label = "Desired power level:",
                                choices = c("0.80" = 0.8, "0.90" = 0.9, "0.95" = 0.95),
                                selected = 0.8,
                                inline = TRUE
                              )
                            ),
                            mainPanel(
                              plotOutput("gridPlot"),
                              p(em('Disclaimer: power lines may fall out of the displayed range. If none is visible, try changing the values'))
                            )
                          )
                 ),
                 
                 # --- Page 2: Effect size vs. Power ---
                 tabPanel("Effect Size vs. Power",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("parameter2", "Parameter:", choices = c("Threshold", "Slope")),
                              sliderInput("subjects2", "N subjects:", min = 15, max = 120, value = 30, step = 1),
                              sliderInput("trials2", "N trials:", min = 30, max = 90, value = 30, step = 1),
                              radioButtons(
                                inputId = "powerlevel2",
                                label = "Desired power level:",
                                choices = c("0.80" = 0.8, "0.90" = 0.9, "0.95" = 0.95),
                                selected = 0.8,
                                inline = TRUE
                              )
                            ),
                            mainPanel(
                              plotOutput("esPowerPlot")
                            )
                          )
                 ),
                 
                 tabPanel("Manual Input",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("parameter3", "Parameter:", choices = c("Threshold", "Slope")),
                              numericInput("manual_subjects", "Number of subjects:", value = 30, min = 15, max = 120, step = 1),
                              numericInput("manual_trials", "Number of trials:", value = 40, min = 30, max = 90, step = 1),
                              numericInput("manual_es", "Effect size (d):", value = 0.5, min = 0, max = 1, step = 0.01),
                              actionButton("manual_submit", "Submit")
                            ),
                            mainPanel(uiOutput("my_ft")
                            )
                          )
                 )
)

server <- function(input, output, session) {
  
  # --- Page 1 Plot ---
  output$gridPlot <- renderPlot({
    df <- read_csv(here::here("results","power analysis","Extracted","all_draws.csv")) %>%
      filter(parameter == input$parameter, es == input$es)
    
    power_target <- as.numeric(input$powerlevel)
    
    df %>%
      mutate(test_type = recode(test_type,
                                "h" = "Hierarchical",
                                "s" = "Simple t-test",
                                "u" = "Uncertainty prop. t-test")) %>%
      ggplot() +
      geom_contour(aes(x = nt, y = ss, z = p, colour = test_type),
                   breaks = power_target, linewidth = 1, linetype = 'dashed') +
      theme_minimal() +
      labs(x = "N trials", y = "N participants", fill = "Power",colour='Test type') +
      coord_cartesian(xlim = c(input$trials[1], input$trials[2]),
                      ylim = c(input$subjects[1], input$subjects[2])) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_color_manual(labels=c('Hierarchical','Simple t-test','Uncertainty prop. t-test'),values =c('#D55E00','#CC79A7','#0072B2'))+
      theme(text = element_text(size = 16))
  })
  
  # --- Page 2 Plot: Effect Size vs. Power ---
  output$esPowerPlot <- renderPlot({
    
    if(input$parameter2 == "Threshold"){
      
      hier_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Hierarchical","thresholdmodel.rds"))
      psi_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi","thresholdmodel.rds"))
      psi_un_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi_un","thresholdmodel.rds"))
      
            
      df = rbind(as_draws_df(hier_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "h"),
                 
                 as_draws_df(psi_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "s"),
                 
                 as_draws_df(psi_un_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "u"))
      
      
      
      
      
    }else{
      
      hier_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Hierarchical","slopemodel.rds"))
      psi_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi","slopemodel.rds"))
      psi_un_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi_un","slopemodel.rds"))
      
      
      df = rbind(as_draws_df(hier_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "h"),
                 
                 as_draws_df(psi_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "s"),
                 
                 as_draws_df(psi_un_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$subjects2) %>% 
                   mutate(trials = input$trials2) %>% 
                   mutate(effect_size = list(seq(0,1,by = 0.05))) %>% 
                   unnest() %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975),
                             
                             ci90_lower = quantile(ps, 0.05),
                             ci90_upper = quantile(ps, 0.95),
                             
                             ci80_lower = quantile(ps, 0.10),
                             ci80_upper = quantile(ps, 0.90)) %>% 
                   mutate(test_type = "u"))
      
      
      
      
    }
    
    desired_power=input$powerlevel2
    
    df %>%
      mutate(test_type = recode(test_type,
                                "h" = "Hierarchical",
                                "s" = "Simple t-test",
                                "u" = "Uncertainty prop. t-test")) %>%
      ggplot(aes(x = effect_size, y = mean, col = test_type, fill = test_type)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = ci95_lower,ymax = ci95_upper), alpha = 0.1)+
      geom_ribbon(aes(ymin = ci90_lower,ymax = ci90_upper), alpha = 0.2)+
      geom_ribbon(aes(ymin = ci80_lower,ymax = ci80_upper), alpha = 0.3)+
      theme_minimal() +
      labs(x = "Effect size (d)", y = "Power",colour='Test type',fill='Test type') +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_color_manual(labels=c('Hierarchical','Simple t-test','Uncertainty prop. t-test'),values =c('#D55E00','#CC79A7','#0072B2'))+
      scale_fill_manual(labels=c('Hierarchical','Simple t-test','Uncertainty prop. t-test'),values =c('#D55E00','#CC79A7','#0072B2'))+
      theme(text = element_text(size = 16))+
      geom_hline(yintercept = 0.05, linetype = 2)+
      geom_hline(yintercept = as.numeric(desired_power), linetype = 2)
  })
  
  
  
  output$my_ft <- renderUI({
    
    
    
    if(input$parameter3 == "Threshold"){
      
      hier_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Hierarchical","thresholdmodel.rds"))
      psi_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi","thresholdmodel.rds"))
      psi_un_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi_un","thresholdmodel.rds"))
      
      df = rbind(as_draws_df(hier_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "h"),
                 
                 as_draws_df(psi_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "s"),
                 
                 as_draws_df(psi_un_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "u"))
      
      
      
      
      
    }else{
      hier_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Hierarchical","slopemodel.rds"))
      psi_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi","slopemodel.rds"))
      psi_un_thresholdmodel <- readRDS(here::here("results","power analysis","Extracted","Psi_un","slopemodel.rds"))
      
      
      
      df = rbind(as_draws_df(hier_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "h"),
                 
                 as_draws_df(psi_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "s"),
                 
                 as_draws_df(psi_un_thresholdmodel$draws(c("intercept_alpha","intercept_beta","expo_alpha","expo_beta"))) %>% 
                   select(-contains(".")) %>% 
                   mutate(intercept_alpha = exp(intercept_alpha),
                          intercept_beta = exp(intercept_beta)) %>% 
                   rename_with(~c("intercept_alpha","intercept_beta","expo_alpha1","expo_alpha2","expo_beta1","expo_beta2")) %>% 
                   mutate(subjects = input$manual_subjects) %>% 
                   mutate(trials = input$manual_trials) %>% 
                   mutate(effect_size = input$manual_es) %>% 
                   mutate(ps = get_p(.)) %>% 
                   group_by(subjects,trials, effect_size) %>% 
                   summarize(mean = mean(ps),
                             ci95_lower = quantile(ps, 0.025),
                             ci95_upper = quantile(ps, 0.975)) %>% 
                   mutate(test_type = "u"))
      
      
      
      
    }
    
    
    
    ft <- flextable::flextable(df %>% ungroup() %>% 
                                 select(mean,ci95_lower,ci95_upper,test_type) %>% 
                                 rename_with(~c("mean_power","q5_power","q95_power","Model")) %>% 
                                 mutate(Model = recode(Model,
                                                       "h" = "Hierarchical",
                                                       "s" = "Simple t-test",
                                                       "u" = "Uncertainty prop. t-test"),
                                        mean_power = round(mean_power, 3),
                                        q5_power   = round(q5_power, 3),
                                        q95_power  = round(q95_power, 3))) %>%  
      flextable::autofit()%>%
      htmltools_value()
    
  })
  
  
  
}

shinyApp(ui, server)
