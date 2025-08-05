# Simple Example: Hierarchical BRMS for Interoception Data
# Based on the BRMS demo but simplified for 2 groups, single condition

library(tidyverse)
library(brms)
library(cmdstanr)

# Utility function
inv_logit <- function(x) {
  y = 1/(1+exp(-x))
  return(y)
}

# Set seed for reproducibility
set.seed(12345)

# ============================================================================
# DATA SIMULATION
# ============================================================================

# Design parameters
n_subj_per_group <- 25  # 25 subjects per group
n_trials_per_intensity <- 20  # 20 trials per intensity level
intensities <- seq(-30, 30, 5)  # Stimulus intensities from -30 to +30
n_groups <- 2

# Group-level parameters (different for each group)
group_params <- list(
  group1 = list(
    mu_alpha = -5.0,    # Threshold for group 1
    sd_alpha = 8.0,
    mu_beta = -2.0,     # Slope for group 1  
    sd_beta = 0.4,
    mu_lambda = -3.5,   # Lapse rate for group 1
    sd_lambda = 1.5
  ),
  group2 = list(
    mu_alpha = 2.0,     # Threshold for group 2 (higher = worse performance)
    sd_alpha = 8.0,
    mu_beta = -1.8,     # Slope for group 2 (steeper = better sensitivity)
    sd_beta = 0.4,
    mu_lambda = -3.0,   # Lapse rate for group 2
    sd_lambda = 1.5
  )
)

# Simulate subject-level parameters
subjects <- expand.grid(
  group = 1:n_groups,
  subj = 1:n_subj_per_group
) %>%
  mutate(
    # Generate subject-specific parameters
    alpha = case_when(
      group == 1 ~ rnorm(n(), group_params$group1$mu_alpha, group_params$group1$sd_alpha),
      group == 2 ~ rnorm(n(), group_params$group2$mu_alpha, group_params$group2$sd_alpha)
    ),
    beta = case_when(
      group == 1 ~ rnorm(n(), group_params$group1$mu_beta, group_params$group1$sd_beta),
      group == 2 ~ rnorm(n(), group_params$group2$mu_beta, group_params$group2$sd_beta)
    ),
    lambda = case_when(
      group == 1 ~ rnorm(n(), group_params$group1$mu_lambda, group_params$group1$sd_lambda),
      group == 2 ~ rnorm(n(), group_params$group2$mu_lambda, group_params$group2$sd_lambda)
    )
  )

# Expand to create trial-level data
sim_data <- subjects %>%
  # Create all intensity levels for each subject
  crossing(intensity = intensities) %>%
  # Create multiple trials per intensity
  crossing(trial = 1:n_trials_per_intensity) %>%
  # Calculate psychometric function
  mutate(
    # Transform parameters to ensure constraints
    beta_transformed = exp(beta),
    lambda_transformed = inv_logit(lambda) / 2,
    
    # Calculate probability of correct response
    theta = lambda_transformed + 
            (1 - 2 * lambda_transformed) * 
            (0.5 + 0.5 * pnorm(beta_transformed * (intensity - alpha) / sqrt(2))),
    
    # Generate binary responses
    response = rbinom(n(), 1, theta)
  ) %>%
  # Create unique subject IDs
  mutate(subj_id = paste0("group", group, "_subj", subj)) %>%
  # Select relevant columns and ensure group is a factor
  select(subj_id, group, intensity, response, theta) %>%
  mutate(group = factor(group, levels = c(1, 2), labels = c("group1", "group2"))) %>%
  # Group by subject and intensity to get summary statistics
  group_by(subj_id, group, intensity, theta) %>%
  summarise(
    n_trials = n(),
    n_correct = sum(response),
    .groups = 'drop'
  )

# ============================================================================
# DATA VISUALIZATION
# ============================================================================

# Plot individual psychometric functions
ggplot(sim_data, aes(x = intensity, y = n_correct/n_trials, color = group)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = theta, group = subj_id), alpha = 0.3) +
  facet_wrap(~group, labeller = labeller(group = c("group1" = "Group 1", "group2" = "Group 2"))) +
  theme_minimal() +
  labs(
    title = "Simulated Psychometric Functions",
    x = "Stimulus Intensity (ΔBPM)",
    y = "P(correct response)",
    color = "Group"
  ) +
  scale_color_manual(values = c("#0072B2", "#E69F00"))

# ============================================================================
# BRMS MODEL FITTING
# ============================================================================

# Define the BRMS formula for hierarchical psychometric function
# Note: Using Phi() for Stan's cumulative normal distribution (equivalent to R's pnorm)
brm_formula <- bf(
  n_correct | trials(n_trials) ~ 
    (inv_logit(lambda) / 2 + 
     (1 - inv_logit(lambda)) * 
     (0.5 + 0.5 * Phi(exp(beta) * (intensity - alpha) / sqrt(2)))),
  alpha ~ group + (1 | subj_id),  # Random intercepts for subjects
  beta ~ group + (1 | subj_id),   # Random intercepts for subjects
  lambda ~ (1 | subj_id),         # Random intercepts for subjects
  nl = TRUE,
  family = binomial(link = "identity")
)

# Use default priors for simplicity (BRMS will automatically set appropriate priors)
priors <- NULL

# Fit the model
cat("Fitting hierarchical BRMS model...\n")
fit <- brm(
  brm_formula,
  data = sim_data,
  prior = priors,
  chains = 4,
  cores = 4,
  seed = 12345,
  backend = "cmdstan",
  warmup = 1000,
  iter = 2000,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 10
  )
)

# ============================================================================
# MODEL DIAGNOSTICS
# ============================================================================

# Check convergence
cat("\nModel Summary:\n")
print(fit)

# Check R-hat values
cat("\nR-hat values (should be ≤ 1.01):\n")
print(rhat(fit))

# Check effective sample sizes
cat("\nEffective sample sizes (should be ≥ 400):\n")
print(neff_ratio(fit))

# ============================================================================
# RESULTS INTERPRETATION
# ============================================================================

# Extract posterior samples
posterior_samples <- as_draws_df(fit)

# Group-level differences - using parameter names from the fitted model
posterior_samples <- as_draws_df(fit)

# Get the actual parameter names from the model
param_names <- names(posterior_samples)
alpha_group_param <- param_names[grepl("b_alpha_group", param_names) & param_names != "b_alpha_Intercept"][1]
beta_group_param <- param_names[grepl("b_beta_group", param_names) & param_names != "b_beta_Intercept"][1]

# Extract the relevant parameters
group_differences <- posterior_samples %>%
  select(
    b_alpha_Intercept, !!sym(alpha_group_param),
    b_beta_Intercept, !!sym(beta_group_param)
  ) %>%
  mutate(
    # Calculate group means
    alpha_group1 = b_alpha_Intercept,
    alpha_group2 = b_alpha_Intercept + !!sym(alpha_group_param),
    beta_group1 = b_beta_Intercept,
    beta_group2 = b_beta_Intercept + !!sym(beta_group_param),
    
    # Calculate differences
    alpha_diff = !!sym(alpha_group_param),
    beta_diff = !!sym(beta_group_param)
  ) %>%
  summarise(
    # Threshold differences
    alpha_diff_mean = mean(alpha_diff),
    alpha_diff_lower = quantile(alpha_diff, 0.025),
    alpha_diff_upper = quantile(alpha_diff, 0.975),
    alpha_diff_p = 2 * min(mean(alpha_diff > 0), mean(alpha_diff < 0)),
    
    # Slope differences
    beta_diff_mean = mean(beta_diff),
    beta_diff_lower = quantile(beta_diff, 0.025),
    beta_diff_upper = quantile(beta_diff, 0.975),
    beta_diff_p = 2 * min(mean(beta_diff > 0), mean(beta_diff < 0))
  )

cat("\nGroup Differences:\n")
cat("Threshold (α) difference [95% CI]:", 
    round(group_differences$alpha_diff_mean, 2), 
    "[", round(group_differences$alpha_diff_lower, 2), ",", 
    round(group_differences$alpha_diff_upper, 2), "]\n")
cat("Pseudo p-value for threshold difference:", 
    round(group_differences$alpha_diff_p, 3), "\n")

cat("Slope (β) difference [95% CI]:", 
    round(group_differences$beta_diff_mean, 2), 
    "[", round(group_differences$beta_diff_lower, 2), ",", 
    round(group_differences$beta_diff_upper, 2), "]\n")
cat("Pseudo p-value for slope difference:", 
    round(group_differences$beta_diff_p, 3), "\n")

# ============================================================================
# VISUALIZATION OF RESULTS
# ============================================================================

# Plot posterior distributions of group differences
posterior_samples %>%
  select(!!sym(alpha_group_param), !!sym(beta_group_param)) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = case_when(
    parameter == alpha_group_param ~ "Threshold Difference (α)",
    parameter == beta_group_param ~ "Slope Difference (β)"
  )) %>%
  ggplot(aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~parameter, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Posterior Distributions of Group Differences",
    x = "Difference (Group 2 - Group 1)",
    y = "Density"
  ) +
  theme(legend.position = "none")

# ============================================================================
# VISUALIZATION OF FITTED FUNCTIONS
# ============================================================================

# Create a grid of intensity values for plotting
intensity_grid <- seq(-30, 30, 0.5)

# Extract posterior samples for group-level parameters
posterior_summary <- posterior_samples %>%
  select(
    b_alpha_Intercept, !!sym(alpha_group_param),
    b_beta_Intercept, !!sym(beta_group_param),
    b_lambda_Intercept
  ) %>%
  mutate(
    # Calculate group means
    alpha_group1 = b_alpha_Intercept,
    alpha_group2 = b_alpha_Intercept + !!sym(alpha_group_param),
    beta_group1 = b_beta_Intercept,
    beta_group2 = b_beta_Intercept + !!sym(beta_group_param),
    lambda_group1 = b_lambda_Intercept,
    lambda_group2 = b_lambda_Intercept  # No group effect on lambda in this model
  )

# Generate fitted functions for each group
fitted_functions <- expand.grid(
  intensity = intensity_grid,
  group = c("group1", "group2"),
  sample = 1:100  # Use first 100 posterior samples for efficiency
) %>%
  left_join(
    posterior_summary %>% 
      slice(1:100) %>%
      mutate(sample = 1:100),
    by = "sample"
  ) %>%
  mutate(
    # Get the appropriate parameters for each group
    alpha = case_when(
      group == "group1" ~ alpha_group1,
      group == "group2" ~ alpha_group2
    ),
    beta = case_when(
      group == "group1" ~ beta_group1,
      group == "group2" ~ beta_group2
    ),
    lambda = case_when(
      group == "group1" ~ lambda_group1,
      group == "group2" ~ lambda_group2
    ),
    
    # Calculate fitted probability
    fitted_prob = inv_logit(lambda) / 2 + 
                  (1 - inv_logit(lambda)) * 
                  (0.5 + 0.5 * pnorm(exp(beta) * (intensity - alpha) / sqrt(2)))
  )

# Calculate mean and credible intervals for each group
fitted_summary <- fitted_functions %>%
  group_by(intensity, group) %>%
  summarise(
    mean_prob = mean(fitted_prob),
    lower_ci = quantile(fitted_prob, 0.025),
    upper_ci = quantile(fitted_prob, 0.975),
    .groups = 'drop'
  ) %>%
  mutate(group_label = case_when(
    group == "group1" ~ "Group 1",
    group == "group2" ~ "Group 2"
  ))

# Plot fitted functions with credible intervals
ggplot() +
  # Add credible intervals
  geom_ribbon(data = fitted_summary, 
              aes(x = intensity, ymin = lower_ci, ymax = upper_ci, fill = group_label), 
              alpha = 0.3) +
  # Add mean fitted functions
  geom_line(data = fitted_summary, 
            aes(x = intensity, y = mean_prob, color = group_label), 
            linewidth = 1.5) +
  # Add observed data points
  geom_point(data = sim_data, 
             aes(x = intensity, y = n_correct/n_trials, color = factor(group, 
                                                                       levels = c("group1", "group2"),
                                                                       labels = c("Group 1", "Group 2"))), 
             alpha = 0.6, size = 2) +
  # Add individual subject lines (faded)
  geom_line(data = sim_data, 
            aes(x = intensity, y = theta, group = subj_id, 
                color = factor(group, 
                              levels = c("group1", "group2"),
                              labels = c("Group 1", "Group 2"))), 
            alpha = 0.2) +
  # Customize appearance
  scale_color_manual(values = c("#0072B2", "#E69F00")) +
  scale_fill_manual(values = c("#0072B2", "#E69F00")) +
  theme_minimal() +
  labs(
    title = "Fitted Psychometric Functions by Group",
    subtitle = "Solid lines = fitted means, Shaded areas = 95% credible intervals, Points = observed data",
    x = "Stimulus Intensity (ΔBPM)",
    y = "P(correct response)",
    color = "Group",
    fill = "Group"
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  ) +
  coord_cartesian(xlim = c(-30, 30), ylim = c(0, 1))

# Print summary of fitted parameters
cat("\nFitted Group-Level Parameters:\n")
cat("Group 1 - Threshold (α):", round(mean(posterior_summary$alpha_group1), 2), 
    "[", round(quantile(posterior_summary$alpha_group1, 0.025), 2), ",", 
    round(quantile(posterior_summary$alpha_group1, 0.975), 2), "]\n")
cat("Group 2 - Threshold (α):", round(mean(posterior_summary$alpha_group2), 2), 
    "[", round(quantile(posterior_summary$alpha_group2, 0.025), 2), ",", 
    round(quantile(posterior_summary$alpha_group2, 0.975), 2), "]\n")
cat("Group 1 - Slope (β):", round(mean(exp(posterior_summary$beta_group1)), 3), 
    "[", round(quantile(exp(posterior_summary$beta_group1), 0.025), 3), ",", 
    round(quantile(exp(posterior_summary$beta_group1), 0.975), 3), "]\n")
cat("Group 2 - Slope (β):", round(mean(exp(posterior_summary$beta_group2)), 3), 
    "[", round(quantile(exp(posterior_summary$beta_group2), 0.025), 3), ",", 
    round(quantile(exp(posterior_summary$beta_group2), 0.975), 3), "]\n")

cat("\nAnalysis complete! Check the plots for visual results.\n") 