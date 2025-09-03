# Hierarchical-Interoception

This repository contains the complete analysis pipeline for hierarchical psychometric function modeling applied to intero & exteroceptive tasks. 
The project implements Bayesian hierarchical models for analyzing Heart Rate Discrimination Task (HRDT) and Respiratory Resistance Sensitivity Task (RRST) data, 
with comprehensive power analysis and interactive visualization tools.

## Table of Contents

1. [Technical Details of the Hierarchical Model](#technical-details-of-the-hierarchical-model)
2. [Using the BRMS Demo](#using-the-brms-demo)
3. [Shiny App Deployment and Usage](#shiny-app-deployment-and-usage)

## Technical Details of the Hierarchical Model

### Model Formulation

The hierarchical model implements a psychometric function with three key parameters:

1. **Threshold (α)**: The stimulus intensity at which the probability of correct response is 0.5
2. **Slope (β)**: The steepness of the psychometric function, indicating sensitivity
3. **Lapse rate (λ)**:  The upper (1-$\lambda$) and lower asymptotes. 

### Mathematical Framework

The psychometric function is defined as:

$$P(response) = \lambda + (1 - 2\lambda) \cdot \left(0.5 + 0.5 \cdot \text{erf}\left(\frac{\beta \cdot (x - \alpha)}{\sqrt{2}}\right)\right)$$

Where:
- $x$ is the stimulus intensity
- $\alpha$ is the threshold
- $\beta$ is the slope (constrained to be positive)
- $\lambda$ is the lapse rate (constrained to [0, 0.5])
- $\text{erf}$ is the error function

### Hierarchical Structure

The model implements a three-level hierarchy:

1. **Group-level parameters**: Population means and variances for each parameter
2. **Subject-level parameters**: Individual deviations from group means
3. **Trial-level observations**: Binary responses to stimulus presentations

### Stan Implementation

The model is implemented in Stan with the following key components:

```
// Group-level parameters
vector[N] gm;  // Group means for all parameters
vector<lower = 0>[N] tau_u;  // Between-participant varability
matrix[N, S] z_expo;  // Participant deviations

// Individual psychometric parameters
vector[S] alpha_int = (gm[1] + (tau_u[1] * z_expo[1,]))';
vector[S] beta_int = (gm[2] + (tau_u[2] * z_expo[2,]))';
vector[S] lapse = (inv_logit(gm[3] + (tau_u[3] * z_expo[3,])) / 2)';
```

### Prior Specifications

We recommend using our empirically informed priors derived from the ~500 subject study these can be added in the stancode as follows::
```
gm[1] ~ normal(-8.67, 5.2);  //group mean threshold
gm[2] ~ normal(-2.3,0.2);    //group mean slope (log)
gm[3] ~ normal(-4.32, 0.29); //group mean lapse rate (logit)

```


## Using the BRMS Demo

The `apps & wrappers/BRMS demo.Rmd` provides a complete workflow for data simulation and analysis:

1. **Data Simulation**: Generate synthetic data with known parameters
2. **Model Specification**: Define the hierarchical psychometric function
3. **Model Fitting**: Fit using BRMS with Stan backend
4. **Model Diagnostics**: Check convergence and sampling quality
5. **Model interpretation**: Extract posterior distributions and credible intervals

Below we show the general workflow of the demo. We first simulate data from the following parameters

```r
# Design parameters
n_subj_HRDT <- 50
n_rep_HRDT <- 30
x_seq_HRDT <- rep(seq(-40, 40, 5), n_rep_HRDT)

# Group-level parameters
mu_alpha_HRDT <- -8.6
sd_alpha_HRDT <- 11.6
mu_beta_HRDT <- -2.3
sd_beta_HRDT <- 0.3
mu_lambda_HRDT <- -4.2
sd_lambda_HRDT <- 2
```

Then specifying the model in brms:


```r
bff_HRDT = bf(
  y|trials(n) ~ (inv_logit(lambda) / 2 + 
                 (1-inv_logit(lambda)) * 
                 (0.5 + 0.5 * erf(exp(beta)*(x-alpha)/(sqrt(2))))),
  alpha ~ condition + (condition | subj),  
  beta ~ condition + (condition | subj),  
  lambda ~ (1 | subj),
  nl = TRUE,
  family = binomial(link = "identity")
)
```

Running this model on the simulated data:

```r
fit_eip_HRDT <- brm(
  bff_HRDT,
  data = sim_data_HRDT,
  chains = 4,
  cores = 4,
  seed = seed,
  backend = "cmdstan",
  warmup = 2000,
  iter = 4000,           
  control = list(
    adapt_delta = 0.9,
    max_treedepth = 10
  ),
  prior = empirically_informed_priors_HRDT
)
```


## Shiny App Deployment and Usage

### App Overview

The Shiny app (`apps & wrappers/shiny app.R`) provides an interactive interface for exploring power analysis results with three main panels:

1. **Power Grid**: Visualize power contours across different sample sizes
2. **Effect Size vs. Power**: Examine power as a function of effect size
3. **Manual Input**: Get specific power estimates for custom parameters

### Deployment

#### Local Deployment

1. **Install Dependencies**:
```r
install.packages(c("shiny", "tidyverse", "flextable", "posterior"))
```

2. **Running the App**:
Either enter this in the R-console
```r
shiny::runApp("apps & wrappers/shiny app.R")
```
or go to the "apps & wrappers" directory and run the `shiny app.R` script 



### App Usage

#### Power Grid Panel

1. **Parameter**: Choose between "Threshold" or "Slope"
2. **Effect Size**: Specify the desired effect size (0.0 to 1.0)
3. **Adjust Sample and trial Ranges**: Use sliders to set subject and trial ranges
4. **Choose Power Level**: Select desired power (0.80, 0.90, or 0.95)
5. **Interpret Results**: Contour lines show combinations achieving target power

#### Effect Size vs. Power Panel

1. **Parameter**: Choose between "Threshold" or "Slope"
2. **Set Sample Size**: Specify number of subjects and trials
3. **View Power Curve**: See how power changes with effect size
3. **Interpret Results**: Psychometric functions show (efffect size by power).

#### Manual Input Panel

1. **Enter Parameters**: Specify exact sample sizes and effect size and parameter
2. **Submit Query**: Get precise power estimates for above
3. **View Results**: See power for hierarchical, simple t-test, and uncertainty propagation methods

## File Structure

```
Hierarchical-Interoception/
├── Analysis/                    # Model comparison for the HRDT
├── apps & wrappers/             # BRMS demo and Shiny app
├── Data/                        # Raw data files
├── figures/                     # Generated figures
├── figures_scripts/             # Figure generation scripts
├── results/                     # Analysis results and power analysis
├── scripts/                     # Utility functions and main scripts
├── simulated_data/              # Simulated datasets (model recovery and power analysis)
└── stanmodels/                  # Stan model definitions (VMP-data refit, poweranalysis and model recovery)
```

## Dependencies

### Stan Installation
```r
# Install cmdstanr
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Check cmdstan installation
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

## Citation

If you use this code in your research, please cite:

```
Courtin, A.S., & Fischer Ehmsen, J. (2025). Hierarchical Bayesian modeling of interoceptive psychometric functions. https://doi.org/10.1101/2025.08.27.672360
```

## Contact

For questions or issues with this repository, please contact the authors or open an issue on the repository page.
