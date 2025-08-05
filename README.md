# Hierarchical-Interoception

This repository contains the complete analysis pipeline for hierarchical psychometric function modeling applied to interoceptive tasks. The project implements Bayesian hierarchical models for analyzing Heart Rate Discrimination Task (HRDT) and Respiratory Resistance Sensitivity Task (RRST) data, with comprehensive power analysis capabilities and interactive visualization tools.

## Table of Contents

1. [Technical Details of the Hierarchical Model](#technical-details-of-the-hierarchical-model)
2. [Data Simulation Workflow](#data-simulation-workflow)
3. [Reproducing Key Figures](#reproducing-key-figures)
4. [Shiny App Deployment and Usage](#shiny-app-deployment-and-usage)

## Technical Details of the Hierarchical Model

### Model Formulation

The hierarchical model implements a psychometric function with three key parameters:

1. **Threshold (α)**: The stimulus intensity at which the probability of correct response is 0.5
2. **Slope (β)**: The steepness of the psychometric function, indicating sensitivity
3. **Lapse rate (λ)**: The asymptotic error rate for extreme stimulus intensities

### Mathematical Framework

The psychometric function is defined as:

$$P(response) = \lambda + (1 - 2\lambda) \cdot \left(0.5 + 0.5 \cdot \text{erf}\left(\frac{\beta \cdot (x - \alpha)}{\sqrt{2}}\right)\right)$$

Where:
- $x$ is the stimulus intensity
- $\alpha$ is the threshold parameter
- $\beta$ is the slope parameter (constrained to be positive)
- $\lambda$ is the lapse rate (constrained to [0, 0.5])
- $\text{erf}$ is the error function

### Hierarchical Structure

The model implements a three-level hierarchy:

1. **Group-level parameters**: Population means and variances for each parameter
2. **Subject-level parameters**: Individual deviations from group means
3. **Trial-level observations**: Binary responses to stimulus presentations

### Stan Implementation

The model is implemented in Stan with the following key components:

```stan
// Group-level parameters
vector[N] gm;  // Group means for all parameters
vector<lower = 0>[N] tau_u;  // Between-participant scales
matrix[N, S] z_expo;  // Participant deviations

// Transformed parameters
vector[S] alpha_int = (gm[1] + (tau_u[1] * z_expo[1,]))';
vector[S] beta_int = (gm[2] + (tau_u[2] * z_expo[2,]))';
vector[S] lapse = (inv_logit(gm[3] + (tau_u[3] * z_expo[3,])) / 2)';
```

### Prior Specifications

The model uses empirically informed priors derived from previous studies:

- **Threshold (α)**: Normal(-8.67, 0.52×10) for intercept, Normal(0, 11.23) for treatment effect
- **Slope (β)**: Normal(-2.3, 0.02×10) for intercept, Normal(0, 0.34) for treatment effect  
- **Lapse rate (λ)**: Normal(-4.32, 0.29)

### Model Comparison

The repository supports both Gaussian and Gumbel cumulative distribution functions:

- **Gaussian**: Uses the error function (erf) for the cumulative normal
- **Gumbel**: Uses the exponential function for the cumulative Gumbel

Model comparison is performed using leave-one-out cross-validation (LOO-CV) with moment matching.

## Data Simulation Workflow

### Using the BRMS Demo

The `apps & wrappers/BRMS demo.Rmd` provides a complete workflow for data simulation and analysis:

1. **Data Simulation**: Generate synthetic data with known parameters
2. **Model Specification**: Define the hierarchical psychometric function
3. **Model Fitting**: Fit using BRMS with Stan backend
4. **Diagnostics**: Check convergence and sampling quality
5. **Inference**: Extract posterior distributions and credible intervals

### Simulation Parameters

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

### Model Specification in BRMS

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

### Running the Analysis

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

## Reproducing Key Figures

### Figure 1: Psychometric Function Forms

**File**: `figures_scripts/Figure1.Rmd`

This figure compares different psychometric function formulations:

```r
PFs <- tibble(x = seq(-6, 11, .001)) %>% 
  mutate(
    Gaussian = pnorm(x - 2),
    Gumbel = 1 - exp(-10^(0.4 * (x - 2.4)))
  ) %>% 
  pivot_longer(cols = !x) %>% 
  rename(Formulation = name)
```

**To reproduce**:
```r
source("figures_scripts/Figure1.Rmd")
```

### Figure 3: Model Recovery

**File**: `figures_scripts/Figure3.Rmd`

This figure shows model recovery performance across different generative models:

```r
# Load model recovery results
load(here::here("results", "model_and_parameter_recovery.RData"))
model_recovery <- results

# Calculate recovery statistics
mod_comp <- tibble(
  winning_model = rep(NaN, 400),
  generative_model = rep(NaN, 400),
  any_pareto_k = rep(NaN, 400),
  significant = rep(NaN, 400),
  meanrhat = rep(NaN, 400),
  sum_div = rep(NaN, 400)
)
```

**To reproduce**:
```r
source("figures_scripts/Figure3.Rmd")
```

### Additional Figures

- **Figure 4**: Parameter recovery plots (`figures_scripts/Figure4.Rmd`)
- **Figure 5**: Power analysis results (`figures_scripts/Figure5.Rmd`) 
- **Figure 6**: Hierarchical vs. simple t-test comparison (`figures_scripts/Figure6.Rmd`)

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

2. **Run the App**:
```r
shiny::runApp("apps & wrappers/shiny app.R")
```

#### Server Deployment

1. **Prepare for Deployment**:
```r
# Ensure all data files are accessible
# Check that results/power analysis/Extracted/all_draws.csv exists
```

2. **Deploy to Shiny Server**:
```bash
# Copy files to server
scp -r "apps & wrappers/" user@server:/path/to/shiny/apps/
scp -r results/ user@server:/path/to/shiny/data/
```

### App Usage

#### Power Grid Panel

1. **Select Parameter**: Choose between "Threshold" or "Slope"
2. **Set Effect Size**: Specify the desired effect size (0.0 to 1.0)
3. **Adjust Sample Ranges**: Use sliders to set subject and trial ranges
4. **Choose Power Level**: Select desired power (0.80, 0.90, or 0.95)
5. **Interpret Results**: Contour lines show combinations achieving target power

#### Effect Size vs. Power Panel

1. **Set Sample Size**: Specify number of subjects and trials
2. **View Power Curve**: See how power changes with effect size
3. **Identify Minimum Effect**: Find the smallest detectable effect size

#### Manual Input Panel

1. **Enter Parameters**: Specify exact sample sizes and effect size
2. **Submit Query**: Get precise power estimates
3. **View Results**: See power for hierarchical, simple t-test, and uncertainty propagation methods

### Data Requirements

The app requires the following data files:
- `results/power analysis/Extracted/all_draws.csv`: Power analysis results
- `results/power analysis/Extracted/all_data.csv`: Summary statistics

### Customization

To modify the app for different analyses:

1. **Update Data Source**: Modify the `all_data` loading in the server function
2. **Add New Parameters**: Extend the parameter selection options
3. **Customize Plots**: Modify the ggplot aesthetics and themes
4. **Add New Panels**: Create additional tabPanels for new functionality

## File Structure

```
Hierarchical-Interoception/
├── Analysis/                    # Main analysis scripts
├── apps & wrappers/            # BRMS demo and Shiny app
├── Data/                       # Raw data files
├── figures/                    # Generated figures
├── figures_scripts/            # Figure generation scripts
├── results/                    # Analysis results and power analysis
├── scripts/                    # Utility functions and main scripts
├── simulated_data/             # Simulated datasets
└── stanmodels/                # Stan model definitions
```

## Dependencies

### Required R Packages
```r
install.packages(c(
  "cmdstanr", "tidyverse", "posterior", "bayesplot", 
  "tidybayes", "rstan", "brms", "pracma", "shiny",
  "flextable", "loo", "furrr"
))
```

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
Courtin, A.S., & Fischer Ehmsen, J. (2025). Hierarchical Bayesian modeling of interoceptive psychometric functions. [Paper DOI to be added]
```

## Contact

For questions or issues with this repository, please contact the authors or open an issue on the repository page.