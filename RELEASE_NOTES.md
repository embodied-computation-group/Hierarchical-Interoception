# Release Notes

## v0.9.0-beta (2025-05-09)

### Pre-publication Release

⚠️ **Important**: The methods and results are under review. Use with appropriate caution for research applications.

### What's New

This is the initial release of the Hierarchical-Interoception toolkit, providing comprehensive tools for hierarchical Bayesian modeling of interoceptive psychophysics data.

### Features

- **Hierarchical Psychometric Models**: Complete Stan implementations for HRDT and RRST data
- **Parameter Recovery Validation**: Extensive validation of model parameter recovery
- **Power Analysis Tools**: Interactive Shiny app for power analysis exploration
- **Educational Resources**: Complete BRMS demo with step-by-step workflow

### Components

- **Stan Models**: Population fitting, power analysis, and parameter recovery models
- **Analysis Scripts**: Complete analysis pipeline from data preparation to visualization
- **Shiny App**: Interactive power analysis explorer with three main panels
- **BRMS Demo**: Educational R Markdown with complete workflow example
- **Raw Data**: HRDT and RRST datasets for analysis and validation

### Getting Started

1. Clone the repository
2. Run `source(here::here("setup.R"))` to install dependencies (see below)
3. Follow the BRMS demo for basic usage
4. Use the Shiny app for power analysis exploration

### Dependencies

- R (>= 4.0.0)
- Stan/CmdStan via cmdstanr
- brms, tidyverse, posterior, bayesplot, tidybayes
- shiny, flextable, here, loo, pracma, furrr

### Known Issues


### Citation

If you use this software in your research, please cite:

```
Courtin, A.S., Fischer Ehmsen, J., Banellis, L., Fardo, F., & Allen, M. (2025). 
Hierarchical Bayesian Modelling of Interoceptive Psychophysics. 
bioRxiv. https://doi.org/10.1101/2025.08.27.672360
```

### Roadmap

- **v1.0.0**: Planned after peer review acceptance
- **v0.9.1+**: Bug fixes and minor improvements before v1.0.0

### Support

For questions or issues, please:
1. Check the README.md for usage instructions
2. Open an issue on the GitHub repository
3. Contact the authors directly

---

*This software is released under the MIT License. See LICENSE file for details.*
