# Hierarchical-Interoception Setup Script
# This script helps set up the environment for using this research software

cat("Setting up Hierarchical-Interoception v0.9.0-beta...\n")

# Check if renv is available, if not install it
if (!require("renv", quietly = TRUE)) {
  cat("Installing renv for dependency management...\n")
  install.packages("renv")
}

# Initialize renv if not already done
if (!file.exists("renv.lock")) {
  cat("Initializing renv...\n")
  renv::init()
}

# Restore dependencies
cat("Installing/updating dependencies...\n")
renv::restore()

# Check for cmdstanr and Stan installation
if (require("cmdstanr", quietly = TRUE)) {
  cat("Checking Stan installation...\n")
  if (!cmdstanr::cmdstan_version()) {
    cat("Installing CmdStan...\n")
    cmdstanr::install_cmdstan()
  }
} else {
  cat("Warning: cmdstanr not available. Please install manually.\n")
}

cat("Setup complete! You can now use the Hierarchical-Interoception tools.\n")
cat("Start with the BRMS demo: app & demo/BRMS demo.Rmd\n")
cat("Or run the Shiny app: shiny::runApp('app & demo/shiny app.R')\n")
