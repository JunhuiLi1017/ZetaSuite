# Install Dependencies for ZetaSuite Shiny App
# This script installs all required packages for the Shiny application

# Function to check if a package is installed
is_package_installed <- function(package_name) {
  return(package_name %in% installed.packages()[,"Package"])
}

# Function to install package if not already installed
install_if_missing <- function(package_name, source = "CRAN") {
  if (!is_package_installed(package_name)) {
    cat("Installing", package_name, "...\n")
    if (source == "CRAN") {
      install.packages(package_name, dependencies = TRUE)
    } else if (source == "BiocManager") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package_name)
    }
    cat("✓", package_name, "installed successfully\n")
  } else {
    cat("✓", package_name, "already installed\n")
  }
}

# Main installation function
install_shiny_dependencies <- function() {
  cat("Installing dependencies for ZetaSuite Shiny App...\n\n")
  
  # Core Shiny packages
  cat("1. Installing core Shiny packages...\n")
  shiny_packages <- c("shiny", "shinydashboard", "DT", "plotly", "shinyjs")
  for (pkg in shiny_packages) {
    install_if_missing(pkg)
  }
  
  # ZetaSuite package dependencies
  cat("\n2. Installing ZetaSuite package dependencies...\n")
  zetasuite_packages <- c(
    "RColorBrewer",
    "Rtsne", 
    "e1071",
    "ggplot2",
    "reshape2",
    "gridExtra",
    "mixtools"
  )
  for (pkg in zetasuite_packages) {
    install_if_missing(pkg)
  }
  
  # Additional useful packages
  cat("\n3. Installing additional useful packages...\n")
  additional_packages <- c(
    "dplyr",
    "tidyr",
    "readr",
    "writexl"
  )
  for (pkg in additional_packages) {
    install_if_missing(pkg)
  }
  
  cat("\n✓ All dependencies installed successfully!\n")
  cat("\nYou can now run the ZetaSuite Shiny app using:\n")
  cat("  shiny::runApp('app.R')\n")
}

# Function to check if all dependencies are available
check_dependencies <- function() {
  cat("Checking ZetaSuite Shiny App dependencies...\n\n")
  
  required_packages <- c(
    "shiny", "shinydashboard", "DT", "plotly", "shinyjs",
    "RColorBrewer", "Rtsne", "e1071", "ggplot2", "reshape2", 
    "gridExtra", "mixtools", "dplyr", "tidyr", "readr", "writexl"
  )
  
  missing_packages <- c()
  installed_packages <- c()
  
  for (pkg in required_packages) {
    if (is_package_installed(pkg)) {
      installed_packages <- c(installed_packages, pkg)
      cat("✓", pkg, "\n")
    } else {
      missing_packages <- c(missing_packages, pkg)
      cat("✗", pkg, "(missing)\n")
    }
  }
  
  cat("\nSummary:\n")
  cat("- Installed packages:", length(installed_packages), "/", length(required_packages), "\n")
  cat("- Missing packages:", length(missing_packages), "/", length(required_packages), "\n")
  
  if (length(missing_packages) > 0) {
    cat("\nMissing packages:\n")
    for (pkg in missing_packages) {
      cat("  -", pkg, "\n")
    }
    cat("\nRun install_shiny_dependencies() to install missing packages.\n")
    return(FALSE)
  } else {
    cat("\n✓ All dependencies are installed and ready!\n")
    return(TRUE)
  }
}

# Function to load all required packages
load_dependencies <- function() {
  cat("Loading ZetaSuite Shiny App dependencies...\n\n")
  
  required_packages <- c(
    "shiny", "shinydashboard", "DT", "plotly", "shinyjs",
    "RColorBrewer", "Rtsne", "e1071", "ggplot2", "reshape2", 
    "gridExtra", "mixtools", "dplyr", "tidyr", "readr", "writexl"
  )
  
  failed_packages <- c()
  loaded_packages <- c()
  
  for (pkg in required_packages) {
    tryCatch({
      library(pkg, character.only = TRUE)
      loaded_packages <- c(loaded_packages, pkg)
      cat("✓", pkg, "loaded\n")
    }, error = function(e) {
      failed_packages <- c(failed_packages, pkg)
      cat("✗", pkg, "failed to load\n")
    })
  }
  
  cat("\nSummary:\n")
  cat("- Successfully loaded:", length(loaded_packages), "/", length(required_packages), "\n")
  cat("- Failed to load:", length(failed_packages), "/", length(required_packages), "\n")
  
  if (length(failed_packages) > 0) {
    cat("\nFailed packages:\n")
    for (pkg in failed_packages) {
      cat("  -", pkg, "\n")
    }
    cat("\nSome packages may need to be reinstalled.\n")
    return(FALSE)
  } else {
    cat("\n✓ All dependencies loaded successfully!\n")
    return(TRUE)
  }
}

# Function to run the Shiny app
run_shiny_app <- function() {
  cat("Starting ZetaSuite Shiny App...\n")
  
  # Check if app.R exists
  if (!file.exists("app.R")) {
    cat("Error: app.R file not found in current directory.\n")
    cat("Please make sure you're in the correct directory.\n")
    return(FALSE)
  }
  
  # Check dependencies
  if (!check_dependencies()) {
    cat("Please install missing dependencies first.\n")
    return(FALSE)
  }
  
  # Load dependencies
  if (!load_dependencies()) {
    cat("Failed to load some dependencies.\n")
    return(FALSE)
  }
  
  # Run the app
  cat("\nLaunching Shiny app...\n")
  cat("The app will open in your default web browser.\n")
  cat("To stop the app, press Ctrl+C in the R console.\n\n")
  
  tryCatch({
    shiny::runApp("app.R", launch.browser = TRUE)
  }, error = function(e) {
    cat("Error launching Shiny app:", e$message, "\n")
    return(FALSE)
  })
}

# Main execution
if (interactive()) {
  cat("ZetaSuite Shiny App Setup\n")
  cat("========================\n\n")
  cat("Available functions:\n")
  cat("1. install_shiny_dependencies() - Install all required packages\n")
  cat("2. check_dependencies() - Check if all packages are installed\n")
  cat("3. load_dependencies() - Load all required packages\n")
  cat("4. run_shiny_app() - Start the Shiny application\n")
  cat("5. generate_example_data() - Generate test data (from generate_example_data.R)\n")
  cat("\nTo get started, run: install_shiny_dependencies()\n")
} else {
  cat("To use this script interactively, run it in R console.\n")
} 