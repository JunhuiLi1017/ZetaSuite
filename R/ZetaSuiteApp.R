#' Launch ZetaSuite Shiny Application
#'
#' @description
#' Launches the ZetaSuite Shiny web application for interactive analysis of 
#' high-throughput screening data and single-cell RNA-seq quality control.
#'
#' @details
#' The Shiny app provides a user-friendly interface for:
#' \itemize{
#'   \item Quality Control Analysis
#'   \item Z-score Normalization
#'   \item Event Coverage Analysis
#'   \item Zeta Score Calculation
#'   \item SVM-based Background Correction
#'   \item Screen Strength Analysis
#'   \item Single Cell Quality Control
#'   \item Interactive visualizations and data export
#' }
#'
#' @param launch.browser Logical. Should the app launch in the default browser?
#'   Default is TRUE.
#' @param port Integer. Port number for the Shiny app. Default is NULL (random port).
#' @param host Character. Host address. Default is "127.0.0.1" (localhost).
#'
#' @return Launches the Shiny application in a web browser.
#'
#' @examples
#' \dontrun{
#' # Launch the ZetaSuite Shiny app
#' ZetaSuiteApp()
#' 
#' # Launch without opening browser automatically
#' ZetaSuiteApp(launch.browser = FALSE)
#' 
#' # Launch on a specific port
#' ZetaSuiteApp(port = 3838)
#' }
#'
#' @export
ZetaSuiteApp <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {
  
  # Check if required packages are installed
  required_packages <- c("shiny", "shinydashboard", "DT", "plotly", "shinyjs")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("The following packages are required to run the ZetaSuite Shiny app but are not installed:\n",
         paste(missing_packages, collapse = ", "), "\n\n",
         "Please install them using:\n",
         "install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))")
  }
  
  # Load required packages
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(plotly)
  library(shinyjs)
  
  # Get the path to the Shiny app
  app_path <- system.file("shiny", package = "ZetaSuite")
  
  if (app_path == "") {
    stop("Shiny app files not found. Please ensure the ZetaSuite package is properly installed.")
  }
  
  # Launch the Shiny app
  shiny::runApp(
    appDir = app_path,
    launch.browser = launch.browser,
    port = port,
    host = host
  )
} 