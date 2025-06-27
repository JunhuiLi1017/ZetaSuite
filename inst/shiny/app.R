# ZetaSuite Shiny Application
# A web interface for high-throughput screening data analysis

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(shinyjs)

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "ZetaSuite Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome", tabName = "welcome", icon = icon("home")),
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Example Data", tabName = "example", icon = icon("database")),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle")),
      menuItem("Z-score Analysis", tabName = "zscore", icon = icon("chart-line")),
      menuItem("Event Coverage", tabName = "eventcoverage", icon = icon("layer-group")),
      menuItem("Zeta Score", tabName = "zeta", icon = icon("calculator")),
      menuItem("SVM Analysis", tabName = "svm", icon = icon("brain")),
      menuItem("FDR Cutoff", tabName = "fdr", icon = icon("cut")),
      menuItem("Single Cell QC", tabName = "singlecell", icon = icon("mobile")),
      menuItem("Results", tabName = "results", icon = icon("download")),
      menuItem("Help & Support", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tabItems(
      
      # Welcome Tab
      tabItem(tabName = "welcome",
        fluidRow(
          box(title = "Welcome to ZetaSuite", width = 12,
            h3("Multi-dimensional High-throughput Data Analysis"),
            p("ZetaSuite is an R package designed for analyzing multi-dimensional high-throughput screening data, particularly two-dimensional RNAi screens and single-cell RNA sequencing data."),
            br(),
            h4("Key Features:"),
            tags$ul(
              tags$li("Quality Control Analysis"),
              tags$li("Z-score Normalization"),
              tags$li("Event Coverage Analysis"),
              tags$li("Zeta Score Calculation"),
              tags$li("SVM-based Background Correction"),
              tags$li("Screen Strength Analysis"),
              tags$li("Single Cell Quality Control")
            ),
            br(),
            h4("Getting Started:"),
            p("1. Use the 'Example Data' tab to explore the built-in dataset"),
            p("2. Or upload your own data in the 'Data Upload' tab"),
            p("3. Follow the analysis workflow through the tabs"),
            br(),
            p("For more information, see the package vignette and documentation.")
          )
        )
      ),
      
      # Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(title = "Upload Data Files", width = 12,
            fileInput("countMat", "Upload Count Matrix (CSV)", accept = ".csv"),
            fileInput("negGene", "Upload Negative Control Genes (CSV)", accept = ".csv"),
            fileInput("posGene", "Upload Positive Control Genes (CSV)", accept = ".csv"),
            fileInput("nonExpGene", "Upload Non-expressed Genes (CSV, optional)", accept = ".csv"),
            actionButton("loadData", "Load Data", class = "btn-primary")
          )
        ),
        fluidRow(
          box(title = "Data Preview", width = 12,
            tabsetPanel(
              tabPanel("Count Matrix", DT::dataTableOutput("countMatPreview")),
              tabPanel("Negative Controls", DT::dataTableOutput("negGenePreview")),
              tabPanel("Positive Controls", DT::dataTableOutput("posGenePreview"))
            )
          )
        )
      ),
      
      # Example Data Tab
      tabItem(tabName = "example",
        fluidRow(
          box(title = "Example Dataset", width = 12,
            h4("HTS2 Screening Dataset"),
            p("This example dataset contains:"),
            tags$ul(
              tags$li("1,609 genes × 100 alternative splicing events"),
              tags$li("30 negative control genes (non-specific siRNAs)"),
              tags$li("20 positive control genes (PTB-targeting siRNAs)"),
              tags$li("50 non-expressed genes (RPKM < 1 in HeLa cells)")
            ),
            br(),
            actionButton("loadExampleData", "Load Example Data", class = "btn-success")
          )
        ),
        fluidRow(
          box(title = "Example Data Summary", width = 12,
            verbatimTextOutput("exampleDataSummary")
          )
        )
      ),
      
      # Quality Control Tab
      tabItem(tabName = "qc",
        fluidRow(
          box(title = "Quality Control Analysis", width = 12,
            actionButton("runQC", "Run Quality Control", class = "btn-success"),
            br(), br(),
            textOutput("qcStatus")
          )
        ),
        fluidRow(
          box(title = "Score Distribution", width = 6,
            plotlyOutput("scoreQCPlot")
          ),
          box(title = "t-SNE Plot", width = 6,
            plotlyOutput("tsnePlot")
          )
        ),
        fluidRow(
          box(title = "Box Plots", width = 6,
            plotlyOutput("boxPlot")
          ),
          box(title = "SSMD Distribution", width = 6,
            plotlyOutput("ssmdPlot")
          )
        )
      ),
      
      # Z-score Analysis Tab
      tabItem(tabName = "zscore",
        fluidRow(
          box(title = "Z-score Normalization", width = 12,
            actionButton("runZscore", "Calculate Z-scores", class = "btn-success"),
            br(), br(),
            textOutput("zscoreStatus")
          )
        ),
        fluidRow(
          box(title = "Z-score Matrix Preview", width = 12,
            DT::dataTableOutput("zscorePreview")
          )
        )
      ),
      
      # Event Coverage Tab
      tabItem(tabName = "eventcoverage",
        fluidRow(
          box(title = "Event Coverage Parameters", width = 6,
            numericInput("binNum", "Number of Bins", value = 100, min = 10, max = 500),
            checkboxInput("combine", "Combine Directions", value = TRUE),
            actionButton("runEventCoverage", "Calculate Event Coverage", class = "btn-success")
          ),
          box(title = "Event Coverage Status", width = 6,
            textOutput("eventCoverageStatus")
          )
        ),
        fluidRow(
          box(title = "Decrease Direction", width = 6,
            plotlyOutput("ecDecreasePlot")
          ),
          box(title = "Increase Direction", width = 6,
            plotlyOutput("ecIncreasePlot")
          )
        )
      ),
      
      # Zeta Score Tab
      tabItem(tabName = "zeta",
        fluidRow(
          box(title = "Zeta Score Parameters", width = 6,
            checkboxInput("useSVM", "Use SVM Curves", value = FALSE),
            actionButton("runZeta", "Calculate Zeta Scores", class = "btn-success")
          ),
          box(title = "Zeta Score Status", width = 6,
            textOutput("zetaStatus")
          )
        ),
        fluidRow(
          box(title = "Zeta Scores Preview", width = 12,
            DT::dataTableOutput("zetaPreview")
          )
        ),
        fluidRow(
          box(title = "Top Hits by Zeta_D", width = 6,
            DT::dataTableOutput("topDecreaseTable")
          ),
          box(title = "Top Hits by Zeta_I", width = 6,
            DT::dataTableOutput("topIncreaseTable")
          )
        )
      ),
      
      # SVM Analysis Tab
      tabItem(tabName = "svm",
        fluidRow(
          box(title = "SVM Analysis", width = 12,
            actionButton("runSVM", "Run SVM Analysis", class = "btn-success"),
            br(), br(),
            textOutput("svmStatus")
          )
        ),
        fluidRow(
          box(title = "SVM Results", width = 12,
            tabsetPanel(
              tabPanel("Decrease Direction", DT::dataTableOutput("svmDecreaseTable")),
              tabPanel("Increase Direction", DT::dataTableOutput("svmIncreaseTable"))
            )
          )
        )
      ),
      
      # FDR Cutoff Tab
      tabItem(tabName = "fdr",
        fluidRow(
          box(title = "FDR Cutoff Parameters", width = 6,
            checkboxInput("combineFDR", "Combine Directions", value = FALSE),
            actionButton("runFDR", "Calculate FDR Cutoffs", class = "btn-success")
          ),
          box(title = "FDR Cutoff Status", width = 6,
            textOutput("fdrStatus")
          )
        ),
        fluidRow(
          box(title = "Zeta Score Distribution by Type", width = 6,
            plotlyOutput("zetaTypePlot")
          ),
          box(title = "Screen Strength Curves", width = 6,
            plotlyOutput("ssCutoffPlot")
          )
        ),
        fluidRow(
          box(title = "FDR Cutoff Results", width = 12,
            DT::dataTableOutput("fdrTable")
          )
        ),
        fluidRow(
          box(title = "Hit Selection", width = 12,
            numericInput("ssThreshold", "Screen Strength Threshold", value = 0.8, min = 0, max = 1, step = 0.1),
            actionButton("selectHits", "Select Hits", class = "btn-primary"),
            br(), br(),
            verbatimTextOutput("hitSelectionResults")
          )
        )
      ),
      
      # Single Cell QC Tab
      tabItem(tabName = "singlecell",
        fluidRow(
          box(title = "Single Cell Data Upload", width = 12,
            fileInput("countMatSC", "Upload Single Cell Count Matrix (CSV)", accept = ".csv"),
            actionButton("loadSCData", "Load Single Cell Data", class = "btn-primary")
          )
        ),
        fluidRow(
          box(title = "Single Cell QC Parameters", width = 6,
            numericInput("binNumSC", "Number of Bins", value = 10, min = 5, max = 100),
            checkboxInput("filterSC", "Filter Low Count Cells", value = TRUE),
            actionButton("runSingleCell", "Run Single Cell QC", class = "btn-success")
          ),
          box(title = "Single Cell QC Status", width = 6,
            textOutput("singleCellStatus")
          )
        ),
        fluidRow(
          box(title = "Zeta Score Distribution", width = 12,
            plotlyOutput("singleCellPlot")
          )
        ),
        fluidRow(
          box(title = "Single Cell Results", width = 12,
            DT::dataTableOutput("singleCellTable")
          )
        )
      ),
      
      # Results Tab
      tabItem(tabName = "results",
        fluidRow(
          box(title = "Download Results", width = 12,
            downloadButton("downloadZscore", "Download Z-scores"),
            downloadButton("downloadZeta", "Download Zeta Scores"),
            downloadButton("downloadFDR", "Download FDR Results"),
            downloadButton("downloadHits", "Download Selected Hits"),
            downloadButton("downloadSingleCell", "Download Single Cell Results"),
            downloadButton("downloadReport", "Download Analysis Report")
          )
        ),
        fluidRow(
          box(title = "Analysis Summary", width = 12,
            verbatimTextOutput("analysisSummary")
          )
        )
      ),
      
      # Help & Support Tab
      tabItem(tabName = "help",
        fluidRow(
          box(title = "Documentation & Support", width = 12,
            h4("Package Documentation"),
            p("For detailed documentation and examples, see the package vignette:"),
            code("vignette(\"ZetaSuite\")"),
            br(), br(),
            h4("Bug Reports & Feature Requests"),
            p("If you encounter any bugs or have feature requests, please report them on our GitHub issues page:"),
            tags$a(href = "https://github.com/JunhuiLi1017/ZetaSuite/issues", 
                   "Report a Bug or Request Feature", 
                   target = "_blank", 
                   class = "btn btn-warning"),
            br(), br(),
            h4("Citation"),
            p("If you use ZetaSuite in your research, please cite:"),
            p("Hao, Y., Shao, C., Zhao, G., Fu, X.D. (2021). ZetaSuite: A Computational Method for Analyzing Multi-dimensional High-throughput Data, Reveals Genes with Opposite Roles in Cancer Dependency. Forthcoming"),
            br(),
            h4("Contact"),
            p("For questions about the package, contact the maintainer:"),
            p("Junhui Li <ljh.biostat@gmail.com>")
          )
        ),
        fluidRow(
          box(title = "Troubleshooting", width = 12,
            h4("Common Issues"),
            tags$ul(
              tags$li("Make sure all required packages are installed"),
              tags$li("Check that CSV files have correct column headers and data types"),
              tags$li("For large datasets, consider reducing the number of bins"),
              tags$li("Ensure you have sufficient positive and negative control samples")
            ),
            br(),
            h4("Data Format Requirements"),
            p("Count Matrix: Rows = Genes/siRNAs, Columns = Readouts/conditions"),
            p("Control Files: First column should contain gene/siRNA identifiers"),
            p("All files should be in CSV format with proper headers")
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive values to store data and results
  values <- reactiveValues(
    countMat = NULL,
    negGene = NULL,
    posGene = NULL,
    nonExpGene = NULL,
    countMatSC = NULL,
    zscoreVal = NULL,
    ecData = NULL,
    zetaData = NULL,
    svmData = NULL,
    fdrData = NULL,
    singleCellData = NULL,
    qcResults = NULL,
    selectedHits = NULL
  )
  
  # Load example data
  observeEvent(input$loadExampleData, {
    tryCatch({
      # Load ZetaSuite package data
      library(ZetaSuite)
      data(countMat)
      data(negGene)
      data(posGene)
      data(nonExpGene)
      data(ZseqList)
      data(SVMcurve)
      
      values$countMat <- countMat
      values$negGene <- negGene
      values$posGene <- posGene
      values$nonExpGene <- nonExpGene
      
      showNotification("Example data loaded successfully!", type = "success")
    }, error = function(e) {
      showNotification(paste("Error loading example data:", e$message), type = "error")
    })
  })
  
  # Load main data
  observeEvent(input$loadData, {
    req(input$countMat, input$negGene, input$posGene)
    
    tryCatch({
      values$countMat <- read.csv(input$countMat$datapath, row.names = 1)
      values$negGene <- read.csv(input$negGene$datapath)
      values$posGene <- read.csv(input$posGene$datapath)
      
      if (!is.null(input$nonExpGene)) {
        values$nonExpGene <- read.csv(input$nonExpGene$datapath)
      }
      
      showNotification("Data loaded successfully!", type = "success")
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  
  # Load single cell data
  observeEvent(input$loadSCData, {
    req(input$countMatSC)
    
    tryCatch({
      values$countMatSC <- read.csv(input$countMatSC$datapath, row.names = 1)
      showNotification("Single cell data loaded successfully!", type = "success")
    }, error = function(e) {
      showNotification(paste("Error loading single cell data:", e$message), type = "error")
    })
  })
  
  # Example data summary
  output$exampleDataSummary <- renderPrint({
    if (!is.null(values$countMat)) {
      cat("=== Example Dataset Summary ===\n\n")
      cat("Count matrix dimensions:", dim(values$countMat), "\n")
      cat("Negative controls:", nrow(values$negGene), "genes\n")
      cat("Positive controls:", nrow(values$posGene), "genes\n")
      if (!is.null(values$nonExpGene)) {
        cat("Non-expressed genes:", nrow(values$nonExpGene), "genes\n")
      }
      cat("\nData ready for analysis!\n")
    } else {
      cat("Click 'Load Example Data' to load the built-in dataset.\n")
    }
  })
  
  # Data preview outputs
  output$countMatPreview <- DT::renderDataTable({
    req(values$countMat)
    DT::datatable(head(values$countMat, 10), options = list(scrollX = TRUE))
  })
  
  output$negGenePreview <- DT::renderDataTable({
    req(values$negGene)
    DT::datatable(values$negGene, options = list(scrollX = TRUE))
  })
  
  output$posGenePreview <- DT::renderDataTable({
    req(values$posGene)
    DT::datatable(values$posGene, options = list(scrollX = TRUE))
  })
  
  # Quality Control
  observeEvent(input$runQC, {
    req(values$countMat, values$negGene, values$posGene)
    
    tryCatch({
      withProgress(message = "Running Quality Control...", {
        values$qcResults <- QC(values$countMat, values$negGene, values$posGene)
      })
      output$qcStatus <- renderText("Quality Control completed successfully!")
      showNotification("Quality Control completed!", type = "success")
    }, error = function(e) {
      output$qcStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in Quality Control:", e$message), type = "error")
    })
  })
  
  # Z-score Analysis
  observeEvent(input$runZscore, {
    req(values$countMat, values$negGene)
    
    tryCatch({
      withProgress(message = "Calculating Z-scores...", {
        values$zscoreVal <- Zscore(values$countMat, values$negGene)
      })
      output$zscoreStatus <- renderText("Z-score calculation completed successfully!")
      showNotification("Z-scores calculated!", type = "success")
    }, error = function(e) {
      output$zscoreStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in Z-score calculation:", e$message), type = "error")
    })
  })
  
  # Event Coverage
  observeEvent(input$runEventCoverage, {
    req(values$zscoreVal, values$negGene, values$posGene)
    
    tryCatch({
      withProgress(message = "Calculating Event Coverage...", {
        values$ecData <- EventCoverage(values$zscoreVal, values$negGene, values$posGene, 
                                      input$binNum, input$combine)
      })
      output$eventCoverageStatus <- renderText("Event Coverage calculation completed successfully!")
      showNotification("Event Coverage calculated!", type = "success")
    }, error = function(e) {
      output$eventCoverageStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in Event Coverage calculation:", e$message), type = "error")
    })
  })
  
  # Zeta Score
  observeEvent(input$runZeta, {
    req(values$zscoreVal, values$ecData)
    
    tryCatch({
      withProgress(message = "Calculating Zeta Scores...", {
        if (input$useSVM && !is.null(values$svmData)) {
          values$zetaData <- Zeta(values$zscoreVal, values$ecData[[1]]$ZseqList, 
                                 values$svmData, SVM = TRUE)
        } else {
          values$zetaData <- Zeta(values$zscoreVal, values$ecData[[1]]$ZseqList, SVM = FALSE)
        }
      })
      output$zetaStatus <- renderText("Zeta Score calculation completed successfully!")
      showNotification("Zeta Scores calculated!", type = "success")
    }, error = function(e) {
      output$zetaStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in Zeta Score calculation:", e$message), type = "error")
    })
  })
  
  # SVM Analysis
  observeEvent(input$runSVM, {
    req(values$ecData)
    
    tryCatch({
      withProgress(message = "Running SVM Analysis...", {
        values$svmData <- SVM(values$ecData)
      })
      output$svmStatus <- renderText("SVM Analysis completed successfully!")
      showNotification("SVM Analysis completed!", type = "success")
    }, error = function(e) {
      output$svmStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in SVM Analysis:", e$message), type = "error")
    })
  })
  
  # FDR Cutoff
  observeEvent(input$runFDR, {
    req(values$zetaData, values$negGene, values$posGene, values$nonExpGene)
    
    tryCatch({
      withProgress(message = "Calculating FDR Cutoffs...", {
        values$fdrData <- FDRcutoff(values$zetaData, values$negGene, values$posGene, 
                                   values$nonExpGene, input$combineFDR)
      })
      output$fdrStatus <- renderText("FDR Cutoff calculation completed successfully!")
      showNotification("FDR Cutoffs calculated!", type = "success")
    }, error = function(e) {
      output$fdrStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in FDR Cutoff calculation:", e$message), type = "error")
    })
  })
  
  # Hit Selection
  observeEvent(input$selectHits, {
    req(values$zetaData, values$fdrData)
    
    tryCatch({
      fdr_table <- values$fdrData[[1]]
      selected_threshold <- fdr_table[fdr_table$SS >= input$ssThreshold, ]
      
      if (nrow(selected_threshold) > 0) {
        best_threshold <- selected_threshold[which.max(selected_threshold$SS), ]
        combined_zeta <- values$zetaData$Zeta_D + values$zetaData$Zeta_I
        hits <- names(combined_zeta[combined_zeta >= best_threshold$Cut_Off])
        values$selectedHits <- data.frame(
          Gene = hits,
          Zeta_D = values$zetaData[hits, "Zeta_D"],
          Zeta_I = values$zetaData[hits, "Zeta_I"],
          Combined_Zeta = combined_zeta[hits],
          stringsAsFactors = FALSE
        )
        
        output$hitSelectionResults <- renderPrint({
          cat("=== Hit Selection Results ===\n\n")
          cat("Screen Strength threshold:", input$ssThreshold, "\n")
          cat("Selected threshold:", best_threshold$Cut_Off, "\n")
          cat("Screen Strength:", best_threshold$SS, "\n")
          cat("Total hits identified:", length(hits), "\n")
          cat("Apparent FDR:", best_threshold$aFDR, "\n")
        })
        
        showNotification(paste("Selected", length(hits), "hits!"), type = "success")
      } else {
        output$hitSelectionResults <- renderPrint({
          cat("No thresholds found with Screen Strength >=", input$ssThreshold, "\n")
          cat("Try lowering the threshold or running FDR analysis with different parameters.\n")
        })
      }
    }, error = function(e) {
      output$hitSelectionResults <- renderPrint({
        cat("Error in hit selection:", e$message, "\n")
      })
    })
  })
  
  # Single Cell QC
  observeEvent(input$runSingleCell, {
    req(values$countMatSC)
    
    tryCatch({
      withProgress(message = "Running Single Cell QC...", {
        values$singleCellData <- ZetaSuitSC(values$countMatSC, input$binNumSC, input$filterSC)
      })
      output$singleCellStatus <- renderText("Single Cell QC completed successfully!")
      showNotification("Single Cell QC completed!", type = "success")
    }, error = function(e) {
      output$singleCellStatus <- renderText(paste("Error:", e$message))
      showNotification(paste("Error in Single Cell QC:", e$message), type = "error")
    })
  })
  
  # Plot outputs
  output$scoreQCPlot <- renderPlotly({
    req(values$qcResults)
    ggplotly(values$qcResults$score_qc)
  })
  
  output$tsnePlot <- renderPlotly({
    req(values$qcResults)
    ggplotly(values$qcResults$tSNE_QC)
  })
  
  output$boxPlot <- renderPlotly({
    req(values$qcResults)
    ggplotly(values$qcResults$QC_box)
  })
  
  output$ssmdPlot <- renderPlotly({
    req(values$qcResults)
    ggplotly(values$qcResults$QC_SSMD)
  })
  
  output$ecDecreasePlot <- renderPlotly({
    req(values$ecData)
    ggplotly(values$ecData[[2]]$EC_jitter_D)
  })
  
  output$ecIncreasePlot <- renderPlotly({
    req(values$ecData)
    ggplotly(values$ecData[[2]]$EC_jitter_I)
  })
  
  output$zetaTypePlot <- renderPlotly({
    req(values$fdrData)
    ggplotly(values$fdrData[[2]]$Zeta_type)
  })
  
  output$ssCutoffPlot <- renderPlotly({
    req(values$fdrData)
    ggplotly(values$fdrData[[2]]$SS_cutOff)
  })
  
  output$singleCellPlot <- renderPlotly({
    req(values$singleCellData)
    ggplotly(values$singleCellData[[2]])
  })
  
  # Table outputs
  output$zscorePreview <- DT::renderDataTable({
    req(values$zscoreVal)
    DT::datatable(head(values$zscoreVal, 10), options = list(scrollX = TRUE))
  })
  
  output$zetaPreview <- DT::renderDataTable({
    req(values$zetaData)
    DT::datatable(head(values$zetaData, 10), options = list(scrollX = TRUE))
  })
  
  output$topDecreaseTable <- DT::renderDataTable({
    req(values$zetaData)
    top_decrease <- head(values$zetaData[order(values$zetaData$Zeta_D, decreasing = TRUE), ], 10)
    DT::datatable(top_decrease, options = list(scrollX = TRUE))
  })
  
  output$topIncreaseTable <- DT::renderDataTable({
    req(values$zetaData)
    top_increase <- head(values$zetaData[order(values$zetaData$Zeta_I, decreasing = TRUE), ], 10)
    DT::datatable(top_increase, options = list(scrollX = TRUE))
  })
  
  output$svmDecreaseTable <- DT::renderDataTable({
    req(values$svmData)
    DT::datatable(values$svmData$cutOffD, options = list(scrollX = TRUE))
  })
  
  output$svmIncreaseTable <- DT::renderDataTable({
    req(values$svmData)
    DT::datatable(values$svmData$cutOffI, options = list(scrollX = TRUE))
  })
  
  output$fdrTable <- DT::renderDataTable({
    req(values$fdrData)
    DT::datatable(values$fdrData[[1]], options = list(scrollX = TRUE))
  })
  
  output$singleCellTable <- DT::renderDataTable({
    req(values$singleCellData)
    DT::datatable(values$singleCellData[[1]], options = list(scrollX = TRUE))
  })
  
  # Download handlers
  output$downloadZscore <- downloadHandler(
    filename = function() { "zscore_results.csv" },
    content = function(file) {
      req(values$zscoreVal)
      write.csv(values$zscoreVal, file)
    }
  )
  
  output$downloadZeta <- downloadHandler(
    filename = function() { "zeta_scores.csv" },
    content = function(file) {
      req(values$zetaData)
      write.csv(values$zetaData, file)
    }
  )
  
  output$downloadFDR <- downloadHandler(
    filename = function() { "fdr_results.csv" },
    content = function(file) {
      req(values$fdrData)
      write.csv(values$fdrData[[1]], file)
    }
  )
  
  output$downloadHits <- downloadHandler(
    filename = function() { "selected_hits.csv" },
    content = function(file) {
      req(values$selectedHits)
      write.csv(values$selectedHits, file, row.names = FALSE)
    }
  )
  
  output$downloadSingleCell <- downloadHandler(
    filename = function() { "single_cell_results.csv" },
    content = function(file) {
      req(values$singleCellData)
      write.csv(values$singleCellData[[1]], file)
    }
  )
  
  output$downloadReport <- downloadHandler(
    filename = function() { "zetaSuite_analysis_report.txt" },
    content = function(file) {
      cat("=== ZetaSuite Analysis Report ===\n\n", file = file)
      
      if (!is.null(values$countMat)) {
        cat("Data Summary:\n", file = file, append = TRUE)
        cat("Count matrix dimensions:", dim(values$countMat), "\n", file = file, append = TRUE)
        cat("Negative controls:", nrow(values$negGene), "genes\n", file = file, append = TRUE)
        cat("Positive controls:", nrow(values$posGene), "genes\n", file = file, append = TRUE)
        if (!is.null(values$nonExpGene)) {
          cat("Non-expressed genes:", nrow(values$nonExpGene), "genes\n", file = file, append = TRUE)
        }
        cat("\n", file = file, append = TRUE)
      }
      
      if (!is.null(values$zetaData)) {
        cat("Zeta Score Summary:\n", file = file, append = TRUE)
        cat("Number of genes:", nrow(values$zetaData), "\n", file = file, append = TRUE)
        cat("Zeta_D range:", range(values$zetaData$Zeta_D), "\n", file = file, append = TRUE)
        cat("Zeta_I range:", range(values$zetaData$Zeta_I), "\n", file = file, append = TRUE)
        cat("\n", file = file, append = TRUE)
      }
      
      if (!is.null(values$selectedHits)) {
        cat("Hit Selection Results:\n", file = file, append = TRUE)
        cat("Total hits identified:", nrow(values$selectedHits), "\n", file = file, append = TRUE)
        cat("\n", file = file, append = TRUE)
      }
      
      cat("Analysis completed on:", Sys.time(), "\n", file = file, append = TRUE)
    }
  )
  
  # Analysis summary
  output$analysisSummary <- renderPrint({
    cat("=== ZetaSuite Analysis Summary ===\n\n")
    
    if (!is.null(values$countMat)) {
      cat("✓ Count matrix loaded:", nrow(values$countMat), "genes ×", ncol(values$countMat), "samples\n")
    }
    
    if (!is.null(values$zscoreVal)) {
      cat("✓ Z-scores calculated\n")
    }
    
    if (!is.null(values$ecData)) {
      cat("✓ Event coverage calculated\n")
    }
    
    if (!is.null(values$zetaData)) {
      cat("✓ Zeta scores calculated\n")
    }
    
    if (!is.null(values$svmData)) {
      cat("✓ SVM analysis completed\n")
    }
    
    if (!is.null(values$fdrData)) {
      cat("✓ FDR cutoffs calculated\n")
    }
    
    if (!is.null(values$selectedHits)) {
      cat("✓ Hits selected:", nrow(values$selectedHits), "genes\n")
    }
    
    if (!is.null(values$singleCellData)) {
      cat("✓ Single cell QC completed\n")
    }
    
    cat("\nAnalysis ready for download!")
  })
}

# Run the application
shinyApp(ui = ui, server = server) 