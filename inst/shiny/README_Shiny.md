# ZetaSuite Shiny Application

A comprehensive web-based interface for the ZetaSuite R package, providing an intuitive way to analyze high-throughput screening data and perform single-cell RNA-seq quality control.

## Features

- **Interactive Web Interface**: User-friendly dashboard with tabbed navigation
- **Built-in Example Dataset**: Explore the HTS2 screening dataset included with the package
- **Complete Analysis Workflow**: Step-by-step analysis from QC to hit selection
- **Interactive Visualizations**: Plotly-powered interactive plots with zoom, pan, and hover
- **Data Upload Support**: Upload your own CSV files for analysis
- **Results Download**: Export all analysis results as CSV files
- **Real-time Progress**: Progress indicators and status updates
- **Error Handling**: Comprehensive error messages and notifications

## Analysis Modules

### 1. Quality Control Analysis
- Score distribution plots
- t-SNE visualization for sample separation
- Box plots comparing control groups
- SSMD (Strictly Standardized Mean Difference) distribution

### 2. Z-score Normalization
- Standardization using negative controls
- Interactive data table preview
- Download normalized results

### 3. Event Coverage Analysis
- Quantification of regulatory effects across thresholds
- Separate analysis for decrease and increase directions
- Interactive jitter plots

### 4. Zeta Score Calculation
- Area-under-curve based scoring
- Optional SVM background correction
- Top hits identification by direction

### 5. SVM Analysis
- Machine learning-based decision boundaries
- Background correction for improved accuracy
- Separate analysis for decrease and increase directions

### 6. FDR Cutoff Analysis
- Screen Strength calculation
- Optimal threshold selection
- Interactive hit selection with customizable thresholds

### 7. Single Cell Quality Control
- Cell quality assessment for scRNA-seq data
- Gaussian mixture model-based cutoff determination
- Quality score distribution visualization

## Installation

### Prerequisites

Make sure you have R installed (version 3.6 or higher) and the following packages:

```r
# Install required packages
install.packages(c("shiny", "shinydashboard", "DT", "plotly", "shinyjs"))

# Install ZetaSuite package dependencies
install.packages(c("RColorBrewer", "Rtsne", "e1071", "ggplot2", "reshape2", 
                   "gridExtra", "mixtools"))
```

### Setup

1. Install the ZetaSuite package
2. Place the `app.R` file in your working directory
3. Run the Shiny application

## Usage

### Starting the Application

```r
# Load the ZetaSuite package
library(ZetaSuite)

# Run the Shiny app
shiny::runApp("app.R")
```

### Quick Start with Example Data

1. **Load Example Data**: Click the "Example Data" tab and press "Load Example Data"
2. **Quality Control**: Go to the "Quality Control" tab and run QC analysis
3. **Z-score Analysis**: Calculate normalized Z-scores
4. **Event Coverage**: Generate event coverage analysis
5. **Zeta Score**: Calculate regulatory effect scores
6. **FDR Cutoff**: Determine optimal thresholds
7. **Hit Selection**: Select significant hits based on Screen Strength
8. **Download Results**: Export all analysis results

### Using Your Own Data

#### Data Format Requirements

1. **Count Matrix (CSV)**:
   - Rows: Genes/siRNAs
   - Columns: Readouts/conditions
   - First column should contain gene/siRNA identifiers
   - Example:
     ```
     Gene,Readout1,Readout2,Readout3
     Gene1,10.5,12.3,8.9
     Gene2,15.2,14.1,16.7
     ```

2. **Negative Control Genes (CSV)**:
   - First column: Gene/siRNA identifiers matching count matrix
   - Example:
     ```
     Gene
     NegCtrl1
     NegCtrl2
     NegCtrl3
     ```

3. **Positive Control Genes (CSV)**:
   - First column: Gene/siRNA identifiers matching count matrix
   - Example:
     ```
     Gene
     PosCtrl1
     PosCtrl2
     PosCtrl3
     ```

4. **Non-expressed Genes (CSV, optional)**:
   - First column: Gene/siRNA identifiers matching count matrix
   - Used for FDR cutoff analysis

#### Single Cell Data

**Single Cell Count Matrix (CSV)**:
- Rows: Cells
- Columns: Genes
- First column should contain cell identifiers
- Example:
  ```
  Cell,Gene1,Gene2,Gene3
  Cell1,5,0,12
  Cell2,0,8,3
  Cell3,15,2,0
  ```

### Analysis Workflow

1. **Data Upload**: Upload your data files in the "Data Upload" tab
2. **Quality Control**: Run QC analysis to assess data quality
3. **Z-score Analysis**: Normalize your data using negative controls
4. **Event Coverage**: Calculate event coverage across thresholds
5. **Zeta Score**: Compute regulatory effect scores
6. **SVM Analysis** (optional): Generate decision boundaries
7. **FDR Cutoff**: Determine optimal thresholds
8. **Hit Selection**: Select significant hits based on Screen Strength
9. **Single Cell QC** (if applicable): Quality control for single-cell data
10. **Results**: Download all analysis results

### Parameter Settings

#### Event Coverage Parameters
- **Number of Bins**: Number of Z-score thresholds (default: 100)
- **Combine Directions**: Whether to combine decrease and increase directions

#### Zeta Score Parameters
- **Use SVM Curves**: Whether to use SVM curves for background correction

#### FDR Cutoff Parameters
- **Combine Directions**: Whether to combine decrease and increase directions
- **Screen Strength Threshold**: Minimum SS value for hit selection (default: 0.8)

#### Single Cell QC Parameters
- **Number of Bins**: Number of expression thresholds (default: 10)
- **Filter Low Count Cells**: Remove cells with <100 total reads

## Example Dataset

The application includes the HTS2 screening dataset from the ZetaSuite package:

- **1,609 genes × 100 alternative splicing events**
- **30 negative control genes** (non-specific siRNAs)
- **20 positive control genes** (PTB-targeting siRNAs)
- **50 non-expressed genes** (RPKM < 1 in HeLa cells)

This dataset demonstrates the complete analysis workflow and can be used to explore the application's features.

## Output Files

The application generates several output files:

1. **Z-score Results**: Normalized Z-score matrix
2. **Zeta Scores**: Regulatory effect scores for each gene
3. **FDR Results**: Cutoff thresholds and Screen Strength values
4. **Selected Hits**: Genes identified as significant based on Screen Strength
5. **Single Cell Results**: Cell quality scores and cutoff values

## Interactive Features

- **Data Tables**: Sortable and searchable data tables with horizontal scrolling
- **Interactive Plots**: Zoom, pan, and hover over plotly visualizations
- **Progress Indicators**: Real-time progress updates for long-running analyses
- **Error Handling**: Comprehensive error messages and notifications
- **Status Updates**: Real-time status updates for each analysis step

## Troubleshooting

### Common Issues

1. **Package Not Found**: Make sure all required packages are installed
2. **Data Format Errors**: Ensure CSV files have correct column headers and data types
3. **Memory Issues**: For large datasets, consider reducing the number of bins or filtering data
4. **SVM Analysis Fails**: Ensure you have sufficient positive and negative control samples

### Performance Tips

- For large datasets (>10,000 genes), consider preprocessing to reduce computational time
- Use the "Filter Low Count Cells" option for single-cell data to remove low-quality cells
- Adjust the number of bins based on your data size and computational resources
- Use `combine = TRUE` for faster event coverage analysis

### Error Messages

- **"Error loading data"**: Check file format and ensure all required files are uploaded
- **"Error in Quality Control"**: Verify that positive and negative control genes are present in the count matrix
- **"Error in Z-score calculation"**: Ensure negative control genes are properly formatted
- **"Error in Event Coverage"**: Check that Z-scores have been calculated first

## Advanced Usage

### Custom Analysis Workflows

The application supports custom analysis workflows:

1. **SVM Background Correction**: Enable SVM analysis for improved accuracy
2. **Separate Direction Analysis**: Analyze decrease and increase directions separately
3. **Custom Thresholds**: Adjust Screen Strength thresholds for different sensitivity/specificity trade-offs
4. **Batch Processing**: Process multiple datasets by uploading different files

### Integration with R Scripts

The Shiny app can be integrated with custom R scripts:

```r
# Load results from Shiny app
zscore_results <- read.csv("zscore_results.csv", row.names = 1)
zeta_results <- read.csv("zeta_scores.csv", row.names = 1)

# Continue with custom analysis
# ...
```

## Support

For issues related to the ZetaSuite package functionality, please refer to the main package documentation. For Shiny app-specific issues:

1. Check the R console for error messages
2. Ensure all dependencies are properly installed
3. Verify data format requirements
4. Check the troubleshooting section above

## References

- **Software**: Hao, Y., Shao, C., Zhao, G., Fu, X.D. (2021). ZetaSuite: A Computational Method for Analyzing Multi-dimensional High-throughput Data, Reveals Genes with Opposite Roles in Cancer Dependency. *Forthcoming*

- **Dataset**: Shao, C., Hao, Y., Qiu, J., Zhou, B., Li, H., Zhou, Y., Meng, F., Jiang, L., Gou, L.T., Xu, J., Li, Y., Wang, H., Yeo, G.W., Wang, D., Ji, X., Glass, C.K., Aza-Blanc, P., Fu, X.D. (2021). HTS2 Screen for Global Splicing Regulators Reveals a Key Role of the Pol II Subunit RPB9 in Coupling between Transcription and Pre-mRNA Splicing. *Cell. Forthcoming*

## License

This Shiny application is provided under the same license as the ZetaSuite package (MIT + file LICENSE). 