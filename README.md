# ZetaSuite

An R package for analyzing multi-dimensional high-throughput screening data, particularly two-dimensional RNAi screens and single-cell RNA sequencing data.

## Installation

```r
# Install from CRAN
install.packages("ZetaSuite")

# Or install from GitHub
devtools::install_github("username/ZetaSuite")

# Load the package
library(ZetaSuite)
```

## Quick Start

```r
# Load example data
data(countMat)
data(negGene)
data(posGene)
data(nonExpGene)

# Quality Control
qc_results <- QC(countMat, negGene, posGene)

# Z-score normalization
zscore_matrix <- Zscore(countMat, negGene)

# Event coverage analysis
ec_results <- EventCoverage(zscore_matrix, negGene, posGene)

# Zeta score calculation
zeta_scores <- Zeta(zscore_matrix, ec_results[[1]]$ZseqList)

# FDR cutoff analysis
fdr_results <- FDRcutoff(zeta_scores, negGene, posGene, nonExpGene)
```

## Interactive Shiny Application

Launch the interactive web interface for ZetaSuite:

```r
# Launch the Shiny app
ZetaSuiteApp()

# Launch without opening browser automatically
ZetaSuiteApp(launch.browser = FALSE)

# Launch on a specific port
ZetaSuiteApp(port = 3838)
```

The Shiny app provides:
- Interactive data upload and visualization
- Step-by-step analysis workflow
- Real-time results and plots
- Data export capabilities
- Built-in example dataset

## Features

- **Quality Control Analysis**: Comprehensive evaluation of experimental design and data quality
- **Z-score Normalization**: Standardization using negative controls as reference
- **Event Coverage Analysis**: Quantification of regulatory effects across thresholds
- **Zeta Score Calculation**: Area-under-curve based scoring for regulatory effects
- **SVM-based Background Correction**: Machine learning approach to filter noise
- **Screen Strength Analysis**: Optimal threshold selection for hit identification
- **Single Cell Quality Control**: Cell quality assessment for scRNA-seq data

## Documentation

For detailed documentation and examples, see the package vignette:

```r
vignette("ZetaSuite")
```

## Bug Reports

If you encounter any bugs or have feature requests, please report them on our GitHub issues page:

[Report a Bug](https://github.com/JunhuiLi1017/ZetaSuite/issues)

## Citation

If you use ZetaSuite in your research, please cite:

Hao, Y., Zhang, S., Shao, C. et al. ZetaSuite: computational analysis of two-dimensional high-throughput data from multi-target screens and single-cell transcriptomics. Genome Biol 23, 162 (2022). https://doi.org/10.1186/s13059-022-02729-4

## License

This package is licensed under the MIT License - see the LICENSE file for details.
