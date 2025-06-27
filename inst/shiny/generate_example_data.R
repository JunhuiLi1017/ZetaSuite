# Generate Example Data for ZetaSuite Shiny App
# This script creates sample data files that can be used to test the application

# Set seed for reproducibility
set.seed(123)

# Generate example count matrix
generate_count_matrix <- function(n_genes = 1000, n_samples = 50) {
  # Create gene names
  gene_names <- paste0("Gene", 1:n_genes)
  
  # Generate count data with some structure
  # Most genes have moderate expression
  base_expression <- rnorm(n_genes, mean = 10, sd = 3)
  
  # Create sample-specific effects
  sample_effects <- rnorm(n_samples, mean = 0, sd = 2)
  
  # Generate count matrix
  count_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
  
  for (i in 1:n_genes) {
    for (j in 1:n_samples) {
      # Add some noise and ensure positive values
      count_matrix[i, j] <- max(0, base_expression[i] + sample_effects[j] + rnorm(1, 0, 1))
    }
  }
  
  # Add some genes with high expression (potential positive controls)
  high_expr_genes <- sample(1:n_genes, 20)
  for (i in high_expr_genes) {
    count_matrix[i, ] <- count_matrix[i, ] * runif(1, 2, 5)
  }
  
  # Add some genes with low expression (potential negative controls)
  low_expr_genes <- sample(setdiff(1:n_genes, high_expr_genes), 30)
  for (i in low_expr_genes) {
    count_matrix[i, ] <- count_matrix[i, ] * runif(1, 0.1, 0.5)
  }
  
  # Set row and column names
  rownames(count_matrix) <- gene_names
  colnames(count_matrix) <- paste0("Sample", 1:n_samples)
  
  return(count_matrix)
}

# Generate control gene lists
generate_control_genes <- function(count_matrix, n_neg = 30, n_pos = 20, n_non_exp = 50) {
  gene_names <- rownames(count_matrix)
  
  # Calculate mean expression for each gene
  mean_expr <- rowMeans(count_matrix)
  
  # Select negative controls (low expression)
  neg_genes <- names(sort(mean_expr))[1:n_neg]
  
  # Select positive controls (high expression)
  pos_genes <- names(sort(mean_expr, decreasing = TRUE))[1:n_pos]
  
  # Select non-expressed genes (very low expression)
  non_exp_genes <- names(sort(mean_expr))[1:n_non_exp]
  
  return(list(
    neg_genes = data.frame(Gene = neg_genes, stringsAsFactors = FALSE),
    pos_genes = data.frame(Gene = pos_genes, stringsAsFactors = FALSE),
    non_exp_genes = data.frame(Gene = non_exp_genes, stringsAsFactors = FALSE)
  ))
}

# Generate single cell count matrix
generate_single_cell_data <- function(n_cells = 500, n_genes = 2000) {
  # Create cell and gene names
  cell_names <- paste0("Cell", 1:n_cells)
  gene_names <- paste0("Gene", 1:n_genes)
  
  # Generate sparse count matrix (typical for single-cell data)
  count_matrix <- matrix(0, nrow = n_cells, ncol = n_genes)
  
  # Set row and column names
  rownames(count_matrix) <- cell_names
  colnames(count_matrix) <- gene_names
  
  # Generate expression for each cell-gene combination
  for (i in 1:n_cells) {
    # Each cell expresses a subset of genes
    n_expressed <- rpois(1, lambda = 500)  # Average 500 genes per cell
    expressed_genes <- sample(1:n_genes, min(n_expressed, n_genes))
    
    for (j in expressed_genes) {
      # Generate count from negative binomial distribution
      count_matrix[i, j] <- rnbinom(1, mu = 5, size = 1)
    }
  }
  
  # Add some high-quality cells (more genes expressed)
  high_quality_cells <- sample(1:n_cells, 50)
  for (i in high_quality_cells) {
    extra_genes <- sample(setdiff(1:n_genes, which(count_matrix[i, ] > 0)), 200)
    for (j in extra_genes) {
      count_matrix[i, j] <- rnbinom(1, mu = 3, size = 1)
    }
  }
  
  # Add some low-quality cells (fewer genes expressed)
  low_quality_cells <- sample(setdiff(1:n_cells, high_quality_cells), 100)
  for (i in low_quality_cells) {
    # Remove some expressed genes
    expressed <- which(count_matrix[i, ] > 0)
    to_remove <- sample(expressed, length(expressed) * 0.7)
    count_matrix[i, to_remove] <- 0
  }
  
  return(count_matrix)
}

# Main function to generate all example data
generate_example_data <- function() {
  cat("Generating example data for ZetaSuite Shiny App...\n")
  
  # Generate main analysis data
  cat("1. Generating count matrix...\n")
  count_matrix <- generate_count_matrix(n_genes = 1000, n_samples = 50)
  
  cat("2. Generating control gene lists...\n")
  controls <- generate_control_genes(count_matrix)
  
  # Generate single cell data
  cat("3. Generating single cell count matrix...\n")
  single_cell_matrix <- generate_single_cell_data(n_cells = 500, n_genes = 2000)
  
  # Write files
  cat("4. Writing data files...\n")
  
  # Write count matrix
  write.csv(count_matrix, "example_count_matrix.csv")
  cat("   - example_count_matrix.csv\n")
  
  # Write control gene lists
  write.csv(controls$neg_genes, "example_negative_controls.csv", row.names = FALSE)
  cat("   - example_negative_controls.csv\n")
  
  write.csv(controls$pos_genes, "example_positive_controls.csv", row.names = FALSE)
  cat("   - example_positive_controls.csv\n")
  
  write.csv(controls$non_exp_genes, "example_non_expressed_genes.csv", row.names = FALSE)
  cat("   - example_non_expressed_genes.csv\n")
  
  # Write single cell matrix
  write.csv(single_cell_matrix, "example_single_cell_matrix.csv")
  cat("   - example_single_cell_matrix.csv\n")
  
  cat("\nExample data generation complete!\n")
  cat("You can now use these files to test the ZetaSuite Shiny application.\n")
  cat("\nFile descriptions:\n")
  cat("- example_count_matrix.csv: Main count matrix for analysis\n")
  cat("- example_negative_controls.csv: Negative control genes\n")
  cat("- example_positive_controls.csv: Positive control genes\n")
  cat("- example_non_expressed_genes.csv: Non-expressed genes for FDR analysis\n")
  cat("- example_single_cell_matrix.csv: Single cell count matrix\n")
  
  # Return summary statistics
  cat("\nData summary:\n")
  cat("- Count matrix:", nrow(count_matrix), "genes ×", ncol(count_matrix), "samples\n")
  cat("- Negative controls:", nrow(controls$neg_genes), "genes\n")
  cat("- Positive controls:", nrow(controls$pos_genes), "genes\n")
  cat("- Non-expressed genes:", nrow(controls$non_exp_genes), "genes\n")
  cat("- Single cell matrix:", nrow(single_cell_matrix), "cells ×", ncol(single_cell_matrix), "genes\n")
}

# Run the data generation
if (interactive()) {
  generate_example_data()
} else {
  cat("To generate example data, run: generate_example_data()\n")
} 