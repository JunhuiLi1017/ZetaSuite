#' Z-score normalization for high-throughput screening data.
#'
#' This function performs Z-score normalization on high-throughput screening data using negative control samples as reference. The Z-score transformation standardizes the data by centering and scaling each column (readout) based on the mean and standard deviation of negative control samples.
#' 
#' @param countMat A matrix of raw count data where rows represent genes/siRNAs and columns represent readouts/conditions. The matrix should have row names corresponding to gene/siRNA identifiers.
#'
#' @param negGene A data frame or matrix containing negative control gene/siRNA identifiers. The first column should contain the gene/siRNA names that match the row names in countMat.
#'
#' @details The function performs Z-score normalization as follows:
#' \enumerate{
#'   \item Extracts negative control samples from the input matrix using the identifiers provided in negGene
#'   \item For each column (readout), calculates the mean and standard deviation using only the negative control samples
#'   \item Applies Z-score transformation: \eqn{Z_{ij} = (X_{ij} - \mu_j) / \sigma_j}
#'   where \eqn{X_{ij}} is the raw value for gene \eqn{i} in readout \eqn{j}, \eqn{\mu_j} is the mean of negative controls in readout \eqn{j}, and \eqn{\sigma_j} is the standard deviation of negative controls in readout \eqn{j}
#' }
#' 
#' This normalization allows for comparison across different readouts and identifies genes/siRNAs that show significant deviation from the negative control distribution.
#' 
#' @return A Z-score normalized matrix with the same dimensions as the input countMat (excluding the Type column added during processing). Each value represents how many standard deviations away from the negative control mean that particular gene/readout combination is.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' ZscoreVal <- Zscore(countMat, negGene)
#' ZscoreVal[1:5, 1:5]
#'
#' @keywords ZetaSuite normalization Z-score
#'
#' @importFrom stats sd
#'
#' @name Zscore
#' @export Zscore

Zscore <- function (countMat, negGene) {
    countMatNeg <- countMat[rownames(countMat) %in% negGene[, 1], ]
    countMatNeg$Type <- rep("Negative", nrow(countMatNeg))
    Nmix <- countMatNeg
    rownames(Nmix) <- rownames(countMatNeg)
    matrixData <- countMat
    rownames(matrixData) <- rownames(countMat)
    outputdata <- matrix(0, nrow(matrixData), ncol(Nmix) - 1)
    for (i in 1:(ncol(Nmix) - 1)) {
        outputdata[, i] <- (((matrixData[, i]) - mean(Nmix[, i]))/sd(Nmix[, i]))
    }
    colnames(outputdata) <- colnames(Nmix)[-ncol(Nmix)]
    rownames(outputdata) <- rownames(matrixData)
    return(outputdata)
}
