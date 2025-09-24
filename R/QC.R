#' Perform quality control analysis for high-throughput screening data.
#'
#' This function performs comprehensive quality control analysis on high-throughput screening data to evaluate experimental design and data quality. It generates multiple diagnostic plots and calculates SSMD (Strictly Standardized Mean Difference) scores to assess the separation between positive and negative controls.
#'
#' @param countMat A matrix of raw count data where rows represent genes/siRNAs and columns represent readouts/conditions. The matrix should have row names corresponding to gene/siRNA identifiers.
#'
#' @param negGene A data frame or matrix containing negative control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in countMat.
#'
#' @param posGene A data frame or matrix containing positive control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in countMat.
#'
#' @return A list containing four diagnostic plots:
#'   \item{score_qc}{A jitter plot showing the distribution of raw scores across all readouts for positive and negative controls}
#'   \item{tSNE_QC}{A t-SNE plot showing the global separation of positive and negative control samples in 2D space}
#'   \item{QC_box}{Side-by-side boxplots showing the distribution of scores for positive and negative controls across all readouts}
#'   \item{QC_SSMD}{A density plot showing the distribution of SSMD scores across readouts, with a threshold line at SSMD=2 and the percentage of high-quality readouts displayed}
#'
#' @details The function performs the following quality control analyses:
#' \enumerate{
#'   \item Creates jitter plots to visualize score distributions across readouts
#'   \item Performs t-SNE dimensionality reduction to assess global sample separation
#'   \item Generates boxplots to compare score distributions between control groups
#'   \item Calculates SSMD scores for each readout: \eqn{\mathrm{SSMD} = (\mu_{pos} - \mu_{neg}) / \sqrt{\sigma_{pos}^2 + \sigma_{neg}^2}}
#'   \item Reports the percentage of readouts with \eqn{|\mathrm{SSMD}| \ge 2} (considered high quality)
#' }
#' 
#' SSMD scores \eqn{\ge 2} indicate good separation between positive and negative controls, suggesting high-quality readouts.
#'
#' @references
#' Laurens van der Maaten & Geoffrey Hinton: Visualizing Data using t-SNE. Journal of Machine Learning Research 2008, 9(2008):2579-2605.
#'
#' Zhang XD: A pair of new statistical parameters for quality control in RNA interference high-throughput screening assays. Genomics 2007, 89:552-561.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' data(posGene)
#' \donttest{QC(countMat, negGene, posGene)}
#'
#' @keywords ZetaSuite quality control SSMD t-SNE
#'
#' @importFrom ggplot2 ggplot geom_jitter scale_color_manual theme_bw theme xlab ylab geom_point geom_boxplot geom_hline ggtitle geom_density geom_vline annotate ylim aes_string element_blank
#' @importFrom reshape2 melt
#' @importFrom Rtsne Rtsne
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#'
#' @importFrom stats quantile var
#'
#' @name QC
#' @export QC

QC <- function (countMat, negGene, posGene){
    countMatNeg <- countMat[rownames(countMat) %in% negGene[,1], ]
    countMatNeg$Type <- rep("Negative", nrow(countMatNeg))
    countMatPos <- countMat[rownames(countMat) %in% posGene[,1], ]
    countMatPos$Type <- rep("Positive", nrow(countMatPos))
    countMatCom <- rbind(countMatNeg, countMatPos)
    countMatCom <- countMatCom[, -1]
    countMatCom[is.na(countMatCom)] <- 0
    countMatComNoNA <- countMatCom
    meltdata <- melt(countMatComNoNA, id = "Type")
    meltdata$value <- as.numeric(as.character(meltdata$value))
    p_score_qc <- ggplot(meltdata) + geom_jitter(aes_string(x = "variable", y = "value", col = "Type"), size = 0.1) + scale_color_manual(values = c("#5aae61", "#c2a5cf")) + theme_bw() + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("")
    Labels <- countMatComNoNA$Type
    #set.seed(42)
    tsne <- Rtsne(countMatComNoNA[, seq(1, ncol(countMatComNoNA) - 1)], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 10000)
    tSNEdata <- as.data.frame(cbind(tsne$Y, countMatComNoNA$Type),stringsAsFactors = FALSE)
    tSNEdata$V3 <- as.factor(tSNEdata$V3)
    tSNEdata$V1 <- as.numeric(tSNEdata$V1)
    tSNEdata$V2 <- as.numeric(tSNEdata$V2)
    p_tSNE_QC <- ggplot(tSNEdata) + geom_point(aes_string(x = "V1", y = "V2", col = "V3"), size = 1) + theme_bw() + scale_color_manual(labels = c("Negative", "Positive"), values = c("#5aae61", "#c2a5cf")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) + xlab("tSNE-1") + ylab("tSNE-2")
    Negative <- countMatComNoNA[countMatComNoNA$Type == "Negative",]
    meltNegative <- melt(Negative, id = "Type")
    Positive <- countMatComNoNA[countMatComNoNA$Type == "Positive",]
    meltPositive <- melt(Positive, id = "Type")
    p1 <- ggplot2::ggplot(meltNegative, aes_string(x = "variable", y = "value")) + geom_boxplot(outlier.shape = NA, col = "#5aae61") + geom_hline(yintercept = 0, col = "red3", linetype = "dashed") + theme_bw() + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Readouts") + ylab("Input score") + ggtitle("Negative")
    p2 <- ggplot2::ggplot(meltPositive, aes_string(x = "variable", y = "value")) + geom_boxplot(outlier.shape = NA, col = "#c2a5cf") + geom_hline(yintercept = 0, col = "red3", linetype = "dashed") + theme_bw() + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Readouts") + ylab("Input score") + ggtitle("Positive")
    p_QC_box <- gridExtra::grid.arrange(p1, p2, nrow = 1)
    SSMD <- matrix(0, (length(Negative[1, ]) - 1), 1)
    for (i in seq(1, (length(Negative[1, ]) - 1))) {
        SSMD[i, 1] = (mean(Positive[, i]) - mean(Negative[, i]))/sqrt(var(Positive[, i]) + var(Negative[, i]))
    }
    SSMD <- as.data.frame(SSMD,stringsAsFactors = FALSE)
    row.names(SSMD) <- colnames(countMatComNoNA)[1:(length(Negative[1, ]) - 1)]
    percentage <- length(SSMD[abs(SSMD[, 1]) >= 2, ])/length(SSMD[, 1])
    labels = paste("SSMD>=2(%) is", percentage, "")
    p_QC_SSMD <- ggplot2::ggplot(SSMD) + geom_density(aes_string(x = abs(SSMD[, 1])), fill = "#43a2ca", alpha = 0.5) + geom_vline(xintercept = 2, col = "red", linetype = "dashed") + theme_bw() + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("SSMD values") + ylab("Density") + annotate("text", x = quantile(SSMD[, 1], probs = 0.75), y = 0.75, label = labels) + ylim(0, 1)
    resultList <- list()
    resultList["score_qc"] <- list(p_score_qc)
    resultList["tSNE_QC"] <- list(p_tSNE_QC)
    resultList["QC_box"] <- list(p_QC_box)
    resultList["QC_SSMD"] <- list(p_QC_SSMD)
    return(resultList)
}
