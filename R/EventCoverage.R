#' Generate event coverage analysis and visualization for alternative splicing data.
#'
#' This function analyzes event coverage across Z-score thresholds and generates visualizations to compare positive and negative control samples. It calculates the proportion of readouts that exceed different Z-score thresholds for each gene, creating the foundation for zeta score calculations.
#'
#' @param ZscoreVal A matrix of Z-scores where rows represent genes and columns represent readouts/conditions. This is typically the output from the Zscore() function.
#'
#' @param negGene A data frame or matrix containing negative control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in ZscoreVal.
#'
#' @param posGene A data frame or matrix containing positive control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in ZscoreVal.
#'
#' @param binNum The number of bins to divide the Z-score range. Recommended value is 100. The function creates Z-score thresholds from the 0.00001 to 0.99999 quantiles of the data.
#'
#' @param combine Logical. Whether to combine the negative and positive Z-score ranges. Default is TRUE. When TRUE, uses symmetric ranges around zero; when FALSE, uses separate ranges for negative and positive values.
#'
#' @return A list containing two sublists:
#'   \item{ECdataList}{A list with the following components:
#'     \itemize{
#'       \item ZseqList: A data frame with Z-score thresholds for decrease (Zseq_D) and increase (Zseq_I) directions
#'       \item EC_N_I: Event coverage matrix for negative controls in increase direction
#'       \item EC_N_D: Event coverage matrix for negative controls in decrease direction
#'       \item EC_P_I: Event coverage matrix for positive controls in increase direction
#'       \item EC_P_D: Event coverage matrix for positive controls in decrease direction
#'     }
#'   }
#'   \item{ECplotList}{A list with two ggplot objects:
#'     \itemize{
#'       \item EC_jitter_D: Jitter plot showing event coverage for decrease direction
#'       \item EC_jitter_I: Jitter plot showing event coverage for increase direction
#'     }
#'   }
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Determines Z-score thresholds based on data quantiles and binNum
#'   \item For each gene and threshold, calculates the proportion of readouts that exceed (increase) or fall below (decrease) the threshold
#'   \item Separates data into negative and positive control groups
#'   \item Generates jitter plots comparing event coverage between control groups
#'   \item Returns both data matrices and visualization plots
#' }
#' 
#' The event coverage matrices can be used as input for SVM analysis and zeta score calculations.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' data(posGene)
#' ZscoreVal <- Zscore(countMat, negGene)
#' ECList <- EventCoverage(ZscoreVal, negGene, posGene, binNum=100, combine=TRUE)
#'
#' @keywords ZetaSuite event coverage alternative splicing
#'
#' @importFrom ggplot2 ggplot geom_jitter ylab xlab theme theme_bw scale_fill_manual aes element_text element_blank
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#'
#' @importFrom stats quantile
#'
#' @name EventCoverage
#' @export EventCoverage

EventCoverage <- function (ZscoreVal, negGene, posGene, binNum, combine = TRUE){
    ZscoreVal[is.na(ZscoreVal)] <- 0
	old <- options()
	on.exit(options(old))
    options(digits = 15)
    meltdata <- melt(ZscoreVal)
    if (binNum <= 0){
        stop("binNum should be more than 0")
    }
    if (combine == TRUE){
        minV = quantile(meltdata$value, probs = c(1e-05))
        maxV = quantile(meltdata$value, probs = c(0.99999))
        maxValAbs <- max(abs(minV), abs(maxV))
        stepmax = maxValAbs/binNum
        Zseq_D <- seq(maxValAbs * (-1), -1.3, abs(stepmax))
        Zseq_I <- seq(1.3, maxValAbs, abs(stepmax))
    }else {
        minV = quantile(meltdata$value, probs = c(1e-05))
        stepmin = minV/binNum
        Zseq_D <- seq(minV, 0, abs(stepmin))
        maxV = quantile(meltdata$value, probs = c(0.99999))
        stepmax = maxV/binNum
        Zseq_I <- seq(0, maxV, abs(stepmax))
    }
    ZseqDI <- data.frame(Zseq_D, Zseq_I,stringsAsFactors = FALSE)

    nColCM <- ncol(ZscoreVal)
    EC_D <- matrix(NA, nrow(ZscoreVal), length(Zseq_D))
    colnames(EC_D) <- Zseq_D
    rownames(EC_D) <- rownames(ZscoreVal)
    EC_I <- EC_D
    colnames(EC_I) <- Zseq_I
    for (i in 1:length(Zseq_D)){
        EC_D[, i] <- rowSums(ZscoreVal < Zseq_D[i])/nColCM
        EC_I[, i] <- rowSums(ZscoreVal > Zseq_I[i])/nColCM
    }

    EC_I_N <- EC_I[rownames(EC_I) %in% negGene[, 1], ]
    EC_D_N <- EC_D[rownames(EC_D) %in% negGene[, 1], ]
    EC_I_P <- EC_I[rownames(EC_I) %in% posGene[, 1], ]
    EC_D_P <- EC_D[rownames(EC_D) %in% posGene[, 1], ]
    ECDN <- melt(EC_D_N)
    ECDP <- melt(EC_D_P)
    ECIN <- melt(EC_I_N)
    ECIP <- melt(EC_I_P)

    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p_EC_jitter_D <- ggplot() + geom_jitter(aes(x = ECDN[, 2], y = ECDN[, 3]), colour = "#67a9cf", size = 0.1) + geom_jitter(aes(x = ECDP[, 2], y = ECDP[, 3]), colour = "#ef8a62", size = 0.1) + ylab("Events percent") + xlab("Zscore") + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = getPalette(22)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p_EC_jitter_I <- ggplot() + geom_jitter(aes(x = ECIN[, 2], y = ECIN[, 3]), colour = "#67a9cf", size = 0.1) + geom_jitter(aes(x = ECIP[, 2], y = ECIP[,3]), colour = "#ef8a62", size = 0.1) + ylab("Events percent") + xlab("Zscore") + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = getPalette(22)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ECdataList <- list()
    ECdataList["ZseqList"] <- list(ZseqDI)
    ECdataList["EC_N_I"] <- list(EC_I_N)
    ECdataList["EC_N_D"] <- list(EC_D_N)
    ECdataList["EC_P_I"] <- list(EC_I_P)
    ECdataList["EC_P_D"] <- list(EC_D_P)
    ECplotList <- list()
    ECplotList['EC_jitter_D'] <- list(p_EC_jitter_D)
    ECplotList['EC_jitter_I'] <- list(p_EC_jitter_I)
    return(list(ECdataList,ECplotList))
}

