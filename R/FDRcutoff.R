#' Determine optimal cutoff thresholds based on Screen Strength analysis.
#'
#' This function calculates optimal cutoff thresholds for identifying significant hits in high-throughput screening data using Screen Strength (SS) analysis. It evaluates the trade-off between sensitivity and specificity by calculating the ratio of apparent FDR to baseline FDR across different zeta score thresholds.
#'
#' @param zetaData A data frame containing zeta scores calculated by the Zeta() function. Should have columns 'Zeta_D' and 'Zeta_I' representing decrease and increase direction scores, respectively.
#'
#' @param negGene A data frame or matrix containing negative control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in zetaData.
#'
#' @param posGene A data frame or matrix containing positive control gene/siRNA identifiers. The first column should contain gene/siRNA names that match the row names in zetaData.
#'
#' @param nonExpGene A data frame or matrix containing non-expressed gene/siRNA identifiers. These genes are used to estimate the baseline false discovery rate. The first column should contain gene/siRNA names that match the row names in zetaData.
#'
#' @param combine Logical. Whether to combine decrease and increase direction zeta scores. Default is FALSE. When TRUE, uses the sum of Zeta_D and Zeta_I; when FALSE, analyzes each direction separately.
#'
#' @return A list containing:
#'   \item{FDR_cutOff}{A data frame with 6 columns:
#'     \itemize{
#'       \item Cut_Off: Zeta score threshold
#'       \item aFDR: Apparent false discovery rate at this threshold
#'       \item SS: Screen Strength = 1 - (aFDR / bFDR)
#'       \item TotalHits: Total number of hits at this threshold
#'       \item Num_nonExp: Number of non-expressed genes among hits
#'       \item Type: Direction ("Decrease", "Increase", or "Combine")
#'     }
#'   }
#'   \item{plotList}{A list with two ggplot objects:
#'     \itemize{
#'       \item Zeta_type: Jitter plots showing zeta score distributions by gene type
#'       \item SS_cutOff: Screen Strength curves showing SS vs. zeta score threshold
#'     }
#'   }
#'
#' @details The function performs the following analysis:
#' \enumerate{
#'   \item Categorizes genes into types: "Gene" (test genes), "Positive" (positive controls), "NS_mix" (negative controls), and "non_exp" (non-expressed genes)
#'   \item Calculates baseline FDR (bFDR) as the proportion of non-expressed genes in the entire dataset
#'   \item For each zeta score threshold, calculates apparent FDR (aFDR) as the proportion of non-expressed genes among hits
#'   \item Computes Screen Strength: SS = 1 - (aFDR / bFDR)
#'   \item Generates plots showing zeta score distributions and SS curves
#' }
#' 
#' Higher Screen Strength values indicate better separation between true hits and false positives. Users can select appropriate thresholds based on desired sensitivity/specificity trade-offs.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(nonExpGene)
#' data(negGene)
#' data(posGene)
#' data(ZseqList)
#' data(countMat)
#' ZscoreVal <- Zscore(countMat, negGene)
#' zetaData <- Zeta(ZscoreVal, ZseqList, SVM=FALSE)
#' cutoffval <- FDRcutoff(zetaData, negGene, posGene, nonExpGene, combine=TRUE)
#'
#' @keywords ZetaSuite FDR cutoff screen strength
#'
#' @importFrom ggplot2 ggplot geom_jitter theme_bw theme xlab scale_color_manual geom_point geom_smooth ylab aes_string element_blank
#'
#' @importFrom grDevices dev.off pdf
#'
#' @name FDRcutoff
#' @export FDRcutoff

FDRcutoff <- function (zetaData, negGene, posGene, nonExpGene, combine = FALSE) {
  zetaData$type <- rep("Gene", nrow(zetaData))
  zetaData[rownames(zetaData) %in% negGene[, 1], "type"] <- "NS_mix"
  zetaData[rownames(zetaData) %in% posGene[, 1], "type"] <- "Positive"
  zetaData[rownames(zetaData) %in% nonExpGene[, 1], "type"] <- "non_exp"
  nonGene_D <- zetaData[zetaData$type %in% c("non_exp", "Gene"), ]

  if (combine == FALSE) {
    maxD <- sort(nonGene_D[, "Zeta_D"], decreasing = TRUE)[10]
    maxI <- sort(nonGene_D[, "Zeta_I"], decreasing = TRUE)[10]
    minD <- sort(nonGene_D[, "Zeta_D"], decreasing = FALSE)[1]
    minI <- sort(nonGene_D[, "Zeta_I"], decreasing = FALSE)[1]
    stepD <- (maxD - minD)/100
    stepI <- (maxI - minI)/100
    sum1 <- nrow(nonGene_D)
    sum2 <- sum(nonGene_D[, "type"] %in% "non_exp")
    iFDR <- sum2/(sum1 - 1)
    seqD <- seq(minD, maxD, stepD)
    FDR_cutOff_de <- matrix(NA, length(seqD), 6)
    index <- 1
    for (num in seqD) {
      totalNum <- sum(zetaData[, "Zeta_D"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_D"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff_de[index, ] <- c(num, FDR_Nexp,screen_Stress,totalNum, numNexp, "Decrease")
      index <- index + 1
    }
    seqI <- seq(minI, maxI, stepI)
    FDR_cutOff_in <- matrix(NA, length(seqI), 6)
    index <- 1
    for (num in seqI) {
      totalNum <- sum(zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff_in[index, ] <- c(num, FDR_Nexp, screen_Stress, totalNum, numNexp, "Increase")
      index <- index + 1
    }
    FDR_cutOff <- rbind(FDR_cutOff_de, FDR_cutOff_in)
  } else {
    maxD <- sort(nonGene_D[, 1] + nonGene_D[, 2], decreasing = TRUE)[10]
    minD <- sort(nonGene_D[, 1] + nonGene_D[, 2], decreasing = FALSE)[1]
    stepD <- (maxD - minD)/100
    sum1 <- nrow(nonGene_D)
    sum2 <- sum(nonGene_D[, 3] %in% "non_exp")
    iFDR <- sum2/(sum1 - 1)
    seqD <- seq(minD, maxD, stepD)
    FDR_cutOff <- matrix(NA, length(seqD), 6)
    index <- 1
    for (num in seqD) {
      totalNum <- sum(zetaData[, "Zeta_D"] + zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_D"] + zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff[index, ] <- c(num, FDR_Nexp, screen_Stress, totalNum, numNexp, "Combine")
      index <- index + 1
    }
  }
  FDR_cutOff <- as.data.frame(FDR_cutOff,stringsAsFactors = FALSE)
  colnames(FDR_cutOff) <- c("Cut_Off","aFDR","SS","TotalHits","Num_nonExp","Type")
  FDR_cutOff$SS <- as.numeric(FDR_cutOff$SS)
  FDR_cutOff$Cut_Off <- as.numeric(FDR_cutOff$Cut_Off)
  zetaData_NS <- zetaData[zetaData$type != "NS_mix", ]

  p1 <- ggplot(zetaData_NS) + geom_jitter(aes_string(x = "type", y = "Zeta_D", col = "type")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + scale_color_manual(values = c("#c994c7", "#67a9cf", "#ef8a62"))
  p2 <- ggplot(zetaData_NS) + geom_jitter(aes_string(x = "type", y = "Zeta_I", col = "type")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + scale_color_manual(values = c("#c994c7", "#67a9cf", "#ef8a62"))
  p_Zeta_type <- gridExtra::grid.arrange(p1, p2, nrow = 1)

  if (combine == FALSE) {
    fdtss <- FDR_cutOff[FDR_cutOff$SS < 0.9, ]
    Dec <- fdtss[fdtss$Type == "Decrease", ]
    Inc <- fdtss[fdtss$Type == "Increase", ]

    p1 <- ggplot(Dec, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    p2 <- ggplot(Inc, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    p_SS_cutOff <- gridExtra::grid.arrange(p1, p2, nrow = 1)

  } else {
    fdtss <- FDR_cutOff[FDR_cutOff$SS < 1, ]

    p_SS_cutOff <- ggplot(fdtss, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
  }
  plotList <- list()
  plotList['Zeta_type'] <- list(p_Zeta_type)
  plotList['SS_cutOff'] <- list(p_SS_cutOff)
  return(list(FDR_cutOff,plotList))
}
