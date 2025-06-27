#' Calculation of zeta and weighted zeta score.
#'
#' This function calculates zeta scores for genes based on their Z-score profiles across different thresholds. The zeta score quantifies the regulatory effect of gene knockdown on alternative splicing events by measuring the area under the curve of event coverage across Z-score thresholds.
#'
#' @param ZscoreVal A matrix of Z-scores where rows represent genes and columns represent readouts/conditions. This is typically the output from the Zscore() function.
#'
#' @param ZseqList A list containing two vectors: 'Zseq_D' (decrease direction thresholds) and 'Zseq_I' (increase direction thresholds). These define the Z-score bins for calculating event coverage.
#'
#' @param SVMcurve Optional. A matrix containing SVM curve data for decrease and increase directions. Required only when SVM=TRUE. The matrix should have 4 columns: Z-score and coverage for decrease direction (columns 1-2), and Z-score and coverage for increase direction (columns 3-4).
#'
#' @param SVM Logical. Whether to use SVM curves for background correction. Default is FALSE. When TRUE, the function subtracts SVM-predicted background from the event coverage before calculating zeta scores.
#'
#' @return A data frame with two columns:
#'   \item{Zeta_D}{Zeta score for decrease direction (exon skipping events)}
#'   \item{Zeta_I}{Zeta score for increase direction (exon inclusion events)}
#' 
#' Each row corresponds to a gene, and the zeta scores represent the cumulative regulatory effect across all Z-score thresholds.
#'
#' @details The function calculates zeta scores as follows:
#' \enumerate{
#'   \item For each Z-score threshold, calculates the proportion of readouts that exceed (increase) or fall below (decrease) the threshold
#'   \item Computes the area under the event coverage curve using trapezoidal integration
#'   \item If SVM=TRUE, subtracts SVM-predicted background coverage before area calculation
#'   \item Returns separate scores for decrease (Zeta_D) and increase (Zeta_I) directions
#' }
#' 
#' Higher zeta scores indicate stronger regulatory effects on alternative splicing.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(ZseqList)
#' data(SVMcurve)
#' data(countMat)
#' data(negGene)
#' ZscoreVal <- Zscore(countMat, negGene)
#' zetaData <- Zeta(ZscoreVal, ZseqList, SVM=FALSE)
#'
#' @keywords ZetaSuite zeta score alternative splicing
#'
#' @name Zeta
#' @export Zeta

Zeta <- function(ZscoreVal,
                 ZseqList,
                 SVMcurve=NULL,
                 SVM=FALSE){

  Zseq_D <- ZseqList$Zseq_D
  Zseq_I <- ZseqList$Zseq_I
  ZscoreVal[is.na(ZscoreVal)] <- 0
  nColZ <- ncol(ZscoreVal)
  outputdata_D <- rep(0,nrow(ZscoreVal))
  outputdata_I <- rep(0,nrow(ZscoreVal))

  if(SVM==TRUE){
    if(is.null(SVMcurve)==TRUE){
      stop("SVMcurve should not be NULL when SVM is TRUE")
    }
    SVMcurveD <- SVMcurve[,1:2]
    SVMcurveI <- SVMcurve[,3:4]
    for (j in 1:(length(Zseq_D)-1)){
      lengthUse_D <- rowSums(ZscoreVal < Zseq_D[j])
      lengthUse_D_add <- rowSums(ZscoreVal < Zseq_D[j+1])
      conD <- (lengthUse_D/nColZ+lengthUse_D_add/nColZ-SVMcurveD[j,2]-SVMcurveD[j+1,2])*Zseq_D[j+1]*(Zseq_D[j]-Zseq_D[j+1])/2
      conD[conD < 0] <- 0
      outputdata_D<-outputdata_D+conD
      lengthUse_I <- rowSums(ZscoreVal > Zseq_I[j])
      lengthUse_I_add <- rowSums(ZscoreVal > Zseq_I[j+1])
      conI <- (lengthUse_I/nColZ+lengthUse_I_add/nColZ-SVMcurveI[j,2]-SVMcurveI[j+1,2])*Zseq_I[j+1]*(Zseq_I[j+1]-Zseq_I[j])/2
      conI[conI <0] <- 0
      outputdata_I <- outputdata_I + conI
    }
  }else{
    for (j in 1:(length(Zseq_D)-1)){
      lengthUse_D <- rowSums(ZscoreVal < Zseq_D[j])
      lengthUse_D_add <- rowSums(ZscoreVal < Zseq_D[j+1])
      outputdata_D <- outputdata_D + (lengthUse_D/nColZ+lengthUse_D_add/nColZ)*Zseq_D[j+1]*(Zseq_D[j]-Zseq_D[j+1])/2
      lengthUse_I <- rowSums(ZscoreVal > Zseq_I[j])
      lengthUse_I_add <- rowSums(ZscoreVal > Zseq_I[j+1])
      outputdata_I <- outputdata_I + (lengthUse_I/nColZ+lengthUse_I_add/nColZ)*Zseq_I[j+1]*(Zseq_I[j+1]-Zseq_I[j])/2
    }
  }
  outputdata_D <- data.frame(outputdata_D,stringsAsFactors = FALSE)
  colnames(outputdata_D)<-c("Zeta_D")
  outputdata_I <- data.frame(outputdata_I,stringsAsFactors = FALSE)
  colnames(outputdata_I)<-c("Zeta_I")

  zetaDat <- cbind(outputdata_D,outputdata_I)
  return(zetaDat)
}


