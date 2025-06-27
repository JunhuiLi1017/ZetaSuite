#' Calculate zeta score for single cell RNA-seq quality control.
#'
#' This function evaluates the quality of cells detected in single-cell RNA-seq data by calculating a zeta score for each cell. The zeta score is based on the distribution of gene expression across different expression thresholds. A cutoff value is automatically determined using a two-component Gaussian mixture model to separate high-quality cells from low-quality or damaged cells.
#'
#' @param countMatSC A matrix of single-cell RNA-seq count data where rows represent cells and columns represent genes.
#'
#' @param binNum The number of bins for zeta score calculation. Default is 10. The function creates expression thresholds from 0 to the 80th percentile of non-zero expression values, divided into binNum intervals.
#'
#' @param filter Logical. Whether to filter out cells with total read counts less than 100. Default is TRUE. This helps remove extremely low-quality cells before analysis.
#'
#' @return A list containing:
#'   \item{zetaData}{A data frame with two columns: 'Cell' (cell identifiers) and 'Zeta' (calculated zeta scores)}
#'   \item{p_cutoff}{A ggplot object showing the distribution of log10-transformed zeta scores with fitted Gaussian mixture components and the determined cutoff threshold}
#'
#' @details The function works as follows:
#' \enumerate{
#'   \item Filters cells based on total read count if filter=TRUE
#'   \item Samples a subset of cells and genes for computational efficiency
#'   \item Creates expression thresholds (bins) from 0 to the 80th percentile of non-zero expression values
#'   \item For each cell, counts how many genes exceed each threshold
#'   \item Calculates the zeta score as a weighted sum of these counts
#'   \item Fits a two-component Gaussian mixture model to log10-transformed zeta scores
#'   \item Determines an optimal cutoff to separate high-quality from low-quality cells
#' }
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMatSC)
#' \donttest{zetaDataSC <- ZetaSuitSC(countMatSC, binNum=50, filter=TRUE)}
#'
#' @keywords ZetaSuite single cell quality control
#'
#' @importFrom ggplot2 ggplot geom_density geom_area xlim geom_vline theme_bw geom_text aes
#' @importFrom reshape2 melt
#' @importFrom mixtools normalmixEM
#'
#' @importFrom stats density dnorm
#'
#' @name ZetaSuitSC
#' @export ZetaSuitSC

ZetaSuitSC <- function(countMatSC,binNum=10,filter=TRUE){
  if(filter==TRUE){
    nCount <- which(rowSums(countMatSC)>100)
    inFilter <- countMatSC[nCount,]
  }else{
    inFilter <- countMatSC
  }
  #set.seed(1234)
  minvalue1 <- min(10000,floor(nrow(inFilter)*0.5))
  minvalue2 <- min(10000,floor(ncol(inFilter)*0.5))
  sampleVec <- inFilter[sample(nrow(inFilter),minvalue1),sample(ncol(inFilter),minvalue2)]
  Meltsample <- melt(sampleVec)
  MeltsamMt0 <- Meltsample[Meltsample$value>0,]
  maxVal <- round(sort(MeltsamMt0$value)[round(nrow(MeltsamMt0)*0.8)]/binNum)*binNum
  if(maxVal==0){
    Zseq <- seq(0,(binNum-1),1)
  }else{
    Zseq <- seq(0,maxVal,maxVal/binNum)
  }
  nFeatureDiffCutOff <- rownames(inFilter)
  for(i in Zseq){
    nFeatureDiffCutOff <- cbind(nFeatureDiffCutOff,data.frame(rowSums(inFilter > i),stringsAsFactors=FALSE))
  }
  colnames(nFeatureDiffCutOff) <- c("Cell",Zseq)
  if(binNum==10){
    binNum <- binNum - 1
  }
  allNum <- nFeatureDiffCutOff[,-1]
  zetaVal <- data.frame(rowSums(allNum[,2:binNum])+(allNum[,1]+allNum[,1+binNum])/2,stringsAsFactors=FALSE)
  zetaData <- cbind(nFeatureDiffCutOff[,1],zetaVal)
  colnames(zetaData) <- c("Cell","Zeta")

  #set.seed(1234)
  lowPop<-0.0001
  densityRes<-cbind(density(log10(zetaData$Zeta),bw=0.08)$x,density(log10(zetaData$Zeta),bw=0.08)$y)
  densityRes<-as.data.frame(densityRes,stringsAsFactors=FALSE)
  mid<-quantile(log10(zetaData$Zeta),probs=0.5)
  mean1<-densityRes[densityRes$V2==max(densityRes[densityRes$V1<mid,]$V2),]$V1
  mean2<-densityRes[densityRes$V2==max(densityRes[densityRes$V1>mid,]$V2),]$V1
  out<-normalmixEM(log10(zetaData$Zeta),k=2,epsilon = 1e-03,mean.constr=c(mean1,mean2))
  x <- seq(0,6,length.out = 1000)
  y1 <- dnorm(x, mean1,out$sigma[1])
  y2 <- dnorm(x, mean2,out$sigma[2])
  y1data<-as.data.frame(cbind(x,y1),stringsAsFactors=FALSE)
  cutoff<-round(10^round(y1data[y1data$y1<=lowPop & y1data$x>mean1 & y1data$x<mean2,][1,1],1),0)
  if(is.na(cutoff)){
    mid=quantile(log10(zetaData$Zeta),probs=0.1)
    mean1<-densityRes[densityRes$V2==max(densityRes[densityRes$V1<mid,]$V2),]$V1
    mean2<-densityRes[densityRes$V2==max(densityRes[densityRes$V1>mid,]$V2),]$V1
    out<-normalmixEM(log10(zetaData$Zeta),k=2,epsilon = 1e-03,mean.constr=c(mean1,mean2))
    x <- seq(0,6,length.out = 1000)
    y1 <- dnorm(x, mean1,out$sigma[1])
    y2 <- dnorm(x, mean2,out$sigma[2])
    y1data<-as.data.frame(cbind(x,y1),stringsAsFactors=FALSE)
    cutoff<-round(10^round(y1data[y1data$y1<=lowPop & y1data$x>mean1,][1,1],1),0)
    if(is.na(cutoff)){
      y2data<-as.data.frame(cbind(x,y2),stringsAsFactors=FALSE)
      cutoff<-round(10^round(y2data[y2data$y2>=lowPop,][1,1],1),0)
    }
  }
  string<-paste("cut-off is ", cutoff)
  p_cutoff <- ggplot()+geom_density(aes(log10(zetaData$Zeta)))+geom_area(aes(x=x,y=y1),fill="orange",alpha=0.3)+geom_area(aes(x=x,y=y2),fill="red",alpha=0.3)+xlim(0,6)+geom_vline(xintercept=log10(cutoff),linetype="dashed")+theme_bw()+geom_text(aes(x=log10(cutoff)+0.71, label=string, y=0.8))
  return(list(zetaData,p_cutoff))
}
