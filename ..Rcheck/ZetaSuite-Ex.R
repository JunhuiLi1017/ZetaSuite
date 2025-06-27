pkgname <- "ZetaSuite"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ZetaSuite')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("EventCoverage")
### * EventCoverage

flush(stderr()); flush(stdout())

### Name: EventCoverage
### Title: Generate event coverage analysis and visualization for
###   alternative splicing data.
### Aliases: EventCoverage
### Keywords: ZetaSuite alternative coverage event splicing

### ** Examples

data(countMat)
data(negGene)
data(posGene)
ZscoreVal <- Zscore(countMat, negGene)
ECList <- EventCoverage(ZscoreVal, negGene, posGene, binNum=100, combine=TRUE)




cleanEx()
nameEx("FDRcutoff")
### * FDRcutoff

flush(stderr()); flush(stdout())

### Name: FDRcutoff
### Title: Determine optimal cutoff thresholds based on Screen Strength
###   analysis.
### Aliases: FDRcutoff
### Keywords: FDR ZetaSuite cutoff screen strength

### ** Examples

data(nonExpGene)
data(negGene)
data(posGene)
data(ZseqList)
data(countMat)
ZscoreVal <- Zscore(countMat, negGene)
zetaData <- Zeta(ZscoreVal, ZseqList, SVM=FALSE)
cutoffval <- FDRcutoff(zetaData, negGene, posGene, nonExpGene, combine=TRUE)




cleanEx()
nameEx("QC")
### * QC

flush(stderr()); flush(stdout())

### Name: QC
### Title: Perform quality control analysis for high-throughput screening
###   data.
### Aliases: QC
### Keywords: SSMD ZetaSuite control quality t-SNE

### ** Examples

data(countMat)
data(negGene)
data(posGene)




cleanEx()
nameEx("SVM")
### * SVM

flush(stderr()); flush(stdout())

### Name: SVM
### Title: Generate SVM decision boundaries for positive and negative
###   control separation.
### Aliases: SVM
### Keywords: SVM ZetaSuite boundary decision machine support vector

### ** Examples

data(countMat)
data(negGene)
data(posGene)
ZscoreVal <- Zscore(countMat, negGene)
ECdataList <- EventCoverage(ZscoreVal, negGene, posGene, binNum=10, combine=TRUE)




cleanEx()
nameEx("SVMcurve")
### * SVMcurve

flush(stderr()); flush(stdout())

### Name: SVMcurve
### Title: The SVM curve lines in Zeta-plot.
### Aliases: SVMcurve
### Keywords: datasets

### ** Examples

  data(SVMcurve)



cleanEx()
nameEx("Zeta")
### * Zeta

flush(stderr()); flush(stdout())

### Name: Zeta
### Title: Calculation of zeta and weighted zeta score.
### Aliases: Zeta
### Keywords: ZetaSuite alternative score splicing zeta

### ** Examples

data(ZseqList)
data(SVMcurve)
data(countMat)
data(negGene)
ZscoreVal <- Zscore(countMat, negGene)
zetaData <- Zeta(ZscoreVal, ZseqList, SVM=FALSE)




cleanEx()
nameEx("ZetaSuitSC")
### * ZetaSuitSC

flush(stderr()); flush(stdout())

### Name: ZetaSuitSC
### Title: Calculate zeta score for single cell RNA-seq quality control.
### Aliases: ZetaSuitSC
### Keywords: ZetaSuite cell control quality single

### ** Examples

data(countMatSC)




cleanEx()
nameEx("Zscore")
### * Zscore

flush(stderr()); flush(stdout())

### Name: Zscore
### Title: Z-score normalization for high-throughput screening data.
### Aliases: Zscore
### Keywords: Z-score ZetaSuite normalization

### ** Examples

data(countMat)
data(negGene)
ZscoreVal <- Zscore(countMat, negGene)
ZscoreVal[1:5, 1:5]




cleanEx()
nameEx("ZseqList")
### * ZseqList

flush(stderr()); flush(stdout())

### Name: ZseqList
### Title: The bin size for Zeta calculation.
### Aliases: ZseqList
### Keywords: datasets

### ** Examples

  data(ZseqList)



cleanEx()
nameEx("countMat")
### * countMat

flush(stderr()); flush(stdout())

### Name: countMat
### Title: Subsampled data from in-house HTS2 screening for global splicing
###   regulators.
### Aliases: countMat
### Keywords: datasets

### ** Examples

  data(countMat)



cleanEx()
nameEx("countMatSC")
### * countMatSC

flush(stderr()); flush(stdout())

### Name: countMatSC
### Title: The cell x gene matrix from single-cell RNA-seq.
### Aliases: countMatSC
### Keywords: datasets

### ** Examples

  data(countMatSC)



cleanEx()
nameEx("negGene")
### * negGene

flush(stderr()); flush(stdout())

### Name: negGene
### Title: Input negative file.
### Aliases: negGene
### Keywords: datasets

### ** Examples

  data(negGene)



cleanEx()
nameEx("nonExpGene")
### * nonExpGene

flush(stderr()); flush(stdout())

### Name: nonExpGene
### Title: Input internal negative control file.
### Aliases: nonExpGene
### Keywords: datasets

### ** Examples

  data(nonExpGene)



cleanEx()
nameEx("posGene")
### * posGene

flush(stderr()); flush(stdout())

### Name: posGene
### Title: Input positive file.
### Aliases: posGene
### Keywords: datasets

### ** Examples

  data(posGene)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
