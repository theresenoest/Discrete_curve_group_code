require(lumi)
require(lumiHumanIDMapping)
library(multtest)
library(limma)
library(splines)
library(sva)
library(glmnet)
library(glmnetcr)
library(GSVA)
library(combinat)
library(limma)

source("readCancerData.r")

## ============================================================================================
##
## Compute t statistic
##  - computeTstatistic
##    x0 and x1 are matrices
##
## ============================================================================================

computeTstatistic <- function(x0, x1) {
  n0 <- ncol(x0)
  n1 <- ncol(x1)
  m0 <- matrix(rep(rowMeans(x0),n0),ncol=n0)
  m1 <- matrix(rep(rowMeans(x1),n1),ncol=n1)
  v0 <- rowSums((x0-m0)**2)/(n0-1)
  v1 <- rowSums((x1-m1)**2)/(n1-1)
  
  (m0[,1]-m1[,1])/sqrt(v0/n0+v1/n1)
} 


## ============================================================================================
##
## Compute degrees of freedom for t-statistic
##  - computeDegreesOgFreedomTstat
##    x0 and x1 are matrices
##
## ============================================================================================

computeDfForTstat <- function(x0, x1) {
  
  n0 <- ncol(x0)
  n1 <- ncol(x1)
  m0 <- matrix(rep(rowMeans(x0),n0),ncol=n0)
  m1 <- matrix(rep(rowMeans(x1),n1),ncol=n1)
  v0 <- rowSums((x0-m0)**2)/(n0-1)
  v1 <- rowSums((x1-m1)**2)/(n1-1)

  (v0/n0+v1/n1)*(v0/n0+v1/n1) / ((v0*v0)/(n0*n0*(n0-1)) + (v1*v1)/(n1*n1*(n1-1)))
} 


## ============================================================================================
##
## ComputeStatisticsOneStratum
##
## ============================================================================================

computeStatisticsOneStratum <- function (data) {  

  nofPeriods <- length(data)
  nofGenes <- nrow(data[[1]])  
  meanPerGeneAndPeriod <- matrix(NA, nrow=nofGenes, ncol=nofPeriods)
  for (i in 1:nofPeriods) 
    meanPerGeneAndPeriod[,i] <- rowMeans(data[[i]])

  nuPerGeneAndPeriodComb <- tStatPerGeneAndPeriodComb <-
    matrix(NA, nrow=nofGenes, ncol=nofPeriods*(nofPeriods-1))
  colNames <- c()
  col <- 0
  for (i1 in 1:nofPeriods) {
    for (i2 in (1:nofPeriods)[-i1]) {
      col <- col+1
      colNames <- c(colNames, paste(i1, "<", i2, sep=""))
      tStatPerGeneAndPeriodComb[, col] <-
        computeTstatistic(data[[i1]], data[[i2]])
      nuPerGeneAndPeriodComb[, col] <-
        computeDfForTstat(data[[i1]], data[[i2]])
    }
  }
  colnames(tStatPerGeneAndPeriodComb) <- colNames
  
  list(meanPerGeneAndPeriod=meanPerGeneAndPeriod,
       tStatPerGeneAndPeriodComb=tStatPerGeneAndPeriodComb,
       nuPerGeneAndPeriodComb=nuPerGeneAndPeriodComb)
}


## ============================================================================================
##
## ComputeStatisticsTwoStrata
##
## ============================================================================================

computeStatisticsTwoStrata <- function(data1, data2, testData, constantList) {
  
  limitPval <- constantList$limitPval
  nofSim <- constantList$nofSim
  nofPeriods <- constantList$nofPeriods
  nofGroups <- factorial(nofPeriods)
  curveNames <- defineCurveClasses(nofPeriods)
  
  statOneStratum1 <- computeStatisticsOneStratum(data1)
  statOneStratum2 <- computeStatisticsOneStratum(data2)
  curveCodeAndPval1 <- computeCurveCodeAndPval(statOneStratum1)
  curveCodeAndPval2 <- computeCurveCodeAndPval(statOneStratum2)
  tstatPval1 <- curveCodeAndPval1$pVal 
  tstatPval2 <- curveCodeAndPval2$pVal 
  curveCode1 <- curveCodeAndPval1$curveCode
  curveCode2 <- curveCodeAndPval2$curveCode  
  
  covVal1A <- covVal2A <- covVal1B <- covVal2B <- covValTestA <- covValTestB <-
    list()
  for (p in 1:nofPeriods) {
    covVal1A[[p]] <- covVal1B[[p]] <-
      matrix(NA, nrow=nofGroups, ncol=ncol(data1[[p]]))
    covVal2A[[p]] <- covVal2B[[p]] <-
      matrix(NA, nrow=nofGroups, ncol=ncol(data2[[p]]))
    if (length(testData)>0) {
      covValTestA[[p]] <-  covValTestB[[p]] <- list()
      for (i in 1:length(testData)) {
        covValTestA[[p]][[i]] <-  covValTestB[[p]][[i]] <-
          matrix(NA, nrow=nofGroups, ncol=ncol(testData[[i]][[p]]))
      }
    }
  }
  
  if (FALSE) {
    curveCode1[tstatPval1>=limitPval] <- ""
    curveCode2[tstatPval2>=limitPval] <- ""
  }
  for (k in 1:nofGroups) {
    if (TRUE) {
      n1 <- n2 <- 100
      sortedPvals1 <- sort(tstatPval1[curveCode1==curveNames[k]])
      sortedPvals2 <- sort(tstatPval2[curveCode2==curveNames[k]])
      if (length(sortedPvals1)<100)
        n1 <- length(sortedPvals1)
      if (length(sortedPvals2)<100)
        n2 <- length(sortedPvals2)
      limitPval1 <- sortedPvals1[n1]
      limitPval2 <- sortedPvals2[n2]
      curveCode1[curveCode1==curveNames[k] & tstatPval1>limitPval1] <- ""
      curveCode2[curveCode2==curveNames[k] & tstatPval2>limitPval2] <- ""
    }
    
    geneIsSelA <- curveCode1==curveNames[k] ## & curveCode2!=curveNames[k]
    geneIsSelB <- curveCode2==curveNames[k] ## & curveCode1!=curveNames[k]
    for (p in 1:nofPeriods) {
      covVal1A[[p]][k,] <- colMeans(data1[[p]][geneIsSelA,])
      covVal2A[[p]][k,] <- colMeans(data2[[p]][geneIsSelA,])
      covVal1B[[p]][k,] <- colMeans(data1[[p]][geneIsSelB,])
      covVal2B[[p]][k,] <- colMeans(data2[[p]][geneIsSelB,])
      if (length(testData)>0) 
        for (i in 1:length(testData)) {
          covValTestA[[p]][[i]][k,] <- colMeans(testData[[i]][[p]][geneIsSelA,])
          covValTestB[[p]][[i]][k,] <- colMeans(testData[[i]][[p]][geneIsSelB,])
        }
    }
  }
  
  tStatA <- matrix(NA, nrow=nofGroups, ncol=nofPeriods)
  tStatB <- matrix(NA, nrow=nofGroups, ncol=nofPeriods)
  for (p in 1:nofPeriods) {
    tStatA[,p] <- rowMeans(covVal1A[[p]])-rowMeans(covVal2A[[p]])
    tStatB[,p] <- rowMeans(covVal1B[[p]])-rowMeans(covVal2B[[p]])
  }
  
  list(covVal1A=covVal1A, covVal2A=covVal2B,
       covVal1B=covVal1B, covVal2B=covVal2B,
       covValTestA=covValTestA, covValTestB=covValTestB,
       tStatA=tStatA, tStatB=tStatB)
}


## ============================================================================================
##
## plotPairOfCov
##
## ============================================================================================

plotPairOfCov <- function(covValTestx, xInd, xlabTxt,
                          covValTesty, yInd, ylabTxt,
                          period, strataSign, strataCols,
                          xlim, ylim) {
  xval <- yval <- c()
  for (i in 1:length(strataSign)) {
    xval <- c(xval, covValTestx[[period]][[i]][xInd,])
    yval <- c(yval, covValTesty[[period]][[i]][yInd,])
  }
  
  plot(covValTestx[[period]][[1]][xInd,],
       covValTesty[[period]][[1]][yInd,],
       col=strataCols[1], pch=strataSign[1],
       xlab=xlabTxt, ylab=ylabTxt, 
       xlim=xlim, ylim=ylim)
  text(0, ylim[2], paste("Period",period))
  for (i in 2:length(strataSign)) 
    points(covValTestx[[period]][[i]][xInd,],
           covValTesty[[period]][[i]][yInd,],
           col=strataCols[i], pch=strataSign[i])
}


## ============================================================================================
##
## simulateData
##
## ============================================================================================

simulateData <- function(rData) {
  nofPeriods <- length(rData)
  
  allData <- c()
  startInd <- rep(NA, nofPeriods+1)
  stopInd <- rep(NA, nofPeriods)
  startInd[1] <- 1
  for (period in 1:nofPeriods) {
    stopInd[period] <- startInd[period]+ncol(rData[[period]])-1
    startInd[period+1] <- stopInd[period]+1
    allData <- cbind(allData, rData[[period]])
  }
  nofPairs <-  ncol(allData)
  
  ## Randomize pairs, i.e. columns in the data
  sampleVec <- sample(1:nofPairs, nofPairs)  
  sData <- list()
  for (period in 1:nofPeriods) 
    sData[[period]] <-
      allData[, sampleVec[(startInd[period]):(stopInd[period])]]
  
  sData
}


## ============================================================================================
##
## defineCurveClasses
##
## ============================================================================================

defineCurveClasses <- function(nofPeriods){
  
  codeList <-  unique(permn(1:nofPeriods))
  
  curveNames <- c()
  cNames <- rep("", length(codeList))
  for (k in 1:length(codeList)) {
    curveNames <- c(curveNames, paste(codeList[[k]], collapse=""))
    sortList <- sort.list(codeList[[k]])
    cNames[k] <-  paste(sortList[1],"<", sortList[nofPeriods], sep="")
  }
  
  names(curveNames) <- cNames
  
  curveNames
}


## ============================================================================================
##
## computeCurveCodeAndPval
##
## ============================================================================================

computeCurveCodeAndPval <- function(statisticsOneStratum) {
  
  meanPerGeneAndPeriod <- statisticsOneStratum$meanPerGeneAndPeriod
  tstatPerGeneAndPeriodComb <- statisticsOneStratum$tStatPerGeneAndPeriod
  nuPerGeneAndPeriodComb <- statisticsOneStratum$nuPerGeneAndPeriodComb
  
  nofGenes <- nrow(meanPerGeneAndPeriod)
  nofPeriods <- ncol(meanPerGeneAndPeriod)
  curveNames <- defineCurveClasses(nofPeriods)
  
  nu <- tstat <- curveCode <- rep(NA, nofGenes)
  for (g in 1:nofGenes) {
    curveCodeHelp <- sort.list(sort.list(meanPerGeneAndPeriod[g,]))
    curveCode[g] <- paste(curveCodeHelp, collapse="")
    nameCmp <- (names(curveNames))[curveNames==curveCode[g]]
    tstat[g] <- tstatPerGeneAndPeriodComb[g,colnames(tstatPerGeneAndPeriodComb)==nameCmp]
    nu[g] <- nuPerGeneAndPeriodComb[g,colnames(tstatPerGeneAndPeriodComb)==nameCmp]
  }
  
  pVal <- pt(tstat,nu)
    
  list(curveCode=curveCode, pVal=pVal)
}
