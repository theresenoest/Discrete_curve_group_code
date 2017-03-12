## Make data list
resPath <- "Res/TwoStrata/" # Requires results subfolder in directory
strataNames <- c("With spread", "Without spread")
strataCols <- 2:1 ## color of point in plots is red for with spread and black  for without spread
strataSign <- c(1,1) ## shape of point in plot is a circle both for with and without spread
dataList <- list()
dataList[[1]] <- list(dataPeriod1Spr,
                      dataPeriod2Spr,
                      dataPeriod3Spr) 
dataList[[2]] <- list(dataPeriod1NotSpr,
                      dataPeriod2NotSpr,
                      dataPeriod3NotSpr) 
names(dataList) <- strataNames
rData1 <- dataList[[1]] 
rData2 <- dataList[[2]] 
    
## Test hypothesis for comparing two strata
limitPval <- constantList$limitPval
nofSim <- constantList$nofSim
nofPeriods <- constantList$nofPeriods
nofGenes <- nrow((rData1)[[1]])
nofGroups <- factorial(nofPeriods)
curveNames <- defineCurveClasses(nofPeriods)
  
print(paste(Sys.time(), "-- Real and sim. data"))
  
## Compute covariates and t-statistics for real and simulated data
statTwoStrata <-
  computeStatisticsTwoStrata(rData1, rData2, dataList, constantList)
covVal1A <- statTwoStrata$covVal1A
covVal1B <- statTwoStrata$covVal1B
covVal2A <- statTwoStrata$covVal2A
covVal2B <- statTwoStrata$covVal2B
covValTestA <- statTwoStrata$covValTestA
covValTestB <- statTwoStrata$covValTestB
tStatA <- statTwoStrata$tStatA
tStatB <- statTwoStrata$tStatB

tStatASim <- array(NA, dim=c(nofSim, nofGroups, nofPeriods))
tStatBSim <- array(NA, dim=c(nofSim, nofGroups, nofPeriods))
for (sim in 1:nofSim) {
  if (sim%%10==1)
    print(paste("    ", Sys.time(), "-- Simulation", sim, "of", nofSim, "...."))
  sData1 <- sData2 <- list() 
  for (p in 1:nofPeriods) {
    sData <- simulateData(list(rData1[[p]], rData2[[p]]))
    sData1[[p]] <- sData[[1]]
    sData2[[p]] <- sData[[2]]
  }    
  statTwoStrata <-
    computeStatisticsTwoStrata(sData1, sData2, c(), constantList)
  tStatASim[sim,,] <- statTwoStrata$tStatA
  tStatBSim[sim,,] <- statTwoStrata$tStatB
} 

## ------------------------
## Make table with results
## ------------------------

resMat <- matrix(NA, nrow=nofGroups, ncol=2*nofPeriods)
colnames(resMat) <-
  c(paste("with",nofPeriods:1,sep=""),
    paste("without",nofPeriods:1,sep=""))
for (k in 1:nofGroups) {
  for (p in 1:nofPeriods) {
    pVal <- signif((sum(abs(tStatA[k,p])<=abs(tStatASim[,k,p]))+1)/(nofSim+1),2)
    mu1 <- round(mean((covVal1A[[p]])[k,],na.rm=T),2)
    mu2 <- round(mean((covVal2A[[p]])[k,],na.rm=T),2)
    resMat[k,nofPeriods-p+1] <- paste(pVal, " (", mu1, ",", mu2, ")", sep="")
  }
  for (p in 1:nofPeriods) {
    pVal <- signif((sum(abs(tStatB[k,p])<=abs(tStatBSim[,k,p]))+1)/(nofSim+1),2)
    mu1 <- round(mean((covVal1B[[p]])[k,],na.rm=T),2)
    mu2 <- round(mean((covVal2B[[p]])[k,],na.rm=T),2)
    resMat[k,2*nofPeriods-p+1] <- paste(pVal, " (", mu1, ",", mu2, ")", sep="")
  } 
}
resMat <- cbind(curveNames, resMat)
resMat <- rbind(colnames(resMat), resMat)
write.table(resMat, paste(resPath, "/resTableTestTwoStrataTime.txt",sep=""),
            col.names=F, row.names=F, quote=F, sep="\t")

## ------------------------
## Make plot with results 
## ------------------------

xInd <- c(2,2,4) ## Velg et tall mellom 1 og 6 for hver periode
yInd <- c(3,3,6) ## Velg et tall mellom 1 og 6 for hver periode
AorB <- c("AB","AB","AB") ## Velg AA, AB, BA eller BB for hver periode  

val <- c()
for (p in 1:nofPeriods) 
  for (i in 1:length(strataSign))
  val <- c(val, covValTestA[[p]][[i]][xInd[p],], covValTestB[[p]][[i]][yInd[p],])
xlim <- ylim <- c(min(val),max(val))

## Choose different covariates for the different periods
pdf(paste(resPath, "/plotCovariates.pdf",sep=""))
par(mar=c(3,3,1,1))
par(mgp=c(1.5,0.5,0))
ind <- 1:2
par(mfrow=c(2,2))
for (p in 0:nofPeriods) {
  if (p==0) {
    plot(c(1,1),xlim=c(0,1),ylim=c(0,1),axes=F,type="n",xlab="",ylab="")
    legend(c(0,1),c(0,1),
           legend=strataNames[ind], col=strataCols[ind], pch=strataSign[ind])
  } else {
    if (AorB[p]=="AA")
      plotPairOfCov(covValTestA, xInd[p],
                    paste("With spread", curveNames[xInd[p]], sep=" "),
                    covValTestA, yInd[p],
                    paste("With spread", curveNames[yInd[p]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[p]=="AB")
      plotPairOfCov(covValTestA, xInd[p],
                    paste("With spread", curveNames[xInd[p]], sep=" "),
                    covValTestB, yInd[p],
                    paste("Without spread", curveNames[yInd[p]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[p]=="BA")
      plotPairOfCov(covValTestB, xInd[p],
                    paste("Without spread", curveNames[xInd[p]], sep=" "),
                    covVal1TestA, yInd[p],
                    paste("With spread", curveNames[yInd[p]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[p]=="BB")
      plotPairOfCov(covValTestB, xInd[p],
                    paste("Without spread", curveNames[xInd[p]], sep=" "),
                    covValTestB, yInd[p],
                    paste("Without spread", curveNames[yInd[p]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
  }
} 
dev.off()

## Always choose covariates selected for period 1
pdf(paste(resPath, "/plotCovariatesCommonCov.pdf",sep="")) 
par(mar=c(3,3,1,1))
par(mgp=c(1.5,0.5,0))
ind <- 1:2
par(mfrow=c(2,2))
for (p in 0:nofPeriods) {
  if (p==0) {
    plot(c(1,1),xlim=c(0,1),ylim=c(0,1),axes=F,type="n",xlab="",ylab="")
    legend(c(0,1),c(0,1),
           legend=strataNames[ind], col=strataCols[ind], pch=strataSign[ind])
  } else {
    if (AorB[1]=="AA")
      plotPairOfCov(covValTestA, xInd[1],
                    paste("With spread", curveNames[xInd[1]], sep=" "),
                    covValTestA, yInd[1],
                    paste("With spread", curveNames[yInd[1]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[1]=="AB")
      plotPairOfCov(covValTestA, xInd[1],
                    paste("With spread", curveNames[xInd[1]], sep=" "),
                    covValTestB, yInd[1],
                    paste("Without spread", curveNames[yInd[1]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[1]=="BA")
      plotPairOfCov(covValTestB, xInd[1],
                    paste("Without spread", curveNames[xInd[1]], sep=" "),
                    covVal1TestA, yInd[1],
                    paste("With spread", curveNames[yInd[1]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
    else if (AorB[1]=="BB")
      plotPairOfCov(covValTestB, xInd[1],
                    paste("Without spread", curveNames[xInd[1]], sep=" "),
                    covValTestB, yInd[1],
                    paste("Without spread", curveNames[yInd[1]], sep=" "), p,
                    strataSign[ind], strataCols[ind], xlim, ylim)
  }
} 
dev.off()    


fact <- 1 
pdf(paste(resPath, "/plotOneCovariate.pdf", sep=""))
plot(c(-10,-10),xlim=c(-2.5,1.2)*fact,ylim=c(0.9,5),xaxt="n",yaxt="n",type="n",xlab="Value of variable",ylab="")
abline(v=-1.2*fact)
axis(1, at=round((-6:6)/5*fact,1))
for (p in 1:nofPeriods) {
  yval <- 4.5-(p-1)*1.5
  text(-1.9*fact, yval, paste("Period", p))
  j <- 0
  for (i in 1:2) {
    j <- j+1
    yval <- (length(dataList)-j+1)/2.5+3.5-(p-1)*1.5
    text(-1.9*fact, yval, names(dataList[i]))
    data <- c()
    if (AorB[p]=="AA" || AorB[p]=="AB")
      data <- covValTestA[[p]][[i]][xInd[p],]
    else if (AorB[p]=="BA" || AorB[p]=="BB")
      data <- covValTestB[[p]][[i]][xInd[p],]
    nc <- length(data)
    points(data, rep(yval,nc), pch="|", col=strataCols[i])
  }
}
dev.off()
