## ============================================================
## Aim: Check if there are more genes in a curve group than expected
## ============================================================

strataNames <- c("WithSpread", "WithoutSpread")
strataNames2 <- c("With spread", "Without spread")
realData <- list()
realData[[1]] <-  list(dataPeriod1Spr,
                       dataPeriod2Spr,
                       dataPeriod3Spr) 
realData[[2]] <-  list(dataPeriod1NotSpr,
                       dataPeriod2NotSpr,
                       dataPeriod3NotSpr) 
names(realData) <- strataNames

## Test hypothesis for each stratum
limitPval <- constantList$limitPval
limitPvalVec <- c((1:9)*0.0001, (1:9)*0.001, (1:9)*0.01, (1:9)*0.1, 1.0)
nofSim <- constantList$nofSim
nofPeriods <- constantList$nofPeriods
nofStrata <- length(realData)
nofGenes <- nrow((realData[[1]])[[1]])
nofGroups <- factorial(nofPeriods)
curveNames <- defineCurveClasses(nofPeriods)
resMat <- resMat2 <- matrix(NA, nrow=length(curveNames)+1, ncol=nofStrata)

for (i in 1:nofStrata) {

  print(paste(Sys.time(), "-- Real and sim. data -- stratum", strataNames[i]))
  rData <- realData[[i]]
  statOneStratum <- computeStatisticsOneStratum(rData)
  curveCodeAndPval <- computeCurveCodeAndPval(statOneStratum)
  tstatPval <- curveCodeAndPval$pVal 
  curveCode <- curveCodeAndPval$curveCode
  
  tstatPvalSim <- curveCodeSim <- matrix(NA, ncol=nofGenes, nrow=nofSim)
  for (sim in 1:nofSim) {
    if (sim%%10==1)
      print(paste("    ", Sys.time(), "-- Simulation", sim, "of", nofSim, "...."))
    sData <- simulateData(rData)
    statOneStratum <- computeStatisticsOneStratum(sData)
    curveCodeAndPval <- computeCurveCodeAndPval(statOneStratum)
    tstatPvalSim[sim,] <- curveCodeAndPval$pVal 
    curveCodeSim[sim,] <- curveCodeAndPval$curveCode
  } 
  
  ## Compute p-values and number of genes for resMat and resMat2
  
  ## Global
  pGlobal <-  rep(NA, length(limitPvalVec))
  selInd <- which(limitPvalVec==limitPval)
  Tstat <- sum(tstatPval<limitPval)
  TstatExp <- mean(rowSums(tstatPvalSim<limitPval))
  for (k in 1:length(limitPvalVec)) {
    TstatCur <- sum(tstatPval<limitPvalVec[k])
    TstatSim <- rowSums(tstatPvalSim<limitPvalVec[k])    
    pGlobal[k] <- signif((sum(TstatCur<=TstatSim)+1)/(nofSim+1),2)
  }
  
  ## Each curve group
  pVal <- matrix(NA, nrow=length(limitPvalVec), ncol=nofGroups)
  Gvec <- GvecExp <- rep(NA, nofGroups)
  for (j in 1:nofGroups) {
    Gvec[j] <- sum(tstatPval<limitPval & curveCode==curveNames[j])
    GSim <- rep(NA, nofSim)
    for (sim in 1:nofSim)
      GSim[sim] <- sum(tstatPvalSim[sim,]<limitPval & curveCodeSim[sim,]==curveNames[j])
    GvecExp[j] <- mean(GSim)
    for (k in 1:length(limitPvalVec)) {
      G <- sum(tstatPval<limitPvalVec[k] & curveCode==curveNames[j])
      GSim <- rep(NA, nofSim)
      for (sim in 1:nofSim)
        GSim[sim] <- sum(tstatPvalSim[sim,]<limitPvalVec[k] & curveCodeSim[sim,]==curveNames[j])
      pVal[k,j] <- signif((sum(G<=GSim)+1)/(nofSim+1),2)        
    }
  }
    
  resMat[,i] <- paste(c(pGlobal[selInd], pVal[selInd,]), sep="")
  resMat2[,i] <- paste(c(Tstat, Gvec), " (", round(c(TstatExp, GvecExp)), ")", sep="")
  
  ## Make plots of p-values
  pdf(paste("Res/EachStratum/plotStratum", strataNames[i], ".pdf", sep=""))
  plot(limitPvalVec, pGlobal, type="l", log="xy", ylim=c(0.0001,100),
       xlab="alpha",
       ylab="p-value",
       col=1, main= strataNames2[i], cex.main=2.0, cex.axis=1.6, cex.lab=1.6)
  for (j in 1:nofGroups)
    lines(limitPvalVec, pVal[,j], col=j, lty=2)
  legend(c(0.002,0.1),c(1.05,100),
         legend=c("Global", paste("CurveGroup", curveNames)), col=c(1, 1:nofGroups), lty=c(1,rep(2,nofGroups)))
  abline(h=0.05, lty=3)
  dev.off()
}

resMat <- cbind(c("Global", curveNames), resMat)
resMat <- rbind(c("CurveGroup", strataNames), resMat)
write.table(resMat, "Res/EachStratum/resTableTestEachStratumTime.txt",
            col.names=F, row.names=F, quote=F, sep="\t")

resMat2 <- cbind(c("Global", curveNames), resMat2)
resMat2 <- rbind(c("CurveGroup", strataNames), resMat2)
write.table(resMat2, "Res/EachStratum/resTable2TestEachStratumTime.txt",
            col.names=F, row.names=F, quote=F, sep="\t")

