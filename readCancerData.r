## ============================================================================================
##
## readCancerData
##   - These scripts requests preprocessed gene expression data (log2  normalized data) as well as background information.
##     Log2 differences in gene expression between matched case-control are calculated
##     and information from background are extracted for further use (in this case: metastatic status)
##
##    The gene expression matrix has genes as rows and samples as columns
##    The background information matrix has samples as rows and variables as columns
## ============================================================================================

mapToGenes <- function(data) {
  ## Map to gene symbol
  nuIDs <- rownames(data) # the rownames in the gene expression matrix are nuIDs
  mappingInfo <- nuID2RefSeqID(nuIDs, lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)
  geneName <- as.character(mappingInfo[,3])
  geneName <- geneName[geneName!=""] 

  ## Compute mean of probes for each gene
  exprs0 <- data
  uGeneNames <- unique(geneName)
  exprs <- matrix(NA, ncol=ncol(exprs0), nrow=length(uGeneNames))
  colnames(exprs) <- colnames(exprs0)
  rownames(exprs) <- uGeneNames
  for (j in 1:length(uGeneNames))
    exprs[j,] <- colMeans(exprs0[mappingInfo[,3]==uGeneNames[j], ,drop=F])
  
  exprs
}

makeListWithData <- function(data, withoutSpread, withSpread, followUpTime, isCase) {

  data <- mapToGenes(data)

  ## Compute log2 differences and standardize data (three different standardizations)
  withoutSpread <- withoutSpread[isCase]
  withSpread <- withSpread[isCase]
  followUpTime <- followUpTime[isCase]
  ## 
  dataN1 <- data[,isCase] - (data[,!isCase]) ## case - ctrl, no standardization
  sd <- apply(dataN1, 1, sd)
  mu <- rowMeans(dataN1)
  
  ## Sort data with respect to followUpTime
  sortList <- sort.list(followUpTime)
  withoutSpread <- withoutSpread[sortList]
  withSpread <- withSpread[sortList]
  followUpTime <- followUpTime[sortList]
  dataN1 <- dataN1[,sortList]
  
  list(withoutSpread=withoutSpread,
       withSpread=withSpread,
       followUpTime=followUpTime,
       dataN1=dataN1)
}

readCancerData <- function() {

  ## === Read data and background ===
  load(file="../Data/prospectiveCancerData.R")
  data <- normdata  # normdata  is the gene expression matrix
  ## The backgound matrix:
  ## Columns must be arranged in data and rows in background so that case and ctrl of the same pair are adjacent  
  data <- data[, ind]
  background <- background[ind, ]
  ## Now the case and control from a case-control pair are found in neighboring columns, i.e. log-differences
  ## can be computed in the makeListWithData function
  
  ## === Find followUpTime ===
  n <- length(background[,"diagnosedate"]) # dates are in the format ISOdate
  followUpTime <- rep(NA, n)
  for (i in 1:n)
    if (background[i,"diagnosedate"]!="") {
      date1 <- getDate(background[i,"freezing_date"])
      date2 <- getDate(background[i,"diagnosedate"])
      followUpTime[i] <- round(difftime(date2,date1,units="days") )
    }

  ## === Define strata - withoutSpread and withSpread ===
  withSpread <- background[,"Case_ctrl"]=="case" & background[,"metastasis"]>0 # requires a variable of metastasis status
  withoutSpread <- background[,"Case_ctrl"]=="case" & background[,"metastasis"]==0

  ## === Define isCase ===
  isCase <- background[,"Case_ctrl"]=="case"

  makeListWithData(data, withoutSpread, withSpread, followUpTime, isCase)
}
