source("tools.r")

dataList <- readCancerData()

# Create three periods to compare. 
# Here chosen to be the first year, the second year and all other years
# The duration of the periods should also consider sample size in each period
isPeriod1 <- dataList$followUpTime<=365
isPeriod2 <- dataList$followUpTime<=365*2 & !isPeriod1
isPeriod3 <- !isPeriod1 & !isPeriod2

# In this analyses, two strata are considered: 
# Spr = With spread, NotSpr = without spread
dataPeriod1Spr <- dataList$dataN1[,dataList$withSpread & isPeriod1]
dataPeriod1NotSpr <- dataList$dataN1[,dataList$withoutSpread & isPeriod1]
dataPeriod2Spr <- dataList$dataN1[,dataList$withSpread & isPeriod2]
dataPeriod2NotSpr <- dataList$dataN1[,dataList$withoutSpread & isPeriod2]
dataPeriod3Spr <- dataList$dataN1[,dataList$withSpread & isPeriod3]
dataPeriod3NotSpr <- dataList$dataN1[,dataList$withoutSpread & isPeriod3]

set.seed(100)
constantList <-
  list(nofPeriods=3, ## Number of periods
       nofSim=1000, ## Number of simulations
       limitPval=0.01,
       df=3) ## Degrees og freedom

# The two following scritps analyses curve groups within and across each stratum (here: with and without spread)
source("testEachStratumTime.r")
source("testTwoStrataTime.r")
