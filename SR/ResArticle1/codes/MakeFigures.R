# Makes a fake example for Paul's simulation study: criteria

rm(list=ls())
source('FunctionsCriteria.R'); # palette('R3')

# Dirs
parmDir <- '../parms/'
critDir <- '../criteria-29jun2026/'
figDir <- '../figures-29jun2026/'
exportFig <- FALSE

# Parms
n <- 1e4; d <- 10; B <- 5
dimName <- paste0('-d', d, '-n', n)
load(paste0(parmDir, 'ParmsLists', dimName, '.Rdata'))
methodNb <- length(methodList); rNb <- length(rList)
load(paste0(parmDir, 'ScenarioList', dimName, '.Rdata'))
# critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'auc', 'auc2'); critNb <- length(critList)
# critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'auc'); critNb <- length(critList)
# critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'auc', 'diag'); critNb <- length(critList)
# critList <- c('rmseSigma', 'fn', 'fp', 'auc', 'diag'); critNb <- length(critList)
critList <- c('rmseSigma', 'fn', 'fp', 'roc', 'acc', 'diag'); critNb <- length(critList)

# Plot parms
plotParms <- as.data.frame(rbind(
  c("MCD",                                  2, "blue",      0),
  c("OfflineUswithoutQuantcorr",            1, "pink",    2),
  c("OfflinewithQuantcorr",                 1, "pink",    6),
  c("OGK",                                  2, "blue",     21),
  c("OnlineUsQuantonlinecorr",              1, "orange",      2),
  c("OnlineUswithoutQuantonlinecorr",       1, "orange",      6),
  c("Oracle",                               2, "black",    21),
  c("SampleNaiveQuantonlinecorr",           2, "green",     2),
  c("SampleNaivewithoutonlinequantilecorr", 2, "green",     6),
  c("StreamingUsonlineQuantcorr",           1, "red",       6),
  c("StreamingUswithoutQuantonlinecorr",    1, "red",       2)
  ))
colnames(plotParms) <- c('method', 'lty', 'col', 'pch')
plotParms$lty <- as.numeric(plotParms$lty)
plotParms$pch <- as.numeric(plotParms$pch)
plotParms <- plotParms[order(methodList), ]

# Plot for each scenario
scenarioNb <- length(scenarioList)
for(ss in 1:scenarioNb){ # ss <- 2
  scenario <- scenarioList[[ss]]
  main <- names(scenarioList)[ss]
  figName <- paste0(main, '-B', B)
  cat(names(scenarioList)[ss], ': ')
  if(exportFig){pdf(paste0(figDir, figName, '.pdf'))}
  par(mfrow=c(length(scenario), critNb), mex=.6)
  for(pp in 1:length(scenario)){ # pp <- 2
    cat(names(scenario)[pp], '')
    config <- scenario[[pp]]
    critName <- paste0(dimName, '-k', config[1], '-l', config[2], '-rho', config[3])
    critFile <- paste0(critDir, 'Crit', critName, '-B', B, '.RData')
    if(file.exists(critFile)){
      load(critFile)
      for(cc in 1:critNb){ # cc <- 1
        plotMatrix <- matrix(NA, rNb, methodNb); rownames(plotMatrix) <- rList; colnames(plotMatrix) <- methodList
        yLab <- critList[cc]
        for(mm in 1:methodNb){ # mm <- 1
          # plotArray[, , mm] <- critArray[, , critList[cc]] # Attention: not the same rList
          plotMatrix[, mm] <- colMeans(critArray[, mm, , critList[cc], drop=FALSE], na.rm=TRUE)
        }
        if(sum(is.na(plotMatrix)) < prod(dim(plotMatrix))){
          yLim <- hRef <- NULL; Log <- ''
          if(critList[cc]=='rmseSigma'){yLim <- range(plotMatrix, na.rm=TRUE); hRef <- 0}
          # if(critList[cc]=='rmseSigma'){yLim <- c(0.1, 10); Log='y'}
          if(critList[cc]=='fp'){yLim <- c(0, .2); hRef <- c(0, 0.05, 1)}
          if(critList[cc]=='fn'){yLim <- c(0, 1); hRef <- c(0, 1)}
          if(critList[cc]=='ari'){yLim <- c(min(min(plotMatrix, na.rm=TRUE), 0), 1); hRef <- c(0, 1)}
          if(critList[cc]=='auc'){yLim <- c(0, 1); hRef <- c(0, 0.5, 1)}
          if(critList[cc]=='diag'){yLim <- c(0, 1); hRef <- c(0, 1)}
          PlotMatrix(plotMatrix, plotParms, main, yLab, yLim, Log, hRef)
          if(critList[cc]=='diag'){lines(rList, 1-2e-4*rList*(100-rList), lty=2, col='gray')}
          main <- hRef <- yLim <- c()
        }else{cat(critList[cc], ':', sum(is.na(plotMatrix)), '==', prod(dim(plotMatrix)), '\n')}
      }
      critArray <- c()
    }else{cat(' [Missing ', critName, '] ', sep='')}
  }
  cat('\n')
  if(exportFig){dev.off()}
}

