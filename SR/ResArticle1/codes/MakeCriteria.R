# Makes a fake example for Paul's simulation study: criteria

rm(list=ls())
setwd('/home/robin/Bureau/Paul/ResArticle1/codes')
source('FunctionsCriteria.R')

# Dirs
parmDir <- '../parms/'
simDir <- '../data_simuls-29jun2026/'
resDir <- '../resAlgos-29jun2026/'
critDir <- '../criteria-29jun2026/'

# Parms
n <- 1e4; d <- 10; B <- 20 # Actually: B = 100
dimName <- paste0('-d', d, '-n', n)
load(paste0(parmDir, 'SimParmsGrid', dimName, '.Rdata')) # Caution: wrong rList
load(paste0(parmDir, 'ParmsLists', dimName, '.Rdata'))
mu0 <- rep(0, d)
parm <- list(mu=mu0, Sigma=Sigma0)

# Criteria for all configurations
critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'auc', 'auc2')
critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'auc', 'auc2', 'diag')
critList <- c('rmseSigma', 'fn', 'fp', 'ari', 'roc', 'acc', 'diag')
rNb <- length(rList)
methodNb <- length(methodList)
# mm <- 1; method <- methodList[1]; k <- 0; l <- 0; rho <- 0.3; rr <- 3; b <- 1
for(k in kList){
  for(l in lList){
    for(rho in rhoList){
      cat('k =', k, ' l =', l, ' rho =', rho, ': \n')
      ok <- FALSE
      critName <- paste0(dimName, '-k', k, '-l', l, '-rho', rho)
      critArray <-
        array(NA, dim=c(B, methodNb, rNb, length(critList)), dimnames=list(1:B, methodList, rList, critList))
      for(rr in 1:rNb){
        r <- rList[[rr]]
        cat('- r =', r, ': ')
        for(b in 1:B){ # b <- 1
          simName <- paste0(critName, '-r', r, '-sim', b)
          simFile <- paste0(simDir, 'SimData', simName, '.RData')
          if(file.exists(simFile)){
            # if(b%%round(sqrt(B))==0){cat(b, '')}
            cat(b, '')
            load(simFile)
            for(mm in 1:methodNb){
              method <- methodList[mm]
              resName <- paste0(method, simName)
              # critFile <- paste0(critDir, 'Crit', critName, '-B', B, '.RData')
              resFileList <- list.files(resDir, pattern=paste0(resName, '.R'))
              if(length(resFileList)==1){
                ok <- TRUE
                # cat(method, '')
                load(paste0(resDir, resFileList[[1]]))
                # if(method=='OfflinewithQuantcorr'){resultats <- resOffline; resultats$outliers_labels <- resultats$outlier_labels}
                criteria <- ComputeCriteria(parm=parm, data=data, results=resultats)
                critArray[b, mm, rr, ] <- unlist(sapply(critList, function(crit){criteria[crit]}))
                # cat(method, ':\n'); print(critArrayList[[mm]][b, , ])
                }
            }
          }
          data <- c()
        }
        cat('\n')
      }
      if(ok){
        critFile <- paste0(critDir, 'Crit', critName, '-B', B, '.RData')
        save(critArray, file=critFile)
        print(critName)
        for(mm in 1:methodNb){
          print(methodList[mm]); print(apply(critArray[, mm, , ], c(2, 3), mean, na.rm=TRUE))
        }
      }
    }
  }
}
