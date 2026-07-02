# Check Paul's simulation files

rm(list=ls())
setwd('/home/robin/RECHERCHE/APPRENTISSAGE/P-Guillot/ResArticle1/codes/')

# Dirs
parmDir <- '../parms/'
simDir <- '../data_simuls-29jun2026/'
resDir <- '../resAlgos-29jun2026/'

# Parms
B <- 20

################################################################################
# Data
simList <- list.files(simDir); simNb <- length(simList)
dataList <- c('Z', 'labelsVrais')
# for(ss in 1:simNb){
#   if(ss%%round(sqrt(simNb))==0){cat(s, '')}
#   sim <- simList[[ss]]
#   load(paste0(simDir, sim))
#   if(!identical(names(data), dataList)){
#     cat(sim, ':', names(data), '\n')
#   }
# }
# cat('\n')

# Get parms list
d <- n <- k <- l <- rho <- r <- rep('', simNb)
for(ss in 1:simNb){
  if(ss%%round(sqrt(simNb))==0){cat(ss, '')}
  res <- gsub('.RData', '', simList[[ss]])
  words <- strsplit(res, '-')[[1]]
  d[ss] <- gsub('d', '', words[2])
  n[ss] <- gsub('n', '', words[3])
  k[ss] <- gsub('k', '', words[4])
  l[ss] <- gsub('l', '', words[5])
  rho[ss] <- gsub('rho', '', words[6])
  r[ss] <- gsub('r', '', words[7])
  if((as.numeric(r[ss]) < 1) & (as.numeric(r[ss]) > 0)){print(res)}
}
cat('\n')
unique(d); table(d)
unique(n); table(n)
unique(k); table(k)
unique(l); table(l)
unique(rho); table(rho)
unique(r); table(r)

################################################################################
# Results
resList <- list.files(resDir); resNb <- length(resList)
critList <- c('variance', 'outliers_labels', 'distances', 'temps')
# for(rr in 1:resNb){
#   if(rr%%round(sqrt(resNb))==0){cat(rr, '')}
#   res <- resList[[rr]]
#   load(paste0(resDir, res))
#   if(!identical(names(resultats), critList)){
#     cat(res, ':', names(resultats), '\n')
#   }
# }
# cat('\n')

# Get parms list
method <- d <- n <- k <- l <- rho <- r <- sim <- rep('', resNb)
for(rr in 1:resNb){
  if(rr%%round(sqrt(resNb))==0){cat(rr, '')}
  res <- gsub('.RData', '', resList[[rr]])
  words <- strsplit(res, '-')[[1]]
  method[rr] <- words[2]
  d[rr] <- gsub('d', '', words[3])
  n[rr] <- gsub('n', '', words[4])
  k[rr] <- gsub('k', '', words[5])
  l[rr] <- gsub('l', '', words[6])
  rho[rr] <- gsub('rho', '', words[7])
  r[rr] <- gsub('r', '', words[8])
  sim[rr] <- gsub('sim', '', words[9])
}
cat('\n')
methodList <- sort(unique(method)); table(method);
unique(d); table(d)
unique(n); table(n)
kList <- sort(as.numeric(unique(k))); table(k)
# kList <- sort(as.numeric(names(which(table(k)>=20))))
# lList <- sort(as.numeric(unique(l))); table(l)
lList <- c(0.01, 0.1, 0.5, 1, 2, 10, 100)
rhoList <- sort(as.numeric(unique(rho))); table(rho)
# rhoList <- sort(as.numeric(names(which(table(rho)>=20))))
rList <- sort(as.numeric(unique(r))); table(r)
# rList <- sort(as.numeric(names(which(table(r)>=20))))
unique(sim)

################################################################################
# Parms
d <- unique(d); n <- unique(n)
dimName <- paste0('-d', d, '-n', n)
save(methodList, kList, lList, rhoList, rList, file=paste0(parmDir, 'ParmsLists', dimName, '.Rdata'))

################################################################################
# Check multiple results
for(k in kList){
  for(l in lList){
    for(rho in rhoList){
      for(method in methodList){
        resName <- paste0(method, dimName, '-k', k, '-l', l, '-rho', rho)
        cat('k=', k, ' l=', l, ' rho=', rho, ' ', method, '\n', sep='')
        for(r in rList){
          # method <- methodList[1]; k <- kList[1]; l <- lList[1]; rho <- rhoList[1]; r <- rList[1]
          simName <- paste0(dimName, '-k', k, '-l', l, '-rho', rho, '-r', r)
          simFile <- paste0(simDir, 'SimData', simName, '.RData')
          if(file.exists(simFile)){
            # # 1 simul x 1 run
            # resFileList <- list.files(resDir, pattern=paste0(resName, '-r', r))
            # # if(length(resFileList)==0){cat(simName, ': no result \n')}
            # if(length(resFileList)>1){print(resFileList)}
            # 1 simul x B runs
            for(b in 1:B){
              resFileList <- list.files(resDir, pattern=paste0(resName, '-r', r, '-sim', 'b'))
              if(length(resFileList)>1){print(resFileList)}
            }
          }
        }
      }
    }
  }
}


