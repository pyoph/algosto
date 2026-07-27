# Analysis of the synthetic dataset for the STARRS paper/vignette

rm(list=ls())
library(STARRS); library(mclust); library(ggplot2); palette('R3')
simDir <- './'
figDir <- '../../../STARRS_paper/figures/'
exportFig <- TRUE

# Parms
K <- 3; n <- 1500; d <- 5
contamination <- 'student'
simName <- paste0('contaminatedGMM-K', K, '-d', d, '-n', n, '-', contamination)
outliers <- TRUE
par(mfcol=c(3, 2), lwd=2, pch=20)

#-------------------------------------------------------------------------------
# Comparing kmeans & kmedians
#-------------------------------------------------------------------------------
load(paste0(simDir, simName, '-common.rda'))
# load(paste0(simDir, simName, '-common-seed1.rda'))
# Fit
kmeansK <- lapply(1:length(dataList), function(f){kmeans(dataList[[f]]$X, centers=K, nstart=20, iter.max=100)})
kmedians <- lapply(1:length(dataList), function(f){kmedians(dataList[[f]]$X)})
kmeansK1 <- lapply(1:length(dataList), function(f){kmeans(dataList[[f]]$X, centers=K+1, nstart=20, iter.max=100)})
# ARI
if(!outliers){for(f in 1:length(dataList)){dataList[[f]]$Z[which(dataList[[f]]$C==1)] <- NA}}
kmeansKari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmeansK[[f]]$cluster)
})
# K clusters
plot(fractions, kmeansKari, type='b', col=1, ylim=c(0, 1))
kmediansKari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmedians[[f]]$allresults[[K]]$cluster)
})
points(fractions, kmediansKari, type='b', col=2)
# K+1 clusters
kmeansK1ari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmeansK1[[f]]$cluster)
})
plot(fractions, kmeansK1ari, type='b', col=1, ylim=c(0, 1))
kmediansK1ari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmedians[[f]]$allresults[[K+1]]$cluster)
})
points(fractions, kmediansK1ari, type='b', col=2)
# K free
kmediansAri <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmedians[[f]]$bestresult$cluster)
})
points(fractions, kmediansAri, type='b', col=2, lty=2)
plot(fractions, sapply(1:length(dataList), function(f){kmedians[[f]]$K}), type='b', col=2, ylim=c(1, 2*K))
abline(h=K, col=8, lty=2)

#-------------------------------------------------------------------------------
# Comparing GMM & RGMM
#-------------------------------------------------------------------------------
load(paste0(simDir, simName, '-different.rda'))
# load(paste0(simDir, simName, '-different-seed1.rda'))
# Fit
gmm <- lapply(1:length(dataList), function(f){mclustBIC(dataList[[f]]$X, G=1:5)})
# # Re-run mclustICL: not found how to get it from the mclustBIC output easily
# gmmIcl <- lapply(1:length(dataList), function(f){mclustICL(dataList[[f]]$X, G=1:5)})
rgmm <- lapply(1:length(dataList), function(f){multipleRobustMM(dataList[[f]]$X, nclust=1:5)})
# ARI
if(!outliers){for(f in 1:length(dataList)){dataList[[f]]$Z[which(dataList[[f]]$C==1)] <- NA}}
# K clusters
gmmKari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, apply(mclustModel(dataList[[f]]$X, gmm[[f]], G=K)$z, 1, which.max))
})
plot(fractions, gmmKari, type='b', col=1, ylim=c(0, 1))
rgmmKari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, rgmm[[f]]$allresults[[K]]$cluster)
})
points(fractions, rgmmKari, type='b', col=2)
# K+1 clusters
gmmK1ari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, apply(mclustModel(dataList[[f]]$X, gmm[[f]], G=K+1)$z, 1, which.max))
})
plot(fractions, gmmK1ari, type='b', col=1, ylim=c(0, 1))
rgmmK1ari <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, rgmm[[f]]$allresults[[K+1]]$cluster)
})
points(fractions, rgmmK1ari, type='b', col=2)
# K BIC
gmmBicK <- sapply(1:length(dataList), function(f){mclustModel(dataList[[f]]$X, gmm[[f]])$G})
gmmBicAri <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, apply(mclustModel(dataList[[f]]$X, gmm[[f]])$z, 1, which.max))
})
points(fractions, gmmBicAri, type='b', col=1, lty=2)
rgmmBicK <- sapply(1:length(dataList), function(f){rgmm[[f]]$Kopt})
rgmmBicAri <- sapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, rgmm[[f]]$bestresult$cluster)
})
points(fractions, rgmmBicAri, type='b', col=2, lty=2)
# # K ICL
# gmmIclK <- sapply(1:length(dataList), function(f){mclustModel(dataList[[f]]$X, gmmIcl[[f]])$G})
# gmmIclAri <- sapply(1:length(dataList), function(f){
#   adjustedRandIndex(dataList[[f]]$Z, apply(mclustModel(dataList[[f]]$X, gmmIcl[[f]])$z, 1, which.max))
# })
# points(fractions, gmmIclAri, type='b', col=1, lty=3)
# rgmmIcl <- rgmm
# for(f in 1:length(dataList)){
#   rgmmIcl[[f]]$Kopt <- which.max(rgmmIcl[[f]]$ICL)
#   rgmmIcl[[f]]$bestresult <- rgmmIcl[[f]]$allresults[[rgmmIcl[[f]]$Kopt]]
# }
# rgmmIclAri <- sapply(1:length(dataList), function(f){
#   adjustedRandIndex(dataList[[f]]$Z, rgmmIcl[[f]]$bestresult$cluster)
# })
# points(fractions, rgmmIclAri, type='b', col=2, lty=3)
# rgmmIclK <- sapply(1:length(dataList), function(f){rgmmIcl[[f]]$Kopt})

plot(fractions, gmmBicK, type='b', col=1, ylim=c(1, 2*K), lty=2)
points(fractions, rgmmBicK, type='b', col=2, lty=2)
# points(fractions, gmmIclK, type='b', col=1, lty=23)
# points(fractions, rgmmIclK, type='b', col=2, lty=3)
abline(h=K, col=8, lty=2)

#-------------------------------------------------------------------------------
# Global plot
#-------------------------------------------------------------------------------
# Global
# if(exportFig){png(paste0(simName, '-detail.png'))}
par(mfcol=c(3, 2), lwd=2, pch=20)

plot(fractions, kmeansKari, type='b', col=1, ylim=c(0, 1))
points(fractions, kmediansKari, type='b', col=2)
plot(fractions, kmeansK1ari, type='b', col=1, ylim=c(0, 1))
points(fractions, kmediansK1ari, type='b', col=2)
points(fractions, kmediansAri, type='b', col=2, lty=2)
plot(fractions, sapply(1:length(dataList), function(f){kmedians[[f]]$K}), type='b', col=2, ylim=c(1, 2*K), lty=2)
abline(h=K, col=8, lty=2)

plot(fractions, gmmKari, type='b', col=1, ylim=c(0, 1))
points(fractions, rgmmKari, type='b', col=2)
plot(fractions, gmmK1ari, type='b', col=1, ylim=c(0, 1))
points(fractions, rgmmK1ari, type='b', col=2)
points(fractions, gmmBicAri, type='b', col=1, lty=2)
points(fractions, rgmmBicAri, type='b', col=2, lty=2)
# points(fractions, gmmIclAri, type='b', col=1, lty=3)
# points(fractions, rgmmIclAri, type='b', col=2, lty=3)

plot(fractions, gmmBicK, type='b', col=1, ylim=c(1, 2*K), lty=2)
points(fractions, rgmmBicK, type='b', col=2, lty=2)
# points(fractions, gmmIclK, type='b', col=1, lty=3)
# points(fractions, rgmmIclK, type='b', col=2, lty=3)
abline(h=K, col=8, lty=2)

# if(exportFig){dev.off()}

# k-means / k-medians
if(exportFig){png(paste0(figDir, simName, '-kmeans-kmedians-ARI.png'), width=480, height=360)}
plot(fractions, kmeansKari, type='b', pch=25, col=1, ylim=c(0, 1), lwd=2,
     xlab='contamination fraction', ylab='ARI')
points(fractions, kmeansK1ari, type='b', pch=25, col=1, lty=2, lwd=2)
points(fractions, kmediansKari, type='b', pch=24, col=2, lwd=2)
points(fractions, kmediansK1ari, type='b', pch=24, col=2, lty=2, lwd=2)
points(fractions, kmediansAri, type='b', pch=24, col=2, lty=3, lwd=2)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, simName, '-kmedians-kHat.png'), width=480, height=360)}
plot(fractions, sapply(1:length(dataList), function(f){kmedians[[f]]$K}), type='b',
     pch=24, col=2, lwd=2, lty=3, ylim=c(2, 2*K),
     xlab='contamination fraction', ylab='selected K')
abline(h=K, col=8, lty=2)
if(exportFig){dev.off()}

# gmm / rgmm
if(exportFig){png(paste0(figDir, simName, '-gmm-rgmm-ARI.png'), width=480, height=360)}
plot(fractions, gmmKari, type='b', pch=25, col=1, ylim=c(0, 1), lwd=2,
     xlab='contamination fraction', ylab='ARI')
points(fractions, gmmK1ari, type='b', pch=25, col=1, lty=2, lwd=2)
points(fractions, gmmBicAri, type='b', pch=24, col=1, lty=3, lwd=2)
points(fractions, rgmmKari, type='b', pch=24, col=2, lwd=2)
points(fractions, rgmmK1ari, type='b', pch=24, col=2, lty=2, lwd=2)
points(fractions, rgmmBicAri, type='b', pch=24, col=2, lty=3, lwd=2)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, simName, '-gmm-rgmm-kHat.png'), width=480, height=360)}
plot(fractions, gmmBicK, type='b',
     pch=25, col=1, lwd=2, lty=3, ylim=c(2, 2*K),
     xlab='contamination fraction', ylab='selected K')
points(fractions, rgmmBicK, type='b', pch=24, col=2, lwd=2, lty=3)
abline(h=K, col=8, lty=2)
if(exportFig){dev.off()}

