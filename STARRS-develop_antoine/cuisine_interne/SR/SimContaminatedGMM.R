# Creation of the synthetic dataset for the STARRS paper/vignette

rm(list=ls()); palette('R3')
library(STARRS)
simDir <- './'
seed <- 1; set.seed(seed) # seed = 1, 5

# Parms
variance <- 'common'
# variance <- 'different'
contamination <- 'student'
K <- 3; nk <- 500; d <- 5
nk <- rep(nk, K); n <- sum(nk)
sigma <- sigmaT <- array(dim=c(K, d, d))
sigma[1, , ] <- 2*diag(d); sigma[2, , ] <- diag(1:d); sigma[3, , ] <- diag(1/(1:d))
if(variance=='common'){for(k in 2:K){sigma[k, , ] <- sigma[1, , ]}}
mu <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
# Uniform contamination
lower <- rep(-10, d); upper <- rep(10, d)
# Student contamination
dfT <- 1
muT <- mu; sigmaT <- sigma
# Global parms
parms <- list(mu=mu, sigma=sigma, contamination=contamination)
if(contamination=='student'){parms$dfT <- dfT; parms$muT <- muT; parms$sigmaT <- sigmaT}
if(contamination=='uniform'){parms$lower <- lower; parms$upper <- upper}

# Simulation
fractions <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5); fractionNb <- length(fractions)
inliers <- rGMMcontStudent(nk=nk, muG=mu, sigmaG=sigma, delta=0, dfT=dfT, muT=muT, sigmaT=sigmaT)
head(inliers$X)
head(inliers$Z)
if(contamination=='student'){
  outliers <- rGMMcontStudent(nk=nk, muG=mu, sigmaG=sigma, delta=max(fractions), dfT=dfT, muT=muT, sigmaT=sigmaT)
}
if(contamination=='uniform'){
  outliers <- rGMMcontUniform(nk=nk, muG=mu, sigmaG=sigma, delta=max(fractions), lUnif=lower, uUnif=upper)
}
ordering <- sample(which(outliers$C==1), sum(outliers$C==1))

# Contamination
dataList <- list()
for(f in 1:fractionNb){
  fraction <- fractions[f]
  sim <- inliers
  if(fraction>0){
    out <-ordering[1:round(fractions[f]*n)]
    sim$C[out] <- 1; sim$X[out, ] <- outliers$X[out, ]; sim$Z[out] <- outliers$Z[out]
    }
  dataList[[f]] <- sim
  }
unlist(lapply(dataList, function(sim){mean(sim$C)}))

# Plot
par(mfcol=c(ceiling(sqrt(fractionNb+1)), round(sqrt(fractionNb+1))), pch=20)
for(f in 1:fractionNb){
  color <- dataList[[f]]$Z; color[which(dataList[[f]]$C==1)] <- K+1
  plot(prcomp(dataList[[f]]$X)$x[, 1:2], col=color)
}

# Exportation
simName <- 'contaminatedGMM'
simName <- paste0(simName, '-K', K, '-d', d, '-n', n, '-', contamination, '-', variance)
simFile <- paste0(simDir, simName, '.rda')
# save(parms, fractions, dataList, file=simFile)

# Fast check
library(mclust)
kmeansK <- lapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmeans(dataList[[f]]$X, centers=K, nstart=20, iter.max=100)$cluster)
  })
kmediansK <- lapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, kmedians(dataList[[f]]$X, nclust=K, ninit=20, niter=100)$bestresult$cluster)
})
gmmK <- lapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, Mclust(dataList[[f]]$X, G=K)$classification)
})
gmm <- lapply(1:length(dataList), function(f){
  adjustedRandIndex(dataList[[f]]$Z, Mclust(dataList[[f]]$X)$classification)
})
plot(fractions, kmeansK, type='b', ylim=c(0, 1))
points(fractions, kmediansK, type='b', col=2)
points(fractions, gmmK, type='b', col=3)
points(fractions, gmm, type='b', col=4)
