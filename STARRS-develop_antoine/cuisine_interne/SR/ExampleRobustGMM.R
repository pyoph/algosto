# First simple example for robustGMM

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(STARRS); library(mclust); palette('R3')
simDir <- '../../STARRS/cuisine_interne/SR/'
figDir <- '../../STARRS_paper/figures/'
seed <- 1; set.seed(seed)
exportFig <- TRUE

# Dimensions
K <- 3; nk <- rep(500, K); d <- 5
# Gaussian parameters
sigma <- sigmaT <- array(dim=c(K, d, d))
sigma[1, , ] <- 2*diag(d); sigma[2, , ] <- diag(1:d); sigma[3, , ] <- diag(1/(1:d))
mu <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
dfT <- 1
# Student contamination
fraction <- 0.3
dfT <- 1; muT <- mu; sigmaT <- sigma

# Simul
sample <- rGMMcontStudent(nk=nk, muG=mu, sigmaG=sigma,
                         delta=fraction,
                         dfT=dfT, muT=muT, sigmaT=sigmaT)

# RGMM
rgmm <- multipleRobustMM(sample$X, nclust=1:5)
plot(rgmm, graph='BIC')
plot(rgmm, graph='ICL')
plot(rgmm, graph='Two_Dim')
plot(rgmm, graph='Two_Dim_Uncertainty')
plot(rgmm, graph='Two_Dim', showOutliers=TRUE)

# GMM
gmm <- Mclust(sample$X, G=1:5)
print(gmm$G)
gmmK <- Mclust(sample$X, G=K)

# Classification
table(sample$Z, rgmm$allresult[[K]]$clusters); print(adjustedRandIndex(sample$Z, rgmm$allresult[[K]]$clusters))
# table(sample$Z, rgmm$bestresult$clusters); print(adjustedRandIndex(sample$Z, rgmm$bestresult$clusters))
table(sample$Z, gmmK$classification); print(adjustedRandIndex(sample$Z, gmmK$classification))
table(sample$Z, gmm$classification); print(adjustedRandIndex(sample$Z, gmm$classification))
sample$Zwo <- sample$Z; sample$Zwo[which(sample$C==1)] <- NA
table(sample$Zwo, rgmm$bestresult$clusters); print(adjustedRandIndex(sample$Zwo, rgmm$bestresult$clusters))
table(sample$Zwo, gmmK$classification); print(adjustedRandIndex(sample$Zwo, gmmK$classification))
table(sample$Zwo, gmm$classification); print(adjustedRandIndex(sample$Zwo, gmm$classification))

# Outlier detection
rgmm$bestresult$isOutlier <- 1*((1:sum(nk))%in%rgmm$bestresult$outliers)
table(sample$C, rgmm$bestresult$isOutlier)

# # Fit GMM
# gmmBic <- Mclust(sample$X, G=1:5)
# gmmK <- Mclust(sample$X, G=K)
# gmmK1 <- Mclust(sample$X, G=K+1)
#
# # Fit RGMM
# rgmmBic <- multipleRobustMM(sample$X, nclust=1:5)
#
# # Classification with outliers
# table(sample$Z, gmmBic$classification); print(adjustedRandIndex(sample$Z, gmmBic$classification))
# table(sample$Z, gmmK$classification); print(adjustedRandIndex(sample$Z, gmmK$classification))
# table(sample$Z, gmmK1$classification); print(adjustedRandIndex(sample$Z, gmmK1$classification))
# table(sample$Z, rgmmBic$bestresult$clusters); print(adjustedRandIndex(sample$Z, rgmmBic$bestresult$clusters))
# table(sample$Z, rgmmBic$allresults[[K]]$clusters); print(adjustedRandIndex(sample$Z, rgmmBic$allresults[[K]]$clusters))
# table(sample$Z, rgmmBic$allresults[[K+1]]$clusters); print(adjustedRandIndex(sample$Z, rgmmBic$allresults[[K+1]]$clusters))
#
# # Classification without outliers
# table(sample$Zwo, gmmK$classification); print(adjustedRandIndex(sample$Zwo, gmmK$classification))
# table(sample$Zwo, gmmK1$classification); print(adjustedRandIndex(sample$Zwo, gmmK1$classification))
# table(sample$Zwo, rgmmBic$allresults[[K]]$clusters); print(adjustedRandIndex(sample$Zwo, rgmmBic$allresults[[K]]$clusters))
# table(sample$Zwo, rgmmBic$allresults[[K+1]]$clusters); print(adjustedRandIndex(sample$Zwo, rgmmBic$allresults[[K+1]]$clusters))
#
# # Output
# plot(gmmK, what='uncertainty')
# rgmmK <- rgmmBic; rgmmK$Kopt <- K; rgmmK$bestresult <- rgmmBic$allresults[[K]]
# rgmmK1 <- rgmmBic; rgmmK1$Kopt <- K+1; rgmmK1$bestresult <- rgmmBic$allresults[[K+1]]
# plot(rgmmK, graph=c('Two_Dim'), showOutliers=FALSE)
# plot(rgmmK, graph=c('Two_Dim'), showOutliers=TRUE)
