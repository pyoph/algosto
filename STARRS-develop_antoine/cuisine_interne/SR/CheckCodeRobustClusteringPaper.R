rm(list=ls())
library(mclust)
library(STARRS)

################################################################################
# 4.1 1 first example
################################################################################
set.seed(1)
# Dimensions
K <- 3; nk <- rep(500, K); d <- 5
# Gaussian parameters
mu <- rbind(mu1=rep(0, d), mu2=rep(3, d), mu3=rep(-3, d))
sigma <- array(dim=c(K, d, d))
sigma[1, , ] <- 2*diag(d); sigma[2, , ] <- diag(1:d); sigma[3, , ] <- diag(1/(1:d))
# Student contamination
fraction <- 0.3
dfT <- 1; dfT <- 1; muT <- mu; sigmaT <- sigma
# Simulation
GMMcontStudent <- rGMMcontStudent(nk=nk, muG=mu, sigmaG=sigma,
                                  delta=fraction,
                                  dfT=dfT, muT=muT, sigmaT=sigmaT)

print(head(STARRS::GMMcontStudent$X))
print(head(STARRS::GMMcontStudent$Z))
print(head(STARRS::GMMcontStudent$C))

robustGMMOutput <- multipleRobustMM(GMMcontStudent$X, nclust=1:5)

plot(robustGMMOutput, graph='BIC')

plot(robustGMMOutput, graph='ICL')

plot(robustGMMOutput, graph='Two_Dim_Uncertainty')

################################################################################
# 4.2 Simulation study
################################################################################
# 4.2.2 Distance-based methods
################################################################################
data('homoGMMStudentCont_K3_d5_n1500')
# k-means with K clusters
kmeansK <- lapply(1:length(homoGMMStudentCont_K3_d5_n1500), function(f){
  kmeans(homoGMMStudentCont_K3_d5_n1500[[f]]$X, centers=K, nstart=20, iter.max=100)})
# k-means with K+1 clusters
kmeansK1 <- lapply(1:length(homoGMMStudentCont_K3_d5_n1500), function(f){
  kmeans(homoGMMStudentCont_K3_d5_n1500[[f]]$X, centers=K+1, nstart=20, iter.max=100)})
# k-medians
kmedians <- lapply(1:length(homoGMMStudentCont_K3_d5_n1500), function(f){
  kmedians(homoGMMStudentCont_K3_d5_n1500[[f]]$X)})

################################################################################
# 4.2.3 Model-based approaches
################################################################################
data('heteroGMMStudentCont_K3_d5_n1500')
# Gaussian mixture model fitting via regular EM
gmm <- lapply(1:length(heteroGMMStudentCont_K3_d5_n1500), function(f){
  mclustBIC(heteroGMMStudentCont_K3_d5_n1500[[f]]$X, G=1:5)})
# Gaussian mixture model fitting via robust EM
robustGMMOutput <- lapply(1:length(heteroGMMStudentCont_K3_d5_n1500), function(f){
  multipleRobustMM(heteroGMMStudentCont_K3_d5_n1500[[f]]$X)})
