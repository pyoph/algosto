#############################################################
#####################File for simulation parameters##########
#############################################################


rm(list=ls())

###########################################################################
##########Output A file for paramaters in d = 10 and in d = 100############
###########################################################################



#############################################################
#####################sample size and number of runs###############
############################################################
#d <- 10
rList <- 5*(0:10)
#load(paste0('SimParmsGrid-d', d, '.Rdata'))
simNb <- 100
n <- 1e4
#############################################################

# KL divergence btw F0 and F1


# Dims and KL distance reached. The parameters are chosen such that the KL distance is in {0,1,10,100}
#d <- 100; KLval <- c(0, 1, 10, 100)
d <- 10; KLval <- c(0, 1, 10, 100)

# Functions
KL <- function(parms1, parms2){
  invSigma2 <- solve(parms2$Sigma)
  0.5*(log(det(parms2$Sigma)/det(parms1$Sigma)) - d + sum(diag(invSigma2%*%parms1$Sigma)) +
         t(parms2$mu-parms1$mu)%*%invSigma2%*%(parms2$mu-parms1$mu))[1, 1]
}


# Contamination parms: F1
ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}

#########################################################################################
######################### Null parms for non contamination : F0##########################
#########################################################################################
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)

#######################################################################################
#################################Choice of parameters for F_1##########################
#######################################################################################



# Contamination
m1 <- (-1)^(1:d)/sqrt(d)
#KL(parms1=parms0, parms2=ParmsF1(m1, 1, 1, 0.5))

# Setting simulation parms
k1grid <- c(0, 2^seq(-4, 5, by=.001)); # k1val <- c(0, 2, 5, 10)
KLk1 <- sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))})
k1val <- k1grid[sapply(KLval, function(kl){which.min(abs(KLk1 - kl))})]
k1val
l1grid <- 2^seq(-5, 0, .01); # l1val <- c(1/20, .2, .5, 1, 2, 20)
KLl1 <- sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))})
l1val <- l1grid[sapply(KLval, function(kl){which.min(abs(KLl1 - kl))})]
l1grid <- 2^seq(0, 5, .01); # l1val <- c(1/20, .2, .5, 1, 2, 20)
KLl1 <- sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))})
l1val <- c(l1val, l1grid[sapply(KLval, function(kl){which.min(abs(KLl1 - kl))})])
l1val <- unique(l1val)
l1val
l1grid <- 2^seq(-5, 5, 0.01)
rho1grid <- seq(0.005, 0.995, by=0.005); #= rho1val <- c(0, .3, .6, 0.8, .95)
KLrho1 <- sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))})
rho1val <- rho1grid[sapply(KLval, function(kl){which.min(abs(KLrho1 - kl))})]
rho1val


#####################################################################################################################
################################# Export of the parameter file one for d = 10 and one for d=100######################
#####################################################################################################################
save(d, KLval, 
     mu0, sigmaSq0, Sigma0, rho0, 
     m1, 
     k1val, l1val, rho1val, n,d,rList,
     file=paste0('SimParmsGrid-n', n,'-d',d, '.Rdata'))

########################################################### 
##################################Plots#############################################################################
################################################################
par(mfrow=c(3, 1))
plot(k1grid, sapply(k1grid, function(k1){KL(parms0, ParmsF1(m1, k1, 1, rho0))}), log='xy')
abline(h=KLval); abline(v=k1val)
plot(l1grid, sapply(l1grid, function(l1){KL(parms0, ParmsF1(m1, 0, l1, rho0))}), log='xy')
abline(h=KLval); abline(v=l1val)
plot(rho1grid, sapply(rho1grid, function(rho1){KL(parms0, ParmsF1(m1, 0, 1, rho1))}), log='y')
abline(h=KLval); abline(v=rho1val)