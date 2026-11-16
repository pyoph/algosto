


library(MASS)

###########################################################################
##########Output A file for parame   ters in d = 10 and in d = 100 such 
################
###########################################################################
####################################
#####################sample size and number of runs###############
############################################################
d <- 10
n <- 1e4
rList = rList = 5*0:10
Ninit = 5*1e2
rList <- c(2,4,6,8,rList)
rList <- sort(rList)
sim = 1
simNb <- 1e2

# --- General settings ---
set.seed(12345)
batch = ncol(d)
#############################################################

#############################################################
#####################File for simulation parameters##########
#############################################################

###########################################################################
##########Output A file for parame   ters in d = 10 and in d = 100 such 
################
###########################################################################

#############################################################
#####################sample size and number of runs###############
############################################################
d <- 10
rList = rList = 5*0:10
rList <- c(2,4,6,8,rList)
rList <- sort(rList)
sim = 1
simNb <- 1e2
n <- 1e4
Ninit = 1e2
cm = 2
batch = ncol(d)


#########################################################################################
######################### Null parms for non contamination : F0##########################
#########################################################################################
mu0 <- rep(0, d)
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
rho0 <- 0.3
Sigma0 <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho0^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
parms0 <- list(mu=mu0, Sigma=Sigma0)


#############################################################

# # # Contamination parms: F1

m1 <- (-1)^(1:d)/sqrt(d)


ParmsF1 <- function(m1, k1, l1, rho1){
  d <- length(m1)
  mu1 <- k1*m1
  sigmaSq1 <- l1*sigmaSq0
  Sigma1 <- diag(sqrt(sigmaSq1)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq1))
  return(list(mu1=mu1, Sigma1=Sigma1))
}

#Contamination scenarios
scenarios_1_param = list(
  list(k = 5, l = 1, rho1 = 0.3),
  
  
  list(k = 10, l = 1, rho1 = 0.3),
  
  list(k = 50, l = 1, rho1 = 0.3),
  list(k = 0, l = 0.5, rho1 = 0.3),
  # 
  list(k = 0, l = 0.1, rho1 = 0.3),
  # 
  list(k = 0, l = 0.01, rho1 = 0.3), 
  list(k = 0, l = 2, rho1 = 0.3),
  
  list(k = 0, l = 10, rho1 = 0.3),
  
  list(k = 0, l = 1e2, rho1 = 0.3),
  list(k = 0, l = 1, rho1 = 0.5),
  
  list(k = 0, l = 1, rho1 = 0.7),
  
  list(k = 0, l = 1, rho1 = 0.95)
  
)
