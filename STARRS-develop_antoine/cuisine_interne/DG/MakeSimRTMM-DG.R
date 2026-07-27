# Results from simul for robust TMM
# 5 scenarios : 1 = unif, 2 = 1 student (df 1), 3 = K students (df 1),
# 4 = 1 student (df 2), 5 = K students (df 2)

rm(list=ls()); par(pch=20)
seed <- 1; set.seed(seed)
library(mclust); library(combinat)

################################################################################
# Parms
p <- 5; K <- 3
n <- 300; # 60, 150, 300, 600, 1500
nk <- rep(round(n/K), K)
s <- 3 # scenario <- 1:5s
delta <- 0.0 # c(0, round(10^seq(-1.8, -.3, by=.25), 2))
load(paste0("cuisine_interne/DG/", 'parmRTMM-p', p, '-K', K, '.Rdata'))
parmList

parm <- parmList[[s]]
if(parm$dist=='uniform'){
  sim <- rTMMcontUniform(nk, parm$dfT0, parm$muT0, parm$sigmaT0, delta,
                             parm$lUnif, parm$uUnif)
}else{
  sim <- rTMMcontStudent(nk, parm$dfT0, parm$muT0, parm$sigmaT0, delta,
                             parm$dfT1, parm$muT1, parm$sigmaT1)
}

# Lancer robustMM(sim$X)

res <- robustMM(sim$X, K=3)
res_mul <- multipleRobustMM(sim$X)

adjustedRandIndex(sim$Z, res$clusters)
