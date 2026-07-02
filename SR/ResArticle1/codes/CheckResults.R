# Check results

rm(list=ls())
source('FunctionsCriteria.R')

# Parms, data & result
load("../parms/SimParmsGrid-d10-n10000.Rdata")
mu0 <- rep(0, d)
parm <- list(mu=mu0, Sigma=Sigma0)

method <- 'MCD'; k <- 0; l <- 0.1; rho1=0.7; r <- 25; sim <- 1
load(paste0('../data_simuls-23jun2026/SimData-d10-n10000-k', k, '-l', l, '-rho', rho1, '-r', r, '-sim', sim, '.RData'))
load(paste0('../resAlgos-res3/Fit-', method, '-d10-n10000-k', k, '-l', l, '-rho', rho1, '-r', r, '-sim', sim, '.RData'))
results <- resultats

# Check
ComputeCriteria(parm, data, results)
