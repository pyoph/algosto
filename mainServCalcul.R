library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(mvtnorm)
library(Rcpp)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
#library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
library("easystats")
library("bigutilsr")
library("RobRegression")
#library("rJava")
#library("REPPlab")
#library("RMM")
library("mvnfast")
library(ROCR)
library(pROC)
library(tidyr)
library(dplyr)
library(binom)

#library(knitr)
library("RSpectra")
library("mvoutlier")
library("RMM")
library("matrixcalc")
#library(Rcpp)
#library(RcppArmadillo)
#install.packages("Metrics")
#library(Metrics)
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")

#source("~/algosto/seuils.R")
source("~/algosto/mainNew.R")

#source("~/algosto/shrinkageCabana.R")
#sourceCpp("~/algosto/valeursVecteursPropres.cpp")
#sourceCpp("~/algosto/RobinsMC2CPP.cpp")
#sourceCpp("~/algosto/update_mean_cov.cpp")
sourceCpp("~/algosto/CodesRCpp.cpp")

res1runFarScenario = calcule_tout(nb_runs = 1,contamin = "moyenne_variance" )

res1runNearScenario = calcule_tout(nb_runs = 1,contamin = "moyenne_variance" )

res1rund100FarScenario = calcule_tout(nb_runs = 1,contamin = "moyenne_variance" )

save(res1rund100,file = "res1rund100.RData")

res10run = calcule_tout(nb_runs = 10,contamin = "moyenne" )

res10run = calcule_tout(nb_runs = 10,contamin = "moyenne_variance")
res100runNearesScenario = calcule_tout(nb_runs = 100,contamin = "moyenne_variance")

save(res100runNearesScenario,file = "res100runNeasScenario.RData")

save(res,file= "resMoyenne.Rdata")

# 
# results_metrics <- RMSEAUCFPdataset(contamin = "moyenne")
# save(results_metrics,file = "results_metricsMoyenne100runs1e4d100.Rdata")
# 
# temps_calcul = results_metrics$temps_calculBP
# 
# 
# results_metrics <- round(results_metrics$results_metrics,2)
# save(results_metrics,file = "results_metricsMaronnaZamar100runs1e4.Rdata")
# 
# 
# results_metrics <- RMSEAUCFPdataset(contamin = "UniformeTronquee")
# 
# results_metrics <- round(results_metrics$results_metrics,2)
# save(results_metrics,file = "results_metricsUniformeTronquee100runs1e4.Rdata")
# 
# 
# results_metrics <- RMSEAUCFPdataset(contamin = "studentTronquee")
# 
# results_metrics <- round(results_metrics$results_metrics,2)
# save(results_metrics,file = "results_metricsstudentTronquee100runs1e4.Rdata")
# 
