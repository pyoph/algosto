
library(xtable)
library(reshape2)
library(RobRegression)
library("robustbase")
library(mvtnorm)
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
library("RMM")
library("mvnfast")
library(ROCR)
library(pROC)
library(tidyr)
library(dplyr)
#library(knitr)
library("RSpectra")
library("mvoutlier")
library("RMM")
#library(Rcpp)
#library(RcppArmadillo)
#install.packages("Metrics")
#library(Metrics)
source("~/work/algosto/parametres.R")
source("~/work/algosto/simulations.R")
source("~/work/algosto/algorithmes.R")
source("~/work/algosto/resultats.R")
source("~/work/algosto/Outliers.R")
source("~/work/algosto/computeOutliers.R")
source("~/work/algosto/seuils.R")
source("~/work/algosto/shrinkageCabana.R")
#sourceCpp("~/algosto/valeursVecteursPropres.cpp")
#sourceCpp("~/algosto/RobinsMC2CPP.cpp")



results_metrics <- RMSEAUCFPdataset(contamin = "moyenne")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsMaronnaZamar100runs1e4.Rdata")


results_metrics <- RMSEAUCFPdataset(contamin = "UniformeTronquee")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsUniformeTronquee100runs1e4.Rdata")


results_metrics <- RMSEAUCFPdataset(contamin = "studentTronquee")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsstudentTronquee100runs1e4.Rdata")

