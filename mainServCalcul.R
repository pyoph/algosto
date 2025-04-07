
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
library("RGMM")
library("mvnfast")
library(ROCR)
library(pROC)
library(tidyr)
library(dplyr)
#library(knitr)
library("RSpectra")
library("mvoutlier")
#library(Rcpp)
#library(RcppArmadillo)
#install.packages("Metrics")
#library(Metrics)
source("~/algosto/parametres.R")
source("~/algosto/simulations.R")
source("~/algosto/algorithmes.R")
source("~/algosto/resultats.R")
source("~/algosto/Outliers.R")
source("~/algosto/computeOutliers.R")
source("~/algosto/seuils.R")
source("~/algosto/shrinkageCabana.R")
#sourceCpp("~/algosto/valeursVecteursPropres.cpp")



results_metrics <- RMSEAUCFPdataset(contamin = "moyenne")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsMaronnaZamar100runs1e4.Rdata")


results_metrics <- RMSEAUCFPdataset(contamin = "UniformeTronquee")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsUniformeTronquee100runs1e4.Rdata")


results_metrics <- RMSEAUCFPdataset(contamin = "studentTronquee")

results_metrics <- round(results_metrics$results_metrics,2)
save(results_metrics,file = "results_metricsstudentTronquee100runs1e4.Rdata")

