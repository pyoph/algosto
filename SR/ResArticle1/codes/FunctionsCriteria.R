# Function for Paul's simulation study

library(mclust); library(pROC); 

library(AUC); library(aricode)

ComputeCriteria <- function(parm, data, results){
  # cat(dim(data$Z), ',', length(data$labelsVrais), '/', dim(results$variance), ',', length(results$outliers_labels), ',', length(results$distances), '\n')
  confusion <- table(data$labelsVrais, results$outliers_labels)
  rmseMu <- sqrt(mean((parm$mu - results$m)^2))
  rmseSigma <- sqrt(mean((parm$Sigma - results$variance)^2))
  fp <- sum((1-data$labelsVrais)*results$outliers_labels) / sum((1-data$labelsVrais))
  fn <- sum(data$labelsVrais*(1-results$outliers_labels)) / sum(data$labelsVrais)
  if(sum(data$labelsVrais) > 0){
    ari <- mclust::adjustedRandIndex(data$labelsVrais, results$outliers_labels)
    # aricode::ARI(data$labelsVrais, results$outliers_labels) # same result as mclust
    roc <- as.vector(pROC::roc(data$labelsVrais, results$distance, direction='<', quiet=TRUE)$auc)
    # AUC::auc(AUC::roc(results$distance, as.factor(data$labelsVrais))) # not the same as pROC
    # auc2 <- AUC::auc(AUC::roc(results$distance, as.factor(data$labelsVrais)))
    acc <- AUC::auc(AUC::accuracy(results$distance, as.factor(data$labelsVrais)))
    diag <- sum(diag(confusion)) / sum(confusion)
  }else{ari <- roc <- diag <- acc <- NA}
  return(list(confusion=confusion, rmseMu=rmseMu, rmseSigma=rmseSigma, ari=ari, roc=roc, acc=acc, fp=fp, fn=fn, diag=diag))
}

# PlotArray <- function(plotArray, plotParms, main, yLab, yLim, Log='', hRef=NULL){
#   rList <- as.numeric(dimnames(plotArray)[[2]])
#   plot(mean(yLim), col=0, xlim=range(rList), ylim=yLim, xlab='', ylab=yLab, main=main, log=Log)
#   for(mm in 1:length(dimnames(plotArray)[[3]])){
#     points(rList, colMeans(plotArray[, , mm, drop=FALSE]), type='b',
#            col=plotParms$col[mm], 
#            lty=plotParms$lty[mm], 
#            pch=plotParms$pch[mm])
#   }
# }

PlotMatrix <- function(plotMatrix, plotParms, main, yLab, yLim, Log='', hRef=NULL){
  rList <- as.numeric(dimnames(plotMatrix)[[1]])
  plot(0, col=0, xlim=range(rList), ylim=yLim, xlab='', ylab=yLab, main=main, log=Log)
  if(!is.null(hRef)){abline(h=hRef, col='gray', lty=2)}
  for(mm in 1:ncol(plotMatrix)){
    points(rList, plotMatrix[, mm], type='b',
           col=plotParms$col[mm], 
           lty=plotParms$lty[mm], 
           pch=plotParms$pch[mm])
  }
}
