library(reshape2)
library(RobRegression)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
library("matlib")
library("MASS")
library("corrplot")
library("dplyr")
library(ROCR)
#install.packages("ROCR")



#Tracé de la courbe ROC et calcul AUC pour plusieurs seuils

courbeROC <- function(labelsVrais,distances){
  
  seuils <- seq(1,64,by =1)
  
  #initialisation des vrais positifs et des faux positifs, et du vecteur pour stocker l'AUC à différents seuils
  tpr <- numeric(length(seuils))
  fpr <- numeric(length(seuils))
  auc <- numeric(length(seuils))
  
  seuils <- sort(seuils,decreasing = TRUE)
  
  for (s in seq_along(seuils)){
    #Calcul des outliers à partir des distances pour chaque seuil
    
    outliers_labels <- detectionOutliers(distances, cutoff = s)
    tc <- table(resultsSimul$labelsVrais, outliers_labels)
    tc <- safe_access_tc(tc)
    #tc
    #print(i)
    if((tc["0","1"] + tc["0","0"])!= 0)
    {tpr[s]   <- round((tc["1", "1"]/(tc["1","1"] + tc["1","0"])),2)}
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["1","0"] + tc["1","1"]) != 0){
      fpr[s] <-  round((tc["0","1"]/(tc["0","0"] + tc["0","1"])),2)
    }
    pred  <- prediction(outliers_labels, resultsSimul$labelsVrais)
    
    # Calculating Area under Curve
    perf <- performance(pred,"auc")
    auc[s] <- round(as.numeric(perf@y.values)*100,2)
  }
  #Construction de la courbe ROC 
  roc_df <- data.frame(Seuils = seuils, TPR = tpr, FPR = fpr,auc = auc)
  
  # Create the ROC plot
  roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_point(size = 2, color = "red") +  # Add points (dots)
    geom_line(color = "blue", size = 1) +  # Connect the points with lines
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = "Courbe ROC pour différents seuils",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal()
  
  print(roc_plot)
  
}