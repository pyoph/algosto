
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

#install.packages("pROC")
#install.packages("ROCR")



#Tracé de la courbe ROC, calcule l'AUC pour plusieurs seuils et renvoie le seuil pour lequel l'AUC est maximal

courbeROC <- function(labelsVrais,distances){
  
  seuils <- seq(1,64,by =1)
  
  #initialisation des vrais positifs et des faux positifs, et du vecteur pour stocker l'AUC à différents seuils
  tpr <- numeric(length(seuils))
  fpr <- numeric(length(seuils))
  auc <- numeric(length(seuils))
  
  seuils <- sort(seuils,decreasing = TRUE)
  
  for (s in seq_along(seuils)){
    #Calcul des outliers à partir des distances pour chaque seuil
    #s = 18
    outliers_labels <- detectionOutliers(distances, cutoff = s)
    tc <- table(labelsVrais, outliers_labels)
    tc <- safe_access_tc(tc)
    #tc
    #print(i)
    if((tc["1","1"] + tc["1","0"])!= 0)
    {tpr[s]   <- round((tc["1", "1"]/(tc["1","1"] + tc["1","0"])),2)}
    #else {faux_positifs_maha[i]   <- 0}
    if((tc["0","1"] + tc["0","0"]) != 0){
      fpr[s] <-  round((tc["0","1"]/(tc["0","0"] + tc["0","1"])),2)
    }
    #outliers_labels[1] <- 1
    #pred  <- prediction(outliers_labels,resultsSimul$labelsVrais)
    
    
    roc_obj <- roc(labelsVrais, outliers_labels)
    auc[s] <- round(auc(roc_obj)*100,2)
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
  
  #Renvoi du seuil pour lequel l'AUC est maximal et de l'AUC maximal
  
  seuil_auc_max <- as.numeric(roc_df$Seuils[which.max(roc_df$auc)])
  auc_max <- max(roc_df$auc)
  
  return(list(seuil_auc_max = seuil_auc_max, auc_max = auc_max))  

}
