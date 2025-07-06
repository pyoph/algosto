##############
#Library needed
################

library(MASS)
library(ggplot2)

#######################################################
#Directories paths à adapter selon votre configuration
#########################################################
paramDir = "C:/Users/Paul/Documents/Simus/ParamsSim"
dataDir =  "C:/Users/Paul/Documents/Simus/DataSim"


########################################################
#Loading of the parameters
#######################################################
load(paste0('SimParmsGrid-n',n,'-d', d, '.Rdata'))




####################
#Simulate data with different contamination rates
####################

#Extraction des paramètres
kList = k1val
lList = l1val
rho1List = rho1val


#######################################################################
#Simulation of a dataset with different parameters####################
######################################################################

genererEchantillon <- function(n, d, mu1, mu2,Sigma1, Sigma2,r) {
  # InitialSisation
  n1 <- floor((1 -r/100) * n)  # Taille du groupe non contaminé
  n2 <- n - n1         # Taille du groupe contaminé
  
  labels_mu1 <- rep(0, n1)  # Labels pour les vecteurs avec moyenne mu1
  labels_mu2 <- rep(1, n2)  # Labels pour les vecteurs avec moyenne mu2
  
  if (r > 0) {
    # Générer les vecteurs selon le type de contamination
      vecteurs_mu1 <- mvrnorm(n1, mu1, Sigma1)
      vecteurs_mu2 <- mvrnorm(n2, mu2, Sigma2)
    
      # Combinaison des vecteurs
      Z <- rbind(vecteurs_mu1, vecteurs_mu2)
      labelsVrais <- c(labels_mu1, labels_mu2)  
    }
    
    
    
  else {
    # Pas de contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    labelsVrais <- rep(0, n)
  } # Mélanger aléatoirement les données
    set.seed(123)  # Pour garantir la reproductibilité
    indices <- sample(nrow(Z))
    Z <- Z[indices, ]
    labelsVrais <- labelsVrais[indices]
  
  
  return(list(Z = Z,labelsVrais = labelsVrais))
}


simNb = 10


# Simulation of datas
for(r in rList){
  #print(paste0("r =",r))
  
  for(k in kList){for(l in lList){for(rho1 in rho1List){
  for(sim in 1:simNb){
    set.seed(sim)
    
    parmsFile <- paste0('Parms-d', d, '-n', n,  '-k', k, '-l', l, '-rho1', rho1, '-r', r,'-sim', sim,".RData")
    setwd(paramDir)
    load(parmsFile)
    setwd(dataDir)               
    m1 <- (-1)^(1:d)/sqrt(d)
    sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
    SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
    dataFile <- paste0('SimData-d', d, '-n', n,  '-k', k, '-l', l, '-rho', rho1, '-r',r,'-sim', sim,".RData")
    print(paste0("nb outliers ",sum(data$labelsVrais == 1)))
    if(!file.exists(dataFile)){
      data = genererEchantillon(n,d,mu1 = mu0, mu2 = k*m1,Sigma1 = Sigma0,Sigma2 = l*SigmaContamin,r)
      save(data, file=dataFile)
      print(paste0("n-",n,"d-",d,"k-",k,"l-",l,"r-",r," save ok"))
      #print("save OK")
      }
    else {print(paste0("n-",n,"d-",d,"k-",k,"l-",l,"r-",r," file exists"))
      load(dataFile)}  
    # data <- Simul(d=d, r=r, k=k, rho1=rho1)
    
  }
    
    
}}}}


# 
# ##################################################
# #Plots
# #################################################
# 
# 
# 
# 
# afficherContaminationScenarios = function(k,l,rho,contamination,rate){
#   
#   m1 <- (-1)^(1:d)/sqrt(d)
#   sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
#   SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(rho1^(0:(d-1))) %*% diag(sqrt(sigmaSq0))
#   
#   data = genererEchantillon(n,d,mu1 = mu0, mu2 = k*m1,Sigma1 = Sigma0,Sigma2 = l*SigmaContamin,r)
#   
#   Z = data$Z
#   # Création d'un dataframe pour ggplot
#   df <- data.frame(X1 = data$Z[,1], X2 = data$Z[,2], 
#                    Label = factor(data$labelsVrais, levels = c(0, 1)))
#   #Création du plot
#   # Create the plot (English version)
#   p <- ggplot(df, aes(x = X1, y = X2, color = Label, alpha = Label)) +
#     geom_point(size = 3) +
#     scale_color_manual(
#       values = c("0" = "blue", "1" = "red"),
#       labels = c("Inlier", "Outlier"),
#       name = "Status"
#     ) +
#     scale_alpha_manual(
#       values = c("0" = 0.1, "1" = 1),  # 0.2 transparency for blue, 1 for red
#       guide = "none"  # Hide alpha from legend
#     ) +
#     labs(
#       title = "Contamination Scenario: Mean and Variance",
#       subtitle = paste("Rate:", rate, "%, k =", k, ", l =", l, ", rho =", rho),
#       x = "X1",
#       y = "X2"
#     ) +
#     theme_minimal() +
#     theme(
#       legend.position = "bottom",
#       plot.title = element_text(face = "bold")
#     )
#   return(p)
# }
# 
