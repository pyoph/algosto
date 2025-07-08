

#######################################################
#Directories paths à adapter selon votre configuration
#########################################################
paramDir = "C:/Users/Paul/Documents/Simus/ParamsSim/ParamsSim"
simDir =  "C:/Users/Paul/Documents/Simus/DataSim"


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

simNb = 100


# Simulation of datas

for(r in rList){
  for(k in kList){
    for(l in lList){
      for (rho1 in rho1List){  
        for(sim in 1:simNb){
          
          setwd(simDir)
          contParam = ParmsF1(m1, k, l, rho1)
          dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData")
          
          
          if(!file.exists(dataFile)){   
            data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r )
            save(data,file = dataFile)
            print(paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData"," save OK"))
          }
          else(print(paste0('SimData-d', d, '-n', n, '-k', k, '-l', 1, '-rho', rho0,'-r',r , '-sim', sim,".RData"," File exists")))
          print(paste0("r ",r,"number of outliers ",sum(data$labelsVrais == 1)))
          
          
          
        }}}
  }}


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
