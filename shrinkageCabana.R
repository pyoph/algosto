
#Paramètre dataset Z calcule la comédiane de Z
comedianeSCCM <- function(Z)
{
  
  COM <- matrix(0,ncol(Z),ncol(Z))
  
  for (i in 1:ncol(Z)){
    for (j in 1:ncol(Z))
    {
      COM[i,j] <- median((Z[,i] - median(Z[,i])*(Z[,j] - median(Z[,j]))))
    }
  }
  return(2.198*COM)
}

#Renvoie l'intensité de shrinkage et la target
shrinkageIntensity <- function(Sccm,SigmaVrai = Sigma1)
{
 
  nu <- sum(diag(Sccm))/ncol(Sccm) 
  
  target <- diag(ncol(Sccm))
  
  target <- nu*target
  
  print(target)
  
  etanum <- abs(sum(diag(Sccm^2)) - sum(diag(SigmaVrai^2)))
  
  etadenom <- nu^2*ncol(Sccm)^2 + sum(diag(Sigma1))^2 - 2*nu*sum(diag(SigmaVrai))
  
  eta <- etanum/etadenom
  
 return (list(nu = nu,eta = eta, target = target))
} 
#Renvoie un estimateur de Sigma à partir d'un dataset Z

SigmaShrink <- function(Z)
{
 
  #Calcul de Sccm 
  
  Sccm <- comedianeSCCM(Z)
  
  #Calcul des paramètres du shrinkage
  paramShrink <- shrinkageIntensity(Sccm) 
  
  nu <- paramShrink$nu
  eta <- paramShrink$eta
  target <- paramShrink$target
  
  sigmaHat <- (1 - eta)*Sccm + eta*diag(ncol(Z))
  
  return(sigmaHat)
}