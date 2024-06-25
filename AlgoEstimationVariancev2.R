###Code pour simuler le deuxième algorithme d'estimation des valeurs propres de la matrice de covariance et comparer les erreurs.

#install.packages("RobRegression")
#install.packages("Gmedian")

library(RobRegression)
library(Gmedian)
library(ggplot2)

#Dimension

d = 5

#Coefficient beta 

beta = 0.3

#Paramètre tau

tau = 4


#Calcul de la fonction h

hh <- function(delta,lambda,U)
{
  
  somme  <- sum((lambda*U - delta)^2) + sum((lambda * U)%*%t((lambda * U))) - sum(((lambda * U)^2))
  
  
  somme <- somme^(-1/2)
  
  return(somme)
}

# Fonction pour créer une matrice de covariance définie positive de taille n x n
create_positive_definite_matrix <- function(n, min_val, max_val) {
  L <- matrix(runif(n^2, min = min_val, max = max_val), nrow = n, ncol = n)
  L[upper.tri(L)] <- 0
  A <- L %*% t(L)
  return(A)
}

# Créer une matrice de covariance 10x10 avec des valeurs entre -10 et 10
A <- create_positive_definite_matrix(10, -0.2,0.2)

# Vérifier que A est définie positive
eigenvalues <- eigen(A)$values
print(eigenvalues)

# Afficher la matrice de covariance
print(A)


print(eigen(A))

#A
#Sigma <- diag(10)

Sigma <- diag(d)

#10000 échantillons tirés selon une loi N(0,Id)

X <- mvtnorm::rmvnorm(n=10000,sigma=Sigma)

dim(X)


t(X) %*% X

#Estimation de la MCM de X selon l'algorithme de Weiszfeld

Xmcm <- WeiszfeldCov(X)

#Valeurs propres de la MCM detoile

detoile <- t(eigen(Xmcm$covmedian)$values)

detoile

#On démarre l'algorithme pas trop loin des valeurs propres de Sigma

lambdaestim <- detoile

lambdaestim

dim(lambdaestim)
n=10000

Un=1*rep(1,d)
Vn=0
lambdaestimlist=c()
for (i in (1:n))
{
  
  #Mise à jour de Un
  U=rnorm(d)
  h=hh(detoile,lambdaestim,U^2)
  Un <- Un + (log(i+1)^tau)*U^2*h
  invUn <- (Un/i)^(-1)
  Vn <- Vn + (log(i+1)^tau)*h
  lambdaestim <- invUn * detoile * Vn/i
  lambdaestimlist=rbind(lambdaestimlist,lambdaestim)
}


lambdaestim

diff <- rep(1,d) - lambdaestim

norm(diff^2)


RobbinsMC=function(mc_sample_size=10000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=detoile)
{
  p=length(vp)
  vp2=init
  lambda=init
  lambdalist=c()
  vplist=c()
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
  slog=1
  for (i in 1:mc_sample_size)
  {
    Z=Y[i,]
    # Z2=X2[i,]
    E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    vp0=vp2
    E2=    (sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    lambda =lambda  - c*i^(-alpha)*lambda*E1 + c*i^(-alpha)* (vp)*E2
    slog=slog+log(i+1)^w
    vp2=vp2+log(i+1)^w *((slog)^(-1)) *(lambda - vp2)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==i) > 0))
    {
      lambdalist=rbind(lambdalist,lambda)
      vplist=rbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=i, lambdalist=lambdalist, vplist=vplist))
}


estimRM= RobbinsMC(n,vp=detoile,samp=1:n)

#Calcul des erreurs pour l'algorithme de Robbins Monro
lambdavrai=rep(1,d)
erreur=estimRM$lambdalist- matrix(rep(lambdavrai,n),nrow=n,byrow=T)
err=rowSums(erreur^2)

#Calcul des erreurs avec l'algorithme se basant sur une inversion de phi
erreur2=lambdaestimlist- matrix(rep(lambdavrai,n),nrow=n,byrow=T)
err2=rowSums(erreur2^2)


plot(1:n,err2,'l',log='xy')
lines(1:n,err,'l',col='red',log='xy')

data <- data.frame(
  x = 1:n,
  err2 = err2,
  err = err
)


# Transformer les données en un format long pour ggplot2
data_long <- data.frame(
  x = rep(data$x, 2),
  value = c(data$err2, data$err),
  algorithm = rep(c("Algorithme 1", "Algorithme 2"), each = n)
)

# Créer le graphique
p <- ggplot(data_long, aes(x = x, y = value, color = algorithm)) +
  geom_line() +        
  scale_x_log10() +                                  # Échelle logarithmique pour l'axe x
  scale_y_log10() +                                  # Échelle logarithmique pour l'axe y
  scale_color_manual(values = c("Algorithme 1" = "red", "Algorithme 2" = "blue"),
                     labels = c("Algorithme 1" = "Courbe rouge : Algorithme avec phi inverse", 
                                "Algorithme 2" = "Courbe bleue : Algorithme Robins Monro")) +
  labs(x = "Itérations", y = "Erreur quadratique", title = "Erreur avec les deux algorithmes : échelle logarithmique") +
  theme_minimal()

# Afficher le graphique
print(p)

err[49500:50000]

# Sélectionner les dernières 1000 valeurs de chaque vecteur
errRM <- err[(length(err) - 999):length(err)]
erreurPhiInv <- err2[(length(err2) - 999):length(err2)]


# Créer un data frame pour ggplot2
data <- data.frame(
  value = c(errRM, erreurPhiInv),
  type = rep(c("Algorithme Robbins Monro", "Algorithme phi inverse"), each = 1000)
)

# Créer le boxplot avec les couleurs spécifiées
p <- ggplot(data, aes(x = type, y = value, fill = type)) +
  geom_boxplot(alpha = 0.5) +  # Opacité réduite pour voir les points individuels
  #scale_y_log10() +            # Échelle logarithmique pour l'axe y
  labs(x = "Type d'erreur", y = "Valeur", title = "Boxplot des erreurs") +
  theme_minimal() +
  theme(legend.position = "none")  

# Définir les couleurs des boxplots
p + scale_fill_manual(values = c("Algorithme Robbins Monro" = "blue", "Algorithme phi inverse" = "red"))

print(p)
