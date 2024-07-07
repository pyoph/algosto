#install.packages("microbenchmark")


library(RobRegression)
library(Gmedian)
library(ggplot2)
library(far)
library(gridExtra)
library(microbenchmark)
### Article 2015

#Nombre d'itérations

n = 1000

#Variables

#Dimension
d = 10

r <- 0.5

#Taux d'apprentissage
c <- 2:20

mvrai <- rep(0,d)


Y <- mvtnorm::rmvnorm(1000000,sigma=Sigma)

Vvrai <- WeiszfeldCov(Y,nitermax = 1000)$covmedian


estimMV <- function(c,r)
{
  
  
  ##Estimation de mn
  
  #Initialisation de m
  
  
  m <- r*rnorm(d)
  
  
  #m <- m/sqrt(sum(m^2))
  
  #m <- m/sqrt(sum(m^2))
  
  Sigma <- diag(d)
  
  #Initialisation de V
  
  V = matrix(0,d,d)
  
  #Calcul de la vraie médiane géométrique
  
  
  
  #Calcul de la vraie MCM par l'algorithme de Weiszfeld
  
  
  
  moyennem = m
  
  VvectproprX <- eigen(cov(Y))$vectors
  
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(d, d, n))
  
  
  
  
  #Tableau de stockage des vecteurs propres estimés orthonormalisés
  
  Uorth <- array(0,dim = c(n,d,d))
  
  #Stockage des vecteurs propres dans un tableau
  
  
  
  U <- array(1, dim = c(n, d, d))
  
  #Initialisation de U avec des vecteurs propres pas trop éloignés de ceux de cov(X)
  
  for (i in (1:n))
  {
    U[i,,] <- VvectproprX + 0.01*matrix(1,d,d)
  }
  
  
  #Stockage des itérations
  miter = matrix(0,n,d)
  
  
  
  
  
  for (i in 1:(n-1))
  {
    
    
    gamma = c/i^(0.6)
    
    #Test pour gérer la non négativité de Vn
    
    #if (gamma >= norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F"))
    #{
    # gamma = norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F")
    #}
    
    #Mise à jour de mn 
    
    m = m + gamma*((Y[i+1,] - m))/sqrt(sum(Y[i+1,]  - m)^2)
    
    miter[i,] = m
    
    if(i == 1)
    {
      moyennem = miter[i,]
    }
    
    if(i >= 2){
      moyennem = colMeans(miter[1:i,])
      
    }
    
    #Mise à jour de V
    
    
    V = V + gamma*((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V)/norm((Y[i+1,] - moyennem) %*% t(Y[i+1,] - moyennem) - V,type = "F")
    
    moyennem = moyennem*i/(i+1) + 1/(i+1)*m
    
    VIter[,,i] = V
    
    #Calcul de la moyenne de V
    
    if (i == 1){
      moyenneV = V 
    }
    
    if(i >= 2) {
      moyenneV <- apply(VIter[,,1:i], c(1, 2), mean)
      
    }
    
    
    #Mise à jour de moyenneV
    
    moyenneV <- i/(i+1)*moyenneV + 1/(i + 1)*V
    
    #Estimation des vecteurs propres de Vt et orthonormalisation
    
    #for (l in (1:d)) 
    #  {
    # Un =  U[i,l,]/sqrt(U[i,l,]^2)
    
    #   U[i+1,l,] <-  U[i,l,] + 1/(i + 1)*(moyenneV %*% Un - U[i,l,])
    # }
    
    
    #Orthonormalisation des vecteurs propres
    #Uorth[i+1,,] <- orthonormalization(U[i+1,,], basis=TRUE, norm=TRUE)
    #Stockage des erreurs quadratiques
    
    erreursM <- rep(0,(n-1))
    
    erreursV <- rep(0,(n-1))
    
    miter[n,] = miter[n-1,]
    VIter[,,n] = VIter[,,n-1]
    
    
    
    
    #Calcul des erreurs quadratiques
    
    for (i in 1:(n-1)) 
    {
      erreursM[i] <- sqrt(sum(mvrai - miter[i])^2)
    }
    
    
    
    for (i in 1:(n-1)) 
    {
      erreursV[i] <- norm(Vvrai - VIter[,,i],type = "2")
    }
    
  }
  return (list(m=m,V=V,moyennem=moyennem,moyenneV=moyenneV,erreursM = erreursM,erreursV = erreursV,miter = miter,VIter = VIter))
}

data <- data.frame(
  x = 1:(n - 1),
  erreursM = erreursM[1:(n-1)],
  erreursV = erreursV[1:(n-1)]

)


# Transformer les données en un format long pour ggplot2
data_long <- data.frame(
  x = rep(data$x, 2),
  value = c(data$erreursM),
  type = rep(c("Erreur m"), each = nrow(data))
)

graphiquesErreurM <- list()


for (c in 5:8)
{
  m <- estimMV(c,0.5)$m
  erreursM <- estimMV(c,1.5)$erreursM
  
  
  #print(m)
  data <- data.frame(
    x = 1:(n - 1),
    erreursM = erreursM
    
  )

  # Graphique de représentation des erreurs d'estimation de m
  p <- ggplot(data_long, aes(x = x, y = value, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur m" = "red"),
                       labels = c(paste("Erreur d'estimation de m c = ", c))) +
    labs(x = "Itérations", y = "Erreur Quadratique", 
         title = c(paste("Erreur d'estimation de m c = ", c)),
         color = "Type d'Erreur") +
    theme_minimal()
  
  # Afficher le graphique
  print(p)
  graphiquesErreurM <- append(graphiquesErreurM,list(p))  
}


tempsCalcul <- microbenchmark(estimMV(20,1.5),times = 5)

tempsCalcul


for (c in 9:11)
{
  m <- estimMV(c,0.5)$m
  erreursM <- estimMV(c,1.5)$erreursM
  
  
  #print(m)
  data <- data.frame(
    x = 1:(n - 1),
    erreursM = erreursM
    
  )
  
  # Graphique de représentation des erreurs d'estimation de m
  p <- ggplot(data_long, aes(x = x, y = value, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur m" = "red"),
                       labels = c(paste("Erreur d'estimation de m c = ", c))) +
    labs(x = "Itérations", y = "Erreur Quadratique", 
         title = c(paste("Erreur d'estimation de m c = ", c)),
         color = "Type d'Erreur") +
    theme_minimal()
  
  # Afficher le graphique
  print(p)
  graphiquesErreurM <- append(graphiquesErreurM,list(p))  
}

all(sapply(graphiquesErreurM, inherits, what = "gg"))

grid.arrange(grobs = graphiquesErreurM, nrow = 4, ncol = 5)

tempsCalcul <- microbenchmark(estimMV(1,19), times = 100)


p <- autoplot(tempsCalcul) +
  ggtitle("Temps d'exécution de la fonction d'estimation de m et de V") +
  xlab("Fonction") +
  ylab("Temps (microsecondes)") +
  theme_minimal()

# Afficher le graphique
print(p)

###Streaming

beta = 1/2

#Taille d'un bloc

#Foncion pour estimer m et V en fonction de la taille des blocs avec une estimation online

oe <- function(t,k)
{
  

  

  #Initialisation de m
  m = r*rnorm(d)
  
  m <- m/sqrt(sum(m^2))
  
  moyennem = m
  
  V <- matrix(0,d,d)
  
  moyenneV <- V
  
  #Stockage des itérations de médiane géométrique 
  
  miter = matrix(0,t,d)
  
  
  #Stockage des itérations de matrices dans un tableau
  VIter <- array(0, dim = c(d, d, t))
  
  
  
  for (i in 1:t) 
  {
    
    #Somme des médianes géométriques
    S <- rep(0,d)
    
    #Somme matrices 
    
    Smatr <- matrix(0,d,d)
    
    gamma = c/i^(0.6) 
    
    S <- rowSums(sapply(1:k, function(j) (Y[(i - 1) * k + j, ] - m) / sqrt(sum((Y[(i - 1) * k + j, ] - m)^2))))
    
    #Mise à jour de m
    
    m <- m + gamma/k^(beta)*S/k
    
    miter[i,] = m
    
    if(i == 1)
    {
      moyennem = miter[i,]
    }
    
    if(i >= 2){
      moyennem = colMeans(miter[1:i,])
      
    }
    
    moyennem = moyennem*i/(i+1) + 1/(i+1)*m
    
    for (j in 1: k) 
    {
      Smatr <- Smatr +  ((Y[(i - 1) * k + j, ] - m) %*% t(Y[(i - 1) * k + j, ] - m))/norm((Y[(i - 1) * k + j, ] - m) %*% t(Y[(i - 1) * k + j, ] - m),type = "F")-V
      
    }
    
    V <- V + gamma/k*1/k^(beta)*Smatr
    
    
    #Stockage des itérations de matrices dans un tableau
    VIter[,,i] <- V
    
    
    
    #Calcul de la moyenne de V
    
    if (i == 1){
      moyenneV = V 
    }
    
    if(i >= 2) {
      moyenneV <- apply(VIter[,,1:i], c(1, 2), mean)
      
    }
    
    
    #Mise à jour de moyenneV
    
    moyenneV <- i/(i+1)*moyenneV + 1/(i + 1)*V
    
    
    
    #Estimation des vecteurs propres de Vt et orthonormalisation
    
    for (l in (1:d)) 
    {
      Un =  U[i,l,]/sqrt(U[i,l,]^2)
      
      U[i+1,l,] <-  U[i,l,] + 1/(i + 1)*(moyenneV %*% Un - U[i,l,])
    }
    
    
    #Orthonormalisation des vecteurs propres
    Uorth[i+1,,] <- orthonormalization(U[i+1,,], basis=TRUE, norm=TRUE)
    
  }  

  erreursM <- rep(0,(t-1))
  
  erreursV <- rep(0,(t-1))
  
  #Calcul des erreurs quadratiques
  
  for (i in 1:(t-1)) 
  {
    erreursM[i] <- sqrt(sum(mvrai - miter[i])^2)
  }
  
  
  
  for (i in 1:(t-1)) 
  {
    erreursV[i] <- norm(Vvrai - VIter[,,i],type = "F")
  }
  
# print(" fin fonction ")
 
  
  
return (list(m=m,V=V,moyennem=moyennem,moyenneV=moyenneV,erreursM = erreursM,erreursV = erreursV,miter = miter,VIter = VIter))
}




graphiquesErreurM <- list()

for (k in list(2,5,10,100))
{
  t <- n/k
  #print(t)
  resultats <- oe(t,k)
  m <- resultats$m
  erreursM <- resultats$erreursM
  
   
  #print(m)
  data <- data.frame(
    x = 1:(t - 1),
    erreursM = erreursM
    
  )
  
  
  # Transformer les données en un format long pour ggplot2
  data_long <- data.frame(
    x = rep(data$x, 1),
    value = c(data$erreursM),
    type = rep(c("Erreur m"), each = nrow(data))
  )
  # Créer le graphique avec ggplot2
  p <- ggplot(data_long, aes(x = x, y = value, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur m" = "red"),
                       labels = c(paste("Erreur d'estimation de m blocs de taille ", k))) +
    labs(x = "Itérations", y = "Erreur Quadratique", 
         title = paste("Erreur d'estimation de m blocs de taille ", k),
         color = "Type d'Erreur") +
    theme_minimal()
  graphiquesErreurM <- append(graphiquesErreurM,list(p))
}

all(sapply(graphiquesErreurM, inherits, what = "gg"))

grid.arrange(grobs = graphiquesErreurM, nrow = 2, ncol = 2)


graphiqueErreurV <- list()

for (k in list(2,5,10,100))
{
  t <- n/k
  print(k)
  resultats <- oe(t,k)
  
  V <- resultats$V
  
  erreursV <- resultats $erreursV
  
  
  #print(m)
  data <- data.frame(
    x = 1:(t - 1),
    erreursV = erreursV[1:t-1]
    
  )
  
  
  # Transformer les données en un format long pour ggplot2
  data_long <- data.frame(
    x = rep(data$x, 1),
    value = c(data$erreursV),
    type = rep(c("Erreur V"), each = nrow(data))
  )
  # Créer le graphique avec ggplot2
  p <- ggplot(data_long, aes(x = x, y = value, color = type)) +
    geom_line(size = 1) +
    scale_x_log10() +   # Échelle logarithmique pour l'axe x
    scale_y_log10() +   # Échelle logarithmique pour l'axe y
    scale_color_manual(values = c("Erreur V" = "red"),
                       labels = c(paste("Erreur d'estimation de V blocs de taille ", k))) +
    labs(x = "Itérations", y = "Erreur Quadratique", 
         title = paste("Erreur d'estimation de V blocs de taille ", k),
         color = "Type d'Erreur") +
    theme_minimal()
  graphiqueErreurV <- append(graphiqueErreurV,list(p))
  #print(p)
}

all(sapply(graphiqueErreurV, inherits, what = "gg"))

grid.arrange(grobs = graphiqueErreurV, nrow = 2, ncol = 2)


tempsCalcul <- microbenchmark(oe(n/2,2), times = 100)


p <- autoplot(tempsCalcul) +
  ggtitle("Temps d'exécution de la fonction d'estimation de m et de V par blocs") +
  xlab("Fonction") +
  ylab("Temps (microsecondes)") +
  theme_minimal()

# Afficher le graphique
print(p)

