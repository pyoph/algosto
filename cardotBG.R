library(RobRegression)
library(Gmedian)
library(ggplot2)


### Article 2015

#Nombre d'échantillons
n = 1000

#Dimension
d = 10 


##Estimation de mn

#Tirage de lois N(0,Id)

m = rep(0,d)

c = 1.5

Sigma <- diag(d)


#Calcul de la vraie médiane géométrique

mvrai <- Gmedian(X)

#Calcul de la vraie MCM par l'algorithme de Weiszfeld

Vvrai <- WeiszfeldCov(X)

#10000 échantillons tirés selon une loi N(0,Id)

X <- mvtnorm::rmvnorm(n,sigma=Sigma)
dim(X)                      
moyennem = rep(0,d)


#Stockage des itérations de matrices dans un tableau
VIter <- array(0, dim = c(d, d, n))

#Initialisation de V

V <- matrix(0,n,d)


dim(moyennem)

dim(X)

V = matrix(0,d,d)

#Stockage des itérations
miter = matrix(0,n,d)



for (i in 1:(n-1))
{

gamma = c/i^(0.8) #pour un premier essai on prend alpha = 0.8

#Mise à jour de mn

m = m + (gamma*(X[i,] - m))/sqrt(sum(X[i,]  - m)^2)

miter[i,] = m

if(i == 1)
{
  moyennem = miter[i,]
}

if(i >= 2){
  moyennem = colMeans(miter[1:i,])
  
}

#Mise à jour de V

V = V + gamma*((X[i+1,] - moyennem) %*% t(X[i+1,] - moyennem) - V)/norm((X[i+1,] - moyennem) %*% t(X[i+1,] - moyennem) - V,type = "F")

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



}


#Stockage des erreurs quadratiques

erreursM <- rep(0,(n-1))

erreursV <- rep(0,(n-1))

#Calcul des erreurs quadratiques

for (i in 1:9999) 
{
  erreursM[i] <- sqrt(sum(mvrai - miter[i])^2)
}



for (i in 1:9999) 
{
  erreursV[i] <- norm(Vvrai$covmedian - VIter[,,i],type = "F")
}

data <- data.frame(
  x = 1000:9999,
  erreursM = erreursM[1000:9999],
  erreursV = erreursV[1000:9999]

)


# Transformer les données en un format long pour ggplot2
data_long <- data.frame(
  x = rep(data$x, 2),
  value = c(data$erreursM, data$erreursV),
  type = rep(c("Erreur M", "Erreur V"), each = nrow(data))
)

# Créer le graphique avec ggplot2
p <- ggplot(data_long, aes(x = x, y = value, color = type)) +
  geom_line(size = 1) +
  scale_x_log10() +   # Échelle logarithmique pour l'axe x
  scale_y_log10() +   # Échelle logarithmique pour l'axe y
  scale_color_manual(values = c("Erreur M" = "red", "Erreur V" = "blue"),
                     labels = c("Erreur M" = "Courbe rouge : Erreur d'estimation de m", 
                                "Erreur V" = "Courbe bleue : Erreur d'estimation de V")) +
  labs(x = "Itérations", y = "Erreur Quadratique", 
       title = "Erreur quadratique des Algorithmes",
       color = "Type d'Erreur") +
  theme_minimal()

# Afficher le graphique
print(p)


###Online learning 

beta = 1/2

#Taille d'un bloc

k = 100

t = n/k

m = rep(0,d)

for (i in 1:t) 
{

S <- rep(0,d)
gamma = c/i^(0.6) 

S <- rowSums(sapply(1:k, function(j) (X[(i - 1) * k + j, ] - m) / sqrt(sum((X[(i - 1) * k + j, ] - m)^2))))

m <- m - gamma/k^(beta)*S/k



}

mvrai = Gmedian(X)
m
mvrai
