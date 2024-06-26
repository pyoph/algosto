library(RobRegression)
library(Gmedian)
library(ggplot2)


### Article 2015

#Nombre d'échantillons
n = 100000

#Dimension
d = 10 


##Estimation de mn

#Tirage de lois N(0,Id)

m = rep(0,d)

c = 1.5


Sigma <- diag(d)

#10000 échantillons tirés selon une loi N(0,Id)

X <- mvtnorm::rmvnorm(n,sigma=Sigma)
dim(X)                      
moyennem = matrix(0,n,d)

dim(X)

#Stockage des itérations
miter = matrix(0,n,d)

#Calcul de mn et de bar(mn)

for (i in 1:n)
{

gamma = c/i^(0.67) #pour un premier essai on prend alpha = 0.67

#Mise à jour de mn

m = m + (gamma*(X[i] - m))/sqrt((X[i]  - m)^2)

miter[i,] = m

moyennem[i,] = moyennem[i,]*i/(i+1) + 1/(i+1)*m

}

miter[100,]
mvrai = Gmedian(X)
sqrt(sum(mvrai^2))
sqrt(sum(m^2)) 

moyennem