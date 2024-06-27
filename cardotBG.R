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

m

Sigma <- diag(d)

#10000 échantillons tirés selon une loi N(0,Id)

X <- mvtnorm::rmvnorm(n,sigma=Sigma)
dim(X)                      
moyennem = matrix(0,n,d)
dim(moyennem)
((X[4,] - moyennem[3,]) %*% t(X[4,] - moyennem[3,]) - V)/norm((X[4,] - moyennem[3,]) %*% t(X[4,] - moyennem[3,]) - V,type = "F")

dim(X)

V = matrix(0,d,d)

#Stockage des itérations
miter = matrix(0,n,d)
miter[1,] = m
m
miter[1,]

sqrt(sum((m - X[100,])^2))


for (i in 1:999)
{

gamma = c/i^(0.67) #pour un premier essai on prend alpha = 0.67

#Mise à jour de mn

m = m + (gamma*(X[i,] - m))/sqrt(sum(X[i,]  - m)^2)

miter[i,] = m

if(i == 1)
{
  moyennem[i,] = miter[i,]
  #print(moyennem[i,])
}

if(i >= 2){
  moyennem[i,] = colMeans(miter[1:i,])
  
}
  
V = V + gamma*((X[i+1,] - moyennem[i,]) %*% t(X[i+1,] - moyennem[i,]) - V)/norm((X[i+1,] - moyennem[i,]) %*% t(X[i+1,] - moyennem[i,]) - V,type = "F")

moyennem = moyennem*i/(i+1) + 1/(i+1)*m


}

