library(RobRegression)
library(Gmedian)
library(ggplot2)


### Article 2015

#Nombre d'échantillons
n = 100

#Dimension
d = 10 


##Estimation de mn

#Tirage de lois N(0,Id)

m = rep(0,d)

c = 1.5

U = rnorm(d)

sommem = 0

#Calcul de mn et de bar(mn)

for (i in 1:n)
{
U = rnorm(d)

#Calcul de la moyenne des m jusqu'à la i-1 ième itération
moyennem = mean(m)

gamma = c/i^(3/2) #pour un premier essai on prend alpha = 3/2

#Mise à jour de mn

m = m + (gamma*(U - m))/sqrt((U  - m)^2)

#Mise à jour de la moyenne des mn

moyennem = moyennem - 1/(i + 1)*(moyennem - m)

}
