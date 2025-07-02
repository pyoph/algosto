r = 25
  
sigmaSq0 <- (1:d); sigmaSq0 <- sigmaSq0 / mean(sigmaSq0)
SigmaContamin <- diag(sqrt(sigmaSq0)) %*% toeplitz(0.995^(0:(d-1))) %*% diag(sqrt(sigmaSq0))

data <- genererEchantillon(n,d,mu1,mu2 = 30*rep(1/sqrt(d), d),p1 = 1- r/100,r/100,Sigma1,Sigma2 = 0.01*SigmaContamin,contamin="moyenne_variance",cluster = "FALSE")

Z = data$Z

resultats = StreamingOutlierDetection(Z,batch = 1)

muhatIter = resultats$miter

SigmaIter = resultats$Sigma

distances = resultats$Sigma

covmed = resultats$moyenneV

P = diag(1/sqrt(colSums(covmed^2)))
  
SigmaF = SigmaIter[nrow(Z),,]
norm(SigmaIter[nrow(Z),,] - Sigma1,"F")

lambda = eigen(SigmaF)$values

if (lambda[1] > 6*lambda[2]){

  print("correction ")
  
  lambdacorr = diag(c(max(lambda[2:d]),lambda[2:d]))
  
P = eigen(covmed)$vectors

Sigmacorr = P %*% lambdacorr %*% t(P)

norm(Sigmacorr - Sigma1,"F")
}

################
#Classification
################

  distances = rep(0, nrow(Z))
  outliers_labels = rep(0,nrow(Z))
  
  for (i in (1:nrow(Z))){
    distances[i] = mahalanobis_generalizedRcpp(Z[i,],muhatIter[i,],eigen(SigmaIter[i,,])$vectors, eigen(SigmaIter[i,,])$values)
    S = distances[i]
    
    cutoff = 1.28*qchisq(.95,df = d)*median(distances[1:i])/qchisq(.5,df = d)    
    if (S > cutoff) {outliers_labels[i] = 1}
    
  }
table(data$labelsVrais,outliers_labels)
  

