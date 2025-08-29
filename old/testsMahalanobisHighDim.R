mahalanobis_high_dim = function(X,mu,Sigma){
  S = 0
  for(i in (1:length(mu))){
    if(Sigma[i,i] != 0){
    S = S + (X[i] - mu[i])^2/(Sigma[i,i])}
  }
  
  return(S)
  
}

distancesHD = rep(0,n)
RIter = array(0,dim= c(n,d,d))
c = rep(0,n)
cutoff_corr = rep(0,n)
outllab = rep(0,n)
v = diag(Sigma0)
R = solve(diag(v,d))%*%Sigma0%*%solve(diag(v,d))
adj_quant = 1 + sum(diag(R%*%R))/d^(3/2)

k = 0;l=1;rho1= 0.3

contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)

Z = data$Z

for(i in (1:n)){
  distancesHD[i] = mahalanobis_high_dim(Z[i,],rep(0,d),Sigma0)
  
  v = diag(Sigma0)
  
  #RIter[i,,] = solve(diag(v,d))%*%Sigma0%*%solve(diag(v,d))
  
  c[i]= median(distancesHD[1:i])/d
  
  
  quantNorm = qnorm(0.975)
  
  cutoff_corr[i] = d + quantNorm*(2*adj_quant*(sum(diag(R%*%R)) -  d^2/n))^(1/2)
  
  if(distancesHD[i] > cutoff_corr[i]){outllab[i] = 1}
  
  }

fp = sum(data$labelsVrais == 0 & outllab == 1)

print(paste0("false positives = ",fp))
