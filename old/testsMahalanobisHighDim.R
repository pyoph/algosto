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


for(i in (1:n)){
  distancesHD[i] = mahalanobis_high_dim(Z[i,],resNaif$meanIter[i,],resNaif$SigmaIter[i,,])
  
  v = diag(resNaif$SigmaIter[i,,])
  
  RIter[i,,] = solve(diag(v,d))%*%resNaif$SigmaIter[i,,]%*%solve(diag(v,d))
  
  c[i]= median(distancesHD[1:i])/d
  
  
  adj_quant = 1 + sum(diag(RIter[i,,]%*%RIter[i,,]))/d^(3/2)
  
  quantNorm = qnorm(0.95)
  
  cutoff_corr[i] = d + quantNorm*(2*adj_quant*(sum(diag(RIter[i,,]%*%RIter[i,,])) -  d^2/n))^(1/2)
  
  if(distancesHD[i] > cutoff_corr[i]/c[i]){outllab[i] = 1}
  
  }


#######Cutoff#####################
