SampleCovOnline = function(Z)
{
  nblignes = nrow(Z)
  
  mean = 1.5*rnorm(ncol(Z))
  Sigma = diag(ncol(Z))
  
  meanIter =   matrix(0,nrow(Z),ncol(Z))
  SigmaIter = array(0, dim = c(nrow(Z),ncol(Z),ncol(Z)))
  
  for (i in (1:(nblignes-1)))
  {
    mean   = mean   + (1.0 / (i + 1)) * (Z[i+1,] - mean);
    Sigma = (i + 1)/i*Sigma + 1/(i + 1)*((Z[i+1,] - meanIter[i,])%*%t(Z[i+1,] - meanIter[i,]) - Sigma)
    meanIter[i,] = mean
    SigmaIter[i,,] = Sigma
    

    
      }
  SigmaIter[nrow(Z),,] = Sigma
  return(list(mean = mean, Sigma = Sigma, meanIter = meanIter, SigmaIter = SigmaIter))
}
#Sigma = SampleCovOnline(Z)$Sigma