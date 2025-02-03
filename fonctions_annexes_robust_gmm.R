###### Reconstruction des valeurs propres a partir de la methode de gradient iterative
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour la descente de gradient
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations  pour differentes iteration
#### pas =pas de la descente de gradient


GradMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp)
{
  p=length(vp)
  vp2=init
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  vplist=c()
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=Y[i,]
      E1= E1 +  Z^2*(sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
      E2= E2 +  (sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
    }
    vp0=vp2
    vp2=vp2 - pas[l]*mc_sample_size^(-1)*vp2*E1 + pas[l]*mc_sample_size^(-1) * (vp)*E2
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==l)) >0)
    {
      vplist=cbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=l,vplist=vplist))
}


###### Reconstruction des valeurs propres a partir de la methode de point fixe
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour l'algo
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations de point fixe pour differentes iteration


FixMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp)
{
  p=length(vp)
  vplist=c()
  vp2=init
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=Y[i,]
      E1= E1 +  Z^2*(sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
      E2= E2 +  (sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
    }
    vp0=vp2
    vp2 = E2* (vp)*(E1)^(-1)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==l)) >0)
    {
      vplist=cbind(vplist,vp2)
    }
    #    if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=l))
}


###### Reconstruction des valeurs propres a partir de la methode de gradient recursive
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations pour differentes tailles d'?chantillon



RobbinsMC=function(mc_sample_size=1000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=vp)
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
      lambdalist=cbind(lambdalist,lambda)
      vplist=cbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=i, lambdalist=lambdalist, vplist=vplist))
}




Weiszfeld_init <- function (X,init=rep(0,ncol(X)), weights = NULL, epsilon=1e-08, nitermax = 100)
{
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1,nrow(X))
  
  return(Weiszfeld_init_rcpp(X,init=init, weights= weights, epsilon=epsilon, nitermax=nitermax))
}

###fonction R pour estimer la matrice de covariance mediane via Weiszfeld en reglant l'initialisation


WeiszfeldCov_init <- function(X, init=rep(0,ncol(X)),init_cov=(diag(ncol(X))), weights=NULL, scores=2, epsilon=1e-08, nitermax = 100){
  ### Computation of the Geometric covariation matrix
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)
  ### require library(rARPACK)
  X <- as.matrix(X)
  n <- nrow(X)
  if (is.null(weights)){ weights <- rep(1,n)}
  Wmed.est <- Weiszfeld_init_rcpp(X,init=init,weights = weights,epsilon=epsilon,nitermax=nitermax)
  WMCM.est <- WeiszfeldCovMat_init_rcpp(X,init_cov=init_cov, median_est = Wmed.est$median,weights=weights,epsilon=epsilon,nitermax=nitermax)
  if (scores==FALSE){
    return(list(median = Wmed.est$median, covmedian=WMCM.est$median, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
  else {
    ### Computation of the eigenvectors and scores
    vectors <- RSpectra::eigs_sym(WMCM.est$median, scores)$vectors
    vscores = sweep(X,2,Wmed.est$median)%*%vectors
    return(list(median=Wmed.est$median, covmedian=WMCM.est$median, scores=vscores, vectors=vectors, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
}

###fonction RCPP pour estimer la median via la descente de gradient en reglant l'initialisation





Gmedian_init <- function (X, init = rep(0,ncol(X)), gamma = 2, alpha = 0.75, nstart=1, epsilon=1e-08,weights=NULL)
{
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1,nrow(X))
  med.X = Gmedianrowvec_init_rcpp(X,init=init,weights=weights, gamma=gamma, alpha=alpha, nstart=nstart, epsilon=epsilon)
  return(med.X)
}


GmedianCov_init <- function(X, init=rep(0,ncol(X)),init_cov=diag(ncol(X)),weights=NULL, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1){
  ### Computation of the Geometric covariation matrix
  ### with averaged stochastic gradient algorithms
  ### input : X   (n x p matrix, n observations, in dimension p)
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)
  ### require library(rARPACK)
  if (is.null(weights)) weights <- rep(1,nrow(X))
  Gmed.est = Gmedian_init(X,init=init,weights=weights,gamma=gamma,alpha=alpha,nstart=nstart)
  GMCM.est = MedianCovMatRow_init_rcpp(X,init_cov=init_cov,Gmedian=Gmed.est,weights=weights,gamma=gc,alpha=alpha,nstart=nstart)
  if (scores==FALSE){
    return(list(median = Gmed.est,covmedian=GMCM.est))
  }
  else {
    ### Computation of the eigenvectors and scores
    vectors <- RSpectra::eigs_sym(GMCM.est, scores)$vectors
    scores = sweep(X,2,Gmed.est)%*%vectors
    return(list(median=Gmed.est,covmedian=GMCM.est,scores=scores,vectors=vectors))
  }
}



MCMrebuild=function(MCM,d,eigenvalues,eigenvectors)
{
  q=ncol(MCM)
  eigenvectors=eigenvectors[,1:d]
  eigenvaluesrebuilt=eigenvalues
  eigenvaluesrebuilt[(d+1):length(eigenvalues)] = rep(median(eigenvalues[(d+1):length(eigenvalues)]),length(eigenvalues) - d)
  
  MCMrebuilt=matrix(0 ,  nrow=q,ncol=q)
  MCMrebuilt=  (eigenvectors[,1:d])%*% diag(eigenvaluesrebuilt[1:d]-eigenvaluesrebuilt[d+1]) %*% t(eigenvectors[,1:d])
  diag(MCMrebuilt)=diag(MCMrebuilt)+eigenvalues[d+1]
  return(list(MCMrebuilt=MCMrebuilt ,error=sum((MCMrebuilt- MCM)^2),d=d))
}


redMCM=function(X,K=1:20,par=T,method_MCM='Weiszfeld'){
  cap='The length must be larger than 10'
  cutoff='The cutoff method is only used if the length of K is smaller than 10'
  ####Attention il faudra utiliser WeiszfeldCovInit
  if (method_MCM=='Weiszfeld')
  {
    MCM=WeiszfeldCov_init(X,scores =0)
  }
  if (method_MCM=='ASGD')
  {
    MCM=GmedianCov_init(X,scores =0)
  }
  MCM=MCM$covmedian
  q=ncol(MCM)
  if (length(which(K >= q)) >0)
  {
    resultat=list(MCMrebuilt=MCM,d=q)
    bestresult=list(MCMrebuilt=MCM,d=q)
    eig=eigen(MCM)
    eigenvalues=eig$values
    eigenvectors=eig$vectors
  }
  if  (length(which(K >= q)) == 0)
  {
    #    eig <- RSpectra::eigs_sym(MCM, max(K)+1+floor(sqrt(max(0,q-max(K)^2))))
    eig <- RSpectra::eigs_sym(MCM, max(K)+1)
    eigenvalues =eig$values
    eigenvectors=eig$vectors
    if (length(K)>1){
      if (par==T){
        numCores = min(parallel::detectCores()-2,length(K))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=K,.multicombine = TRUE, .inorder = T,.export = 'MCMrebuild',
                                  .packages = c( "mvtnorm","foreach","Rcpp","RcppArmadillo"),.combine='list')  %dopar% {
                                    resultatk=MCMrebuild(MCM=MCM,d=dpar,eigenvalues=eigenvalues,eigenvectors=eigenvectors[,1:dpar])
                                    return(resultatk)
                                  }
        stopCluster(cl)
      }
      if (par==F){
        resultat=foreach::foreach(dpar=K,.multicombine = TRUE, .inorder = T,.export = 'MCMrebuild',
                                  .packages = c( "mvtnorm","foreach","Rcpp","RcppArmadillo"),.combine='list')  %do% {
                                    resultatk=MCMrebuild(MCM=MCM,d=dpar,eigenvalues=eigenvalues,eigenvectors=eigenvectors[,1:dpar])
                                    return(resultatk)
                                  }
      }
      err=c()
      for (i in 1:length(K))
      {
        err=c(err,resultat[[i]]$error)
      }
    }
    if (length(K)==1)
    {
      bestresult=MCMrebuild(MCM=MCM,d=K,eigenvalues=eigenvalues,eigenvectors=eigenvectors[,1:K])
      resultat=bestresult
    }
    
    if (length(K)>1)
    {
      if (length(K)<10)
      {
        cutoff=KneeArrower::findCutoff(K, err,method="curvature")
        Kopt=min(max(K),floor(cutoff$x)+1)
        I=which(K==Kopt)
        bestresult=resultat[[I]]
      }
    }
    if (length(K) >=10 )
    {
      dim=(K) *(1+q) - (K)*(K+1)/2+1
      capu=cbind(K,sqrt(dim),dim,err)
      cap=capushe(capu)
      Kopt=as.numeric(cap@DDSE@model)
      I=which(K==Kopt)
      bestresult=resultat[[I]]
    }
  }
  return(list(result=bestresult,allresults=resultat,cutoff=cutoff,cap=cap,MCM=MCM,eigenvalues=eigenvalues,eigenvectors=eigenvectors))
}





RobLassoOnline=function(X,Y,invSigma=diag(rep(1,ncol(Y))),
                        lambda=0,c='default',init=matrix(runif(ncol(X)*ncol(Y))-0.5,nrow=ncol(X),ncol=ncol(Y)),
                        sigid=T,alphalasso=0.66,w=2,ridge=1)
{
  betachap=init
  betabar=init
  taun=0
  #  p=ncol(X)
  #  H=diag(rep(1,p))
  if (sigid==T)
  {
    if (c=='default')
    {
      c=(ncol(Y)*ncol(X))
    }
    for (j in 1:nrow(Y))
    {
      a= (X[j,])%*%((Y[j,]- X[j,]%*%betachap))
      b=  sum((Y[j,]- X[j,]%*%betachap)^2)
      betachap=betachap+ c*j^(-0.66)*( a/sqrt(b) - lambda*betachap/sqrt(sum(betachap^2)))
      taun=taun+log(j+1)^w
      betabar=betabar+log(j+1)^w/taun*(betachap - betabar)
      #      grad=grad- a + lambda*beta/sqrt(sum(beta^2))
    }
  }
  if (sigid ==F){
    if (c=='default')
    {
      c=((1/ncol(Y))*sum(invSigma^2))^(-1/2)
    }
    for (j in 1:nrow(Y))
    {
      #      H=H-as.numeric((1+X[j,]%*%(H%*%(X[j,])))^(-1))*(H%*%X[j,])%*%(t(X[j,])%*%H)
      a= (X[j,])%*%((Y[j,]- X[j,]%*%betachap)%*%invSigma)
      b=  as.numeric((Y[j,]- X[j,]%*%betachap)%*%(invSigma%*%t(Y[j,]- X[j,]%*%betachap)))
      betachap=betachap+ c*j^(-alphalasso)*( a/sqrt(b) - lambda*betachap*(sqrt(sum(betachap^2))^(ridge-2)))
      taun=taun+log(j+1)^w
      betabar=betabar+log(j+1)^w/taun*(betachap - betabar)
      #      grad=grad- a + lambda*beta/sqrt(sum(beta^2))
    }
  }
  return(betabar)
}



RobLassoOffline=function(X,Y,invSigma=diag(rep(1,ncol(Y))),niter=50,ridge=1,lambda=0,
                         init=matrix(runif(ncol(X)*ncol(Y))-0.5,nrow=ncol(X),ncol=ncol(Y)),sigid=T,epsilon=10^(-8),tol=10^(-3))
{
  betachap=init
  #  p=ncol(X)
  #  H=diag(rep(1,p))
  if (lambda==0){
    if (sigid==T)
    {
      for (i in 1:niter)
      {
        mat=diag(rep(1,ncol(X)))
        #matinv=diag(rep(1,ncol(X)))
        grad=0
        
        
        for (j in 1:nrow(X))
        {
          norminv= as.numeric(sum((Y[j,] - X[j,]%*%betachap)^2))
          if (norminv > tol){
            norminv=1/sqrt(norminv)
          } else {norminv=0}
          grad=grad + (X[j,])%*%t(Y[j,])*norminv
          mat=mat+X[j,]%*%t(X[j,])*norminv
          #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
        }
        grad=grad/nrow(X)
        mat=mat/nrow(X)
        betaold=betachap
        betachap= solve(mat)%*%grad
        if (sum((betachap- betaold)^2)< epsilon)
        {
          break
        }
      }
      #    for (j in 1:nrow(X))
      #    {
      #      norminv= (sqrt(sum(Y[j,] - X[j,]%*%betachap)^2))
      #      if (norminv > tol){
      #        norminv=as.numeric(norminv^(-1))
      #      } else {norminv=0}
      #      grad=grad + (X[j,])%*%(Y[j,] - X[j,]%*%betachap)*norminv
      #      mat=mat+X[j,]%*%t(X[j,])*norminv
      #      #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
      #    }
      #    betachap=betachap + solve(mat)%*%grad
      #  }
    }
    if (sigid==F)
    {
      for (i in 1:niter)
      {
        mat=diag(rep(1,ncol(X)))
        #matinv=diag(rep(1,ncol(X)))
        grad=0
        
        
        for (j in 1:nrow(X))
        {
          norminv=stats::mahalanobis(x = Y[j,] ,center = X[j,]%*%betachap,cov = invSigma,inverted = TRUE)
          #          norminv= as.numeric( (Y[j,] - X[j,]%*%betachap)%*%invSigma%*%t(Y[j,] - X[j,]%*%betachap))
          if (norminv > tol){
            norminv=1/sqrt(norminv)
          } else {norminv=0}
          grad=grad + (X[j,])%*%t(Y[j,])*norminv
          mat=mat+X[j,]%*%t(X[j,])*norminv
          #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
        }
        grad=grad/nrow(X)
        mat=mat/nrow(X)
        betaold=betachap
        betachap= solve(mat)%*%grad
        if (sum((betachap- betaold)^2)< epsilon)
        {
          break
        }
      }
      #    for (j in 1:nrow(X))
      #    {
      #      norminv= (sqrt(sum(Y[j,] - X[j,]%*%betachap)^2))
      #       if (norminv > tol){
      #        norminv=as.numeric(norminv^(-1))
      #       } else {norminv=0}
      #      grad=grad + (X[j,])%*%(Y[j,] - X[j,]%*%betachap)*norminv
      #      mat=mat+X[j,]%*%t(X[j,])*norminv
      #      #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
      #    }
      #    betachap=betachap + solve(mat)%*%grad
      #  }
    }
  }
  if (lambda >0){
    betalin=c(betachap)
    if (sigid==T)
    {
      for (i in 1:niter)
      {
        mat=diag(rep(1,ncol(X)))
        #matinv=diag(rep(1,ncol(X)))
        grad=0
        
        
        for (j in 1:nrow(X))
        {
          norminv= as.numeric(sum((Y[j,] - X[j,]%*%betachap)^2))
          if (norminv > tol){
            norminv=1/sqrt(norminv)
          } else {norminv=0}
          grad=grad + (X[j,])%*%((Y[j,] - X[j,]%*%betachap))*norminv
          #                    grad=grad + (X[j,])%*%((Y[j,]))*norminv
          mat=mat+X[j,]%*%t(X[j,])*norminv
          #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
        }
        if (sum(betachap^2)> 0){
          grad=c(grad/nrow(X) - lambda*betachap/as.numeric(sqrt(sum(betachap^2))) )
        }
        if (sum(betachap^2)==0)
        {
          grad=c(grad)/nrow(X)
        }
        #                grad=c(grad/nrow(X))
        mat=mat/nrow(X)
        if (sum(betachap^2)>0){
          #          hess=fastmatrix::kronecker.prod(mat,diag(ncol(Y)))+ lambda * diag(ncol(X)*ncol(Y)) / as.numeric(sqrt(sum(betachap^2)))
          hess=fastmatrix::kronecker.prod(diag(ncol(Y)),mat)+ lambda * diag(ncol(X)*ncol(Y)) / as.numeric(sqrt(sum(betachap^2)))
          hess=solve(hess)
        }
        if (sum(betachap^2)==0)
        {
          #          hess=fastmatrix::kronecker.prod(solve(mat),diag(ncol(Y)))
          hess=fastmatrix::kronecker.prod(diag(ncol(Y)),solve(mat))
        }
        betaold=betachap
        #betalin= betalin+  (hess)%*%grad
        betalin=betalin+grad%*%hess
        betachap= matrix(betalin,ncol=ncol(Y),byrow=F)
        if (sum((betachap- betaold)^2)< epsilon)
        {
          break
        }
      }
      #    for (j in 1:nrow(X))
      #    {
      #      norminv= (sqrt(sum(Y[j,] - X[j,]%*%betachap)^2))
      #      if (norminv > tol){
      #        norminv=as.numeric(norminv^(-1))
      #       } else {norminv=0}
      #      grad=grad + (X[j,])%*%(Y[j,] - X[j,]%*%betachap)*norminv
      #      mat=mat+X[j,]%*%t(X[j,])*norminv
      #      #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
      #    }
      #    betachap=betachap + solve(mat)%*%grad
      #  }
    }
    if (sigid==F)
    {
      for (i in 1:niter)
      {
        mat=diag(rep(1,ncol(X)))
        #matinv=diag(rep(1,ncol(X)))
        grad=0
        
        
        for (j in 1:nrow(X))
        {
          norminv=stats::mahalanobis(x = Y[j,] ,center = X[j,]%*%betachap,cov = invSigma,inverted = TRUE)
          #          norminv= as.numeric( (Y[j,] - X[j,]%*%betachap)%*%invSigma%*%t(Y[j,] - X[j,]%*%betachap))
          if (norminv > tol){
            norminv=1/sqrt(norminv)
          } else {norminv=0}
          grad=grad + (X[j,])%*%((Y[j,] - X[j,]%*%betachap)%*%invSigma)*norminv
          #                    grad=grad + (X[j,])%*%((Y[j,])%*%invSigma)*norminv
          mat=mat+X[j,]%*%t(X[j,])*norminv
          #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
        }
        if (sum(betachap^2)> 0){
          grad=c(grad/nrow(X) - lambda*betachap/as.numeric(sqrt(sum(betachap^2))) )
        }
        if (sum(betachap^2)==0)
        {
          grad=c(grad)/nrow(X)
        }
        #        grad=c(grad/nrow(X))
        mat=mat/nrow(X)
        if (sum (betachap^2)>0){
          #       hess=fastmatrix::kronecker.prod(mat,invSigma)+ lambda * diag(ncol(X)*ncol(Y)) / as.numeric(sqrt(sum(betachap^2)))
          hess=fastmatrix::kronecker.prod(invSigma,mat)+ lambda * diag(ncol(X)*ncol(Y)) / as.numeric(sqrt(sum(betachap^2)))
          hess=solve(hess)
        }
        if (sum(betachap^2)==0)
        {
          #          hess=fastmatrix::kronecker.prod((mat),invSigma)
          hess=fastmatrix::kronecker.prod(invSigma,(mat))
          hess=solve(hess)
        }
        
        betaold=betachap
        #        betalin= betalin+ hess%*%grad
        betalin= betalin+ grad%*%hess
        betachap= matrix(betalin,ncol=ncol(Y),byrow=F)
        if (sum((betachap- betaold)^2)< epsilon)
        {
          break
        }
      }
      #    for (j in 1:nrow(X))
      #    {
      #      norminv= (sqrt(sum(Y[j,] - X[j,]%*%betachap)^2))
      #      if (norminv > tol){
      #        norminv=as.numeric(norminv^(-1))
      #      } else {norminv=0}
      #      grad=grad + (X[j,])%*%(Y[j,] - X[j,]%*%betachap)*norminv
      #      mat=mat+X[j,]%*%t(X[j,])*norminv
      #      #matinv= matinv - norminv * as.numeric((1+  norminv * as.numeric(X[j,]%*%matinv%*%(X[j,]))))^(-1)*(matinv%*%X[j,])%*%(t(X[j,])%*%matinv)
      #    }
      #    betachap=betachap + solve(mat)%*%grad
      #  }
    }
  }
  
  return( betachap )
}
