####full method
utils::globalVariables(c("Sigma", "dpar"))


Robust_Variance=function(X,K=ncol(X),par=TRUE,alphaRM=0.75,c='default',w=2,mc_sample_size='default',
                         methodMC='Robbins',niterMC=50,method_MCM='Weiszfeld',
                         eps_vp=10^(-6))
{
  q=ncol(X)
  a=redMCM(X,K,par=par,method_MCM=method_MCM)
  d=a$result$d
  eigenvalues=a$eigenvalues
  if (a$result$d < length(eigenvalues))
  {
    eigenvalues=c(eigenvalues[1:(d-1)],rep(median(eigenvalues[(d):length(eigenvalues)]),q-d+1))
  }
  eigenvalues=eigenvalues[1:(d)]
  eigenvectors=a$eigenvectors
  if (c=='default')
  {
    c=max(eigenvalues)
  }
  if (mc_sample_size=='default')
  {
    mc_sample_size=floor(nrow(X)*q)
    #    mc_sample_size=floor(nrow(X)*q^2/d)
  }
  if (methodMC=='Robbins'){
    eigSigma=RobbinsMC(mc_sample_size=mc_sample_size,vp=eigenvalues[1:(d)],init=eigenvalues,
                       epsilon=10^(-8),alpha=alphaRM,c=c,w=w,samp=mc_sample_size)
  }
  if (methodMC=='Grad'){
    eigSigma=GradMC(mc_sample_size=mc_sample_size,vp=eigenvalues[1:(d)],epsilon=10^(-8),init = eigenvalues,
                    samp=mc_sample_size,niter=niterMC,pas=rep(c,niterMC))
  }
  if (methodMC=='Fix'){
    eigSigma=FixMC(mc_sample_size=mc_sample_size,vp=eigenvalues[1:(d)],epsilon=10^(-8),init = eigenvalues,
                   samp=mc_sample_size,niter=niterMC,pas=rep(c,niterMC))
  }
  #  eigSigma=FixMC(mc_sample_size=1000,niter=50,vp=eigenvalues,epsilon=10^(-8),pas=rep(1,50),samp=50,init=rep(0,length(eigenvalues)))
  eigSigma=eigSigma$vp
  eigSigma=c(eigSigma,rep(eigSigma[d],q-(d)))
  if (length(which(K>=q))>0)
  {
    eigSigma=c(eigSigma,0)
  }
  eigSigma=  apply(cbind(eigSigma,rep(eps_vp,length(eigSigma))),1,FUN=max)
  invsigma=matrix(0 ,  nrow=q,ncol=q)
  if (length(which(K>=q))>0)
  {
    eigSigma=c(eigSigma,0)
    invsigma= (eigenvectors[,1:d])%*% diag(eigSigma[1:d]^(-1)) %*% t(eigenvectors[,1:d])
  }
  if (length(which(K>=q)) ==0)
  {
    invsigma= (eigenvectors[,1:d])%*% diag(eigSigma[1:d]^(-1)-eigSigma[d+1]^(-1)) %*% t(eigenvectors[,1:d])
    diag(invsigma)=diag(invsigma)+eigSigma[d]^(-1)
  }
  Sigmarebuilt=matrix(0 ,  nrow=q,ncol=q)
  Sigmarebuilt=  (eigenvectors[,1:d])%*% diag(eigSigma[1:d]-eigSigma[d+1]) %*% t(eigenvectors[,1:d])
  diag(Sigmarebuilt)=diag(Sigmarebuilt)+eigSigma[d+1]
  return(list(Sigma=Sigmarebuilt,invSigma=invsigma,MCM=a$MCM,
              eigenvalues=eigSigma,eigenvectors=a$eigenvectors,
              MCM_eigenvalues=a$eigenvalues,
              cap=a$cap,reduction_results=a$allresults))
}



Robust_Mahalanobis_regression=function(X,Y,alphaRM=0.66,alphareg=0.66,w=2, lambda=0,
                                       creg='default', ##anciennement sqrt(ncol(X)*ncol(Y))
                                       K=2:30,par=TRUE,epsilon=10^(-8),method_regression='Offline',niter_regression=50,
                                       cRM='default',mc_sample_size='default',method_MCM='Weiszfeld',
                                       methodMC='Robbins',niterMC=50,ridge=1,eps_vp=10^(-4),
                                       nlambda=50,scale='none',tol=10^(-3))

{
  if (scale=='robust')
  {
    Y=DescTools::RobScale(Y,center = FALSE)
    scales=attr(Y,"scaled:scale")
    c=1
    #    svd_results=svd(X)
    #    X=svd_results$u
  }

  if (method_regression=='Online'){
    betachap=RobLassoOnline(X=X,Y=Y,sigid=T,alphalasso=alphareg,w=w, lambda=0,
                            c=creg)
  }
  if (method_regression=='Offline'){
    betachap=RobLassoOffline(X=X,Y=Y,sigid=T, lambda=0,epsilon=epsilon,
                             niter = niter_regression,tol=tol)
  }
  if (is.character(lambda))
  {
    rss <- sqrt(sum((Y - X%*%betachap)^2))
    lambdaMax <- rss / sqrt(sum(betachap^2))
    lambda <- 10^seq(-6, log10(lambdaMax), length.out=nlambda)
  }

  #betachap=RobLassoOnline(X,Y,lambda=0,Sigma = a$Sigma,invSigma = a$invSigma,c=sqrt(p*q),halfinvSigma = a$halfinvSigma)

  Ytilde=Y-X%*%betachap
  a=Robust_Variance(X=Ytilde,K=K,par=par,alphaRM=alphaRM,c=cRM,
                    w=w,mc_sample_size=mc_sample_size,method_MCM=method_MCM,
                    methodMC=methodMC,niterMC=niterMC,eps_vp=eps_vp)
  Residual_Variance=a$Sigma
  if (method_regression=='Online'){
    if (length(lambda) ==1)
    {
      betachap=RobLassoOnline(X=X,Y=Y,lambda=lambda,invSigma = a$invSigma,
                              alphalasso=alphareg,w=w,
                              c=creg,sigid=F,ridge=ridge,
                              init=betachap)
    }
    if (length(lambda)>1){
      if (par==T){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOnline')  %dopar% {
          resultatk=RobLassoOnline(X=X,Y=Y,lambda=dpar,invSigma = a$invSigma,
                                   alphalasso=alphareg,w=w,
                                   c=creg,sigid=F,ridge=ridge,
                                   init=betachap)
          return(resultatk)
        }
        stopCluster(cl)
      }
      if (par==F){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOnline')  %do% {
          resultatk=RobLassoOnline(X=X,Y=Y,lambda=dpar,invSigma = a$invSigma,
                                   alphalasso=alphareg,w=w,
                                   c=creg,sigid=F,ridge=ridge,
                                   init=betachap)
          return(resultatk)
        }
        stopCluster(cl)
      }
    }
  }
  if (method_regression=='Offline'){
    if (length(lambda)==1){
      betachap=RobLassoOffline(X=X,Y=Y,lambda=lambda,invSigma = a$invSigma,
                               sigid=F,ridge=ridge,
                               init=betachap,niter=niter_regression,epsilon=epsilon,tol=tol)
    }
    if (length(lambda)>1){
      if (par==T){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOffline')  %dopar% {
          resultatk=RobLassoOffline(X=X,Y=Y,lambda=dpar,invSigma = a$invSigma,
                                    sigid=F,ridge=ridge,niter=niter_regression,epsilon=epsilon,
                                    init=betachap,tol=tol)
          return(resultatk)
        }
        stopCluster(cl)
      }
      if (par==F){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOffline')  %do% {
          resultatk=RobLassoOffline(X=X,Y=Y,lambda=dpar,invSigma = a$invSigma,
                                    sigid=F,ridge=ridge,niter=niter_regression,epsilon=epsilon,
                                    init=betachap,tol=tol)
          return(resultatk)
        }
        stopCluster(cl)
      }
    }

  }
  if (length(lambda)>1)
  {
    criterion=c()
    for (i in 1:length(lambda))
    {
      criterion=c(criterion,sqrt(sum((Y-X%*%resultat[[i]])^2)))
    }
    critmin=which.min(criterion)
    all_beta=resultat
    lambdaopt=lambda[[critmin]]
    betachap=resultat[[critmin]]
  }
  if (length(lambda)==1)
  {
    all_beta=betachap
    lambdaopt=lambda
    criterion=sqrt(sum((Y-X%*%betachap)^2))
  }
  if (scale=='robust')
  {
    Residual_Variance=diag(scales)%*%Residual_Variance%*%diag(scales)
    all_beta2=all_beta
    #    betachap=solve(t(svd_results$v))%*%diag((svd_results$d)^(-1))%*%betachap%*%diag(scales)
    betachap= betachap%*%diag(scales)
    # all_beta=betachap
    if (length(lambda)>1)
    {
      for (i in 1:length(lambda))
      {
        #      all_beta[[i]]=solve(t(svd_results$v))%*%diag((svd_results$d)^(-1))%*%all_beta2%*%diag(scales)
        all_beta[[i]]= all_beta2[[i]]%*%diag(scales)
      }
    }
  }
  return(list(beta=betachap,Residual_Variance=Residual_Variance,variance_results=a,criterion=criterion,all_beta=all_beta,lambda_opt=lambdaopt))
}

#### Direct regression, with or without Mahalanobis distance



Robust_regression=function(X,Y, Mat_Mahalanobis=diag(rep(1,ncol(Y))),
                           niter=50,lambda=0,c='default',method='Offline',
                           alpha=0.66,w=2,ridge=1,nlambda=50,
                           init=matrix(runif(ncol(X)*ncol(Y))-0.5,nrow=ncol(X),ncol=ncol(Y)),
                           epsilon=10^(-8), Mahalanobis_distance = FALSE,
                           par=TRUE,scale='none', tol=10^(-3))
{
  if (is.character(lambda))
  {
    rss <- sqrt(sum((Y - X%*%betachap)^2))
    lambdaMax <- rss / sqrt(sum(betachap^2))
    lambda <- 10^seq(-6, log10(lambdaMax), length.out=nlambda)
  }
  if (scale == 'robust')
  {
    Y=DescTools::RobScale(Y,center = FALSE)
    scales=attr(Y,"scaled:scale")
    c=1
    #   svd_results=svd(X)
    #    X=svd_results$u
    Mat_Mahalanobis=diag(scales)%*%Mat_Mahalanobis%*%diag(scales)
  }
  if (Mahalanobis_distance==T)
  {
    sigid=FALSE
  }
  if (Mahalanobis_distance==F)
  {
    sigid=T
  }
  if (method=='Online')
  {
    if (length(lambda) ==1){
      betachap=RobLassoOnline(X=X,Y=Y,invSigma=Mat_Mahalanobis,
                              lambda=lambda,c=c,init=init,alphalasso=alpha,w=w,
                              sigid=sigid,ridge=ridge)
    }
    if (length(lambda) >1){
      if (par==T){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOnline')  %dopar% {
          resultatk=RobLassoOnline(X=X,Y=Y,lambda=dpar,invSigma=Mat_Mahalanobis,
                                   alphalasso=alpha,w=w,
                                   c=c,sigid=sigid,ridge=ridge)
          return(resultatk)
        }
        stopCluster(cl)
      }
      if (par==F){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOnline')  %do% {
          resultatk=RobLassoOnline(X=X,Y=Y,lambda=dpar,invSigma=Mat_Mahalanobis,
                                   alphalasso=alpha,w=w,
                                   c=c,sigid=sigid,ridge=ridge)
          return(resultatk)
        }
        stopCluster(cl)
      }

    }

  }
  if (method=='Offline')
  {
    if (length(lambda)==1){
      betachap=RobLassoOffline(X=X,Y=Y,invSigma=Mat_Mahalanobis,
                               niter=niter,ridge=ridge,lambda=lambda,
                               init=init,sigid=sigid,epsilon=epsilon,tol=tol)
    }
    if (length(lambda)>1){
      if (par==T){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,
                                  .export = c('RobLassoOffline'),.packages=c('fastmatrix'))  %dopar% {
                                    resultatk=RobLassoOffline(X=X,Y=Y,lambda=dpar,invSigma=Mat_Mahalanobis,
                                                              sigid=sigid,ridge=ridge,niter=niter,
                                                              epsilon=epsilon,tol=tol)
                                    return(resultatk)
                                  }
        stopCluster(cl)
      }
      if (par==F){
        numCores = min(parallel::detectCores()-2,length(lambda))
        cl = parallel::makeCluster(numCores  )
        doParallel::registerDoParallel(cl)
        resultat=foreach::foreach(dpar=lambda,.multicombine = TRUE, .inorder = T,.export = 'RobLassoOffline')  %do% {
          resultatk=RobLassoOffline(X=X,Y=Y,lambda=dpar,invSigma=Mat_Mahalanobis,
                                    sigid=sigid,ridge=ridge,niter=niter,
                                    epsilon=epsilon,tol=tol)
          return(resultatk)
        }
        stopCluster(cl)
      }
    }
  }
  if (length(lambda)>1)
  {
    criterion=c()
    for (i in 1:length(lambda))
    {
      criterion=c(criterion,sqrt(sum((Y-X%*%resultat[[i]])^2)))
    }
    critmin=which.min(criterion)
    all_beta=resultat
    lambdaopt=lambda[[critmin]]
    betachap=resultat[[critmin]]
  }
  if (length(lambda)==1)
  {
    all_beta=betachap
    lambdaopt=lambda
    criterion=sqrt(sum((Y-X%*%betachap)^2))
  }
  if (scale == 'robust')
  {
    all_beta2=all_beta
    #    betachap=solve(t(svd_results$v))%*%diag((svd_results$d)^(-1))%*%betachap%*%diag(scales)
    betachap= betachap%*%diag(scales)
    # all_beta=betachap
    if (length(lambda)>1)
    {
      for (i in 1:length(lambda)){
        #        all_beta[[i]]=solve(t(svd_results$v))%*%diag((svd_results$d)^(-1))%*%all_beta2%*%diag(scales)
        all_beta[[i]]= all_beta2[[i]]%*%diag(scales)
      }
    }
  }

  return(list(beta=betachap,criterion=criterion,all_beta=all_beta,lambda_opt=lambdaopt))
}
