#######################Génération des données#############################
for (sim in 1:1e2){
  for(sc in scenarios){
    
    k = sc$k
    l = sc$l
    rho1 = sc$rho1
    
    
    id_pool <- sample(1:n)
    
    outlier_sets <- vector("list", length(rList))
    
    id_pool <- sample(1:n)
    
    r_max <- max(rList)
    
    outlier_sets = list()
    
    for (m in seq_along(rList)) {
      n_active = floor(rList[m]/100*n)
      outlier_sets[[m]] = id_pool[1:n_active]
    }
    
    for (m in seq_along(rList)){
      
      r = rList[m]
      dataFile <- paste0('SimData-d', d, '-n', n, '-k', k, '-l', l, '-rho', rho1,'-r',r ,"-sim",sim,".RData")
      
      print(dataFile)
      #contParam = ParmsF1(m1, k, l, rho1)
      setwd(SimDir)
        contParam = ParmsF1(m1, k, l, rho1)
      ok = FALSE
      if(r == 0){
        
        
        if(!file.exists(dataFile)){
        data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r)
        
        save(data,file = dataFile)
          }
        else {
          print("File exists")
          load(dataFile)}
      }
      
      if(r != 0){
        
        if(!(file.exists(dataFile))){
        data = genererEchantillon_new(n,d,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r, id_outliers =  outlier_sets[[m]])
        save(data,file = dataFile)
        }
        else{
          print("File exists")
          load(dataFile)
            }

      
      
    }}
  
}}