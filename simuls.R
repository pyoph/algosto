# Define output directory (relative path)
SimDir = "SimulatedData"
if (!dir.exists(SimDir)) dir.create(SimDir)


# Source parameters (must be in the same directory)
source("parametres.R")

setwd(SimDir)


##################Function to generate F_0 and F_1#################
genererEchantillon_new <- function(n, d, mu1, mu2, Sigma1, Sigma2, r, id_outliers = NULL) {
  
  # Initialisation
  Z <- matrix(0, n, d)
  labelsVrais <- rep(0, n)
  
  if (r == 0) {
    # Cas sans contamination
    Z <- mvrnorm(n, mu1, Sigma1)
    
  } else {
    # Cas avec contamination
    idx_out <- id_outliers
    idx_in  <- setdiff(1:n, idx_out)
    
    n1 <- length(idx_in)
    n2 <- length(idx_out)
    
    # Inliers
    Z[idx_in, ] <- mvrnorm(n1, mu1, Sigma1)
    labelsVrais[idx_in] <- 0
    
    # Outliers
    Z[idx_out, ] <- mvrnorm(n2, mu2, Sigma2)
    labelsVrais[idx_out] <- 1
  }
  
  return(list(Z = Z, labelsVrais = labelsVrais))
}

#######################Data generation#############################
for (sim in 1:simNb){
  for(sc in scenarios_1_param ){
    
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
      #setwd(SimDir)
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
