#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Weiszfeld geometric median computation
 //'
 //' @description Compute the median of a matrix X using the Weiszfeld algorithm
 //'
 //' \eqn{m_{t+1} = \frac{\sum_{k=1}{^n}X_k/\| X_k-m_t\|}{\sum_{k=1}{^n}1/\| X_k-m_t\|}}
 //'
 //' @param X A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are \eqn{n} observations in \eqn{\mathbf R^p}.
 //' @param init_median Initialisation of the median
 //' @param weights Weights
 //' @param epsilon Tolerance
 //' @param nitermax Maximum number of iterations
 //'
 //' @return A list with the median estimator and the number of iterations
 //'
 //' @export
 // [[Rcpp::export]]
 Rcpp::List WeiszfeldMedianRcpp(const arma::mat& X,
                                const arma::rowvec& init_median,
                                const arma::rowvec& weights,
                                double epsilon = 1e-08,
                                int nitermax = 100){
   // X : n * p  matrix
   // Inputs
   const int n = X.n_rows ;
   const int p = X.n_cols ;
   // Containers and intialisation of the algorithm
   arma::rowvec meanvec(p);
   arma::rowvec w(n);
   arma::rowvec medvec(p);
   double diffxn, normxm = 10;
   int iter = 0;
   medvec = init_median;
   meanvec = init_median;
   /* Boucle Weiszfeld */
   while (iter < nitermax and normxm > epsilon)
   {
     for (int it=0 ; it < n ; it++){
       diffxn = norm(X.row(it)-meanvec);
       w(it) = (diffxn > 0) ? weights(it)/diffxn : 0;
     }
     w = w/sum(w); // normalisation
     medvec = w*X; // geometric median update, element by element product
     normxm = norm(medvec-meanvec)/sqrt(double(p)); // stopping criterion
     meanvec = medvec;
     iter++;
   }
   // Returns ;
   Rcpp::List ret ;
   ret["estimator"] = medvec ;
   ret["iter"] = iter ;
   return Rcpp::wrap(ret);
 }

//' @title Weiszfeld geometric median covariation matrix computation
    //'
    //' @description Compute the median covariation matrix of a matrix X using the Weiszfeld algorithm
    //'
    //' TODO comparer avec notations papier --> ajouter les poids w_k comme ds med
    //'
    //' \eqn{V_{t+1} = \frac{\sum_{k=1}{^n}\| (X_k-m^*)(X_k-m^*)^T-V_t\|_F^{-1}(X_k-m^*)(X_k-m^*)^T}{\sum_{k=1}{^n}\| (X_k-m^*)(X_k-m^*)^T-V_t\|_F^{-1}}}
    //'
    //' @param X A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are \eqn{n} observations in \eqn{\mathbf R^p}.
    //' @param median_est Estimation of the median \eqn{m^*}, typically an output of \code{\link{WeiszfeldMedianRcpp}}
    //' @param init_median_cov Initialisation of the median covariation matrix
    //' @param weights Weights
    //' @param epsilon Tolerance
    //' @param nitermax Maximum number of iterations
    //'
    //' @return A list with the median covariation matrix, the number of iterations and the weigths
    //'
    //' @export
    // [[Rcpp::export]]
    Rcpp::List WeiszfeldMedianCovarianceRcpp(const arma::mat& X,
                                             const arma::rowvec& median_est,
                                             const arma::mat& init_median_cov,
                                             const arma::rowvec& weights,
                                             double epsilon=1e-08,
                                             int nitermax=100)
    {
      // X : n * p  matrix
      // Inputs
      const int n = X.n_rows ;
      const int p = X.n_cols ;
      // Containers
      arma::mat Xcent(n,p); // TODO on recentre X ? et pas de normalisation sur la var ?
      for (int it=0 ; it < n ; it++)
      {
        Xcent.row(it) = X.row(it)-median_est;
      }
      // arma::mat medinit = arma::trans(Xcent) * Xcent/n;
      arma::mat medinit = init_median_cov;
      arma::mat  medest(p,p);
      // medest = init_median_cov;
      arma::rowvec w(n);
      double diffxn, normxm = 1;
      int iter = 0;
      
      while (iter < nitermax and normxm > epsilon )
      {
        medest.fill(0.0);
        for (int it=0 ; it < n ; it++){
          diffxn = arma::norm(arma::trans(Xcent.row(it))*Xcent.row(it)-medinit,"fro");
          w(it) = (diffxn > 0) ? weights(it)/diffxn : 0;
        }
        w = w/sum(w);
        for (int it=0 ; it < n ; it++)
        {
          medest += w(it)*arma::trans(Xcent.row(it))*Xcent.row(it);
        }
        normxm = norm(medest-medinit,"fro")/p;
        medinit = medest;
        iter++;
      }
      // Returns ;
      Rcpp::List ret ;
      ret["estimator"] = medest ;
      ret["iter"] = iter ;
      ret["weights"] = w; // TODO wieghts, voir ou on l'utilise
      return Rcpp::wrap(ret);
    }


//' @title Averaged Stochastic Gradient geometric median computation
    //'
    //' @description Compute the median of a matrix X using the Averaged Stochastic Gradient algorithm
    //'
    //' \eqn{m_{k+1} = m_k + \gamma_{k+1}\frac{X_{k+1} − m_k}{\| X_{k+1} − m_k\|}}
    //'
    //' \eqn{\bar m_{k+1} = \bar m_k + \frac 1 {k+1} (m_{k+1} − \bar m_k )}
    //'
    //' @param X  A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are n observations in \eqn{\mathbf R^p}.
    //' @param init_median  Initialisation of the median
    //' @param weights Weights TODO ?? dans quel algo ?
    //' @param gamma TODO
    //' @param alpha TODO --> gamma ds le papier
    //' @param nstart TODO
    //' @param epsilon TODO
    //'
    //' @return A list with the computed median
    //'
    //' @export
    // [[Rcpp::export]]
    Rcpp::List ASGMedianRcpp(const arma::mat& X,
                             const arma::rowvec& init_median,
                             const arma::rowvec& weights,
                             double gamma = 2,
                             double alpha = 0.75,
                             int nstart = 1,
                             double epsilon = 1e-8)
    {
      try {
        if (alpha < 0.5 || alpha > 1) {
          throw std::range_error("alpha must be in [0.5, 1]");
        }
        if (gamma < 0 ) {
          throw std::range_error("gamma must be positive");
        }
      } catch(std::exception &ex) {
        forward_exception_to_r(ex);
      } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
      }
      
      // X : n * p  matrix
      // Inputs
      const int n = X.n_rows ;
      const int p = X.n_cols ;
      // Containers and intialisation of the algorithm
      arma::rowvec medvec = init_median ;
      arma::rowvec medrm = init_median;
      double w, normxm ;
      // Number of replications of the algorithm
      for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
        // Stochastic gradient algorithms
        for (int it = 1 ; it < n ; it++)
        {
          normxm = arma::norm(X.row(it)-medrm);
          if (normxm > epsilon) {
            w = weights(it)*sqrt(double(p))*gamma*pow(double(it+1),-alpha)/normxm;
            medrm += w * (X.row(it)-medrm) ;
          }
          medvec += (medrm-medvec)/(it+1);
        }
      }
      // Returns ;
      Rcpp::List ret ;
      ret["estimator"] = medvec ;
      
      return Rcpp::wrap(ret);
    }

//' @title Averaged Stochastic Gradient geometric median covariation matrix computation
    //'
    //' @description Compute the median covariation matrix of a matrix X using the Averaged Stochastic Gradient algorithm
    //'
    //' TODO -V_n ??
    //' \eqn{V_{k+1} = V_k + \gamma_{k+1}\frac{(X_{k+1} − \bar m_k)(X_{k+1} − \bar m_k)^T-V_n }{\| (X_{k+1} − \bar m_k)(X_{k+1} − \bar m_k)^T-V_n \|_F}}
    //'
    //' \eqn{\bar V_{k+1} = \bar V_k + \frac 1 {k+1} (V_{k+1} − \bar V_k )}
    //'
    //' @param X  A numerical matrix corresponding to the data of size \eqn{n\times p}. The rows are n observations in \eqn{\mathbf R^p}.
    //' @param median_est Estimation of the median \eqn{m^*}, typically an output of \code{\link{ASGMedianRcpp}}
    //' @param init_median_cov  Initialisation of the median covariance matrix
    //' @param weights Weights TODO ?? dans quel algo ?
    //' @param gamma TODO
    //' @param alpha TODO --> gamma ds le papier
    //' @param nstart TODO
    //'
    //' @return A list with the computed median covariation matrix estimator
    //'
    //' @export
    // [[Rcpp::export]]
    Rcpp::List ASGMedianCovarianceRcpp(const arma::mat&  X,
                                       const arma::rowvec&  median_est,
                                       const arma::mat&  init_median_cov,
                                       const arma::rowvec& weights,
                                       double gamma = 2,
                                       double alpha = 0.75,
                                       int nstart = 1){
      // X : n * p  matrix
      // Inputs
      const int n = X.n_rows ;
      const int p = X.n_cols ;
      // Containers
      arma::mat medav(p,p);
      medav=init_median_cov;
      arma::mat  diffmat(p,p), medrm(p,p),diffmed(p,p)  ;
      double nrmrm , w ;
      
      // Initialization of the algorithm
      //    diffmed = X.row(0)-median_est ;
      //    medav =   arma::mat init_median_cov;
      medrm = medav ;
      
      // Number of replications of the algorithm
      for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
        // Stochastic gradient algorithms
        for (int it = 1 ; it < n ; it++)
        {
          diffmed = X.row(it)-median_est ;
          diffmat = arma::trans(diffmed)*diffmed;
          diffmat -= medrm ;
          nrmrm = arma::norm(diffmat,"fro") ; // Frobenius norm divided by the dimension
          //           nrmrm = arma::sqrt(arma::sum(arma::square(diffmat)))/p ;
          w = p *weights(it) * gamma * pow(double(it+1),-alpha) * pow(nrmrm,-1);
          medrm +=   w*diffmat ;
          medav += (medrm-medav)/(it+1);
        }
      }
      
      Rcpp::List ret ;
      ret["estimator"] = medav ;
      
      return Rcpp::wrap(ret);
    }



//' @title Robbins-Monro method for the estimation of the eigen values of the variance-covariance matrix
 //'
 //' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
 //' the eigen values \eqn{\lambda} of a variance-covariance matrix using the Robbins-Monro method.
 //'
 //' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
 //'
 //' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
 //' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
 //' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
 //' @param delta Vector of size p of the eigen values of the median covariation matrix
 //' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta
 //' @param gamma Robbins-Monro parameter (0.75 by default)
 //' @param c Robbins-Monro parameter (2 by default)
 //' @param w Robbins-Monro parameter (2 by default)
 //' @param init_bar, A VOIR AVEC DAPHNE: c'est un extraparamèetre dont on se sert pour l'estimation en ligne, on est obligé de le donner?
 //' @param c_bar, A VOIR AVEC DAPHNE: c'est un extraparamèetre dont on se sert pour l'estimation en ligne, on est obligé de le donner?
 //' @param c_tilde, A VOIR AVEC DAPHNE: c'est un extraparamèetre dont on se sert pour l'estimation en ligne, on est obligé de le donner?
 //' @param c_bar, A VOIR AVEC DAPHNE: c'est un extraparamèetre dont on se sert pour l'estimation en ligne, on est obligé de le donner?
 //' @param sumlog Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
 //' @param out_index Indexes of the iterations for which we want to output the values of the estimates, default is the number of observations N
 //'
 //' @return A list containing the estimated eigen values \eqn{\lambda}, the number of iterations,
 //' and the values of the estimates for different iterations
 
 // [[Rcpp::export]]
 Rcpp::List robbinsMCRcpp(const arma::mat& U,
                          const arma::vec& delta,
                          arma::vec init,
                          arma::vec init_bar,
                          double gamma = 0.75,
                          double c = 2,
                          double w = 2,
                          double c_bar = 0,
                          double c_tilde = 0,
                          double sumlog = 1,
                          double epsilon = 1e-8,
                          IntegerVector out_index = IntegerVector::create()
 ) {
   
   int N = U.n_rows;
   int p = U.n_cols;
   
   arma::vec lambda = init;
   arma::vec lambda_bar = init_bar;
   
   arma::mat lambdalist(p, 0);
   arma::mat lambdabarlist(p, 0);
   
   for (int i = 0; i < N; ++i) {
     arma::rowvec U_i = U.row(i);
     arma::vec U_i_sq = trans(U_i % U_i);  // U_i^2
     
     arma::vec lambda_Ui2 = lambda % U_i_sq;
     arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;
     
     double E2 = std::pow(
       arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2)
       + arma::as_scalar(lambda_Ui2.t() * lambda_Ui2)
       - arma::accu(lambda_Ui2 % lambda_Ui2),
       -0.5);
       
       arma::vec E1 = U_i_sq * E2;
       
       arma::vec lambda_old = lambda;
       
       double i_gamma = std::pow(i + 1 + c_tilde, -gamma);
       lambda = lambda - c * i_gamma * lambda % E1 + c * i_gamma * delta * E2;
       
       sumlog = sumlog + std::pow(std::log(i + 2 + c_bar), w);
       lambda_bar = lambda_bar +  std::pow(std::log(i + 2 + c_bar), w) / sumlog * (lambda - lambda_bar);
       
       double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));
       
       // If index i is in out_index, save
       if (std::find(out_index.begin(), out_index.end(), i + 1) != out_index.end()) {
         lambdalist.insert_cols(lambdalist.n_cols, lambda);
         lambdabarlist.insert_cols(lambdabarlist.n_cols, lambda_bar);
       }
       
       if (eps < epsilon) break;
   }
   
   return List::create(
     Named("vp") = lambda_bar,
     Named("niter") = sumlog,
     Named("lambdalist") = lambdalist,
     Named("vplist") = lambdabarlist,
     Named("lambda") = lambda
   );
 }



//' @title Fix point method for the estimation of the eigen values of the variance-covariance matrix
 //'
 //' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
 //' the eigen values \eqn{\lambda} of a variance-covariance matrix using the fix point method.
 //'
 //' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
 //'
 //' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
 //' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
 //' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
 //' @param delta Vector of size p of the eigen values of the median covariation matrix
 //' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta
 //' @param niter Maximum number of iterations for the gradient descent algorithm, by default 10
 //' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
 //' @param out_index Indexes of the iterations for which we want to output the values of the estimates, default is the last iteration
 //'
 //' @return A list containing the estimated eigen values \eqn{\lambda}, the number of iterations,
 //' and the values of the estimates for different iterations
 //'
 
 
 
 // [[Rcpp::export]]
 Rcpp::List fixMCRcpp(const arma::mat& U,
                      const arma::vec& delta,
                      arma::vec init,
                      int niter = 10,
                      double epsilon = 1e-8,
                      IntegerVector out_index = IntegerVector::create())
 {
   int N = U.n_rows;
   int p = U.n_cols;
   
   arma::vec lambda = init;
   arma::mat lambdalist(p, 0); // Each column will be a lambda
   
   for (int l = 0; l < niter; ++l) {
     arma::vec E1(p, arma::fill::zeros);
     double E2 = 0.0;
     
     for (int i = 0; i < N; ++i) {
       arma::rowvec U_i = U.row(i);
       arma::vec U_i_sq = trans(U_i % U_i); // element-wise square
       
       arma::vec lambda_Ui2 = lambda % U_i_sq;
       arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;
       
       double term1 = arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2);
       double term2 = arma::as_scalar(lambda_Ui2.t() * lambda_Ui2);
       double term3 = arma::accu(lambda_Ui2 % lambda_Ui2);
       
       double Ei = std::pow(term1 + term2 - term3, -0.5);
       
       E1 += U_i_sq * Ei;
       E2 += Ei;
     }
     
     arma::vec lambda_old = lambda;
     lambda = E2 * delta / E1;
     
     double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));
     
     // Save if iteration is in out_index
     if (std::find(out_index.begin(), out_index.end(), l + 1) != out_index.end()) {
       lambdalist.insert_cols(lambdalist.n_cols, lambda);
     }
     
     if (eps < epsilon) break;
   }
   
   return List::create(
     Named("vp") = lambda,
     Named("niter") = niter,
     Named("vplist") = lambdalist
   );
 }








//' @title Gradient descent for the estimation of the eigen values of the variance-covariance matrix
 //'
 //' @description Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
 //' the eigen values \eqn{\lambda} of a variance-covariance matrix using the gradient descent method.
 //'
 //' @param U Matrix of size N*p corresponding to \eqn{\Sigma^{-1/2}(X-\mu)}. The rows are the observations.
 //'
 //' - In a gaussian model typically `U = matrix(rnorm(N*p),ncol=p)`
 //' - In a Student model `U <- matrix(rnorm(N*p)/sqrt(rchisq(1,df=df))*sqrt((df-2)), ncol=p))`
 //' - In a Laplace model `U <- LaplacesDemon::rmvl(N,mu=rep(0,p),Sigma=diag(p))`
 //' @param delta Vector of size p of the eigen values of the median covariation matrix
 //' @param init Initial value of the vector of eigen values of the variance-covariance matrix, by default equal to delta
 //' @param niter Maximum number of iterations for the gradient descent algorithm, by default 10
 //' @param epsilon Stopping criterion: the algorithm stops when the difference between two iterations is less than epsilon, by default 1e-08
 //' @param step Step of the gradient descent, by default 1
 //' @param out_index Indexes of the iterations for which we want to output the values of the estimates, default is the last iteration
 //'
 //' @return A list containing the estimated eigen values \eqn{\lambda}, the number of iterations,
 //' and the values of the estimates for different iterations
 //'
 
 // [[Rcpp::export]]
 Rcpp::List gradMCRcpp(const arma::mat& U,
                       const arma::vec& delta,
                       arma::vec init,
                       int niter = 10,
                       double epsilon = 1e-8,
                       NumericVector step = NumericVector::create(),
                       IntegerVector out_index = IntegerVector::create())
 {
   int N = U.n_rows;
   int p = U.n_cols;
   
   arma::vec lambda = init;
   arma::mat lambdalist(p, 0); // Each column will be a lambda
   
   // Default steps if none provided
   if (step.size() == 0) {
     step = rep(1.0, niter);
   }
   
   for (int l = 0; l < niter; ++l) {
     arma::vec E1(p, arma::fill::zeros);
     double E2 = 0.0;
     
     for (int i = 0; i < N; ++i) {
       arma::rowvec U_i = U.row(i);
       arma::vec U_i_sq = trans(U_i % U_i); // U_i^2
       
       arma::vec lambda_Ui2 = lambda % U_i_sq;
       arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;
       
       double term1 = arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2);
       double term2 = arma::as_scalar(lambda_Ui2.t() * lambda_Ui2);
       double term3 = arma::accu(lambda_Ui2 % lambda_Ui2);
       
       double Ei = std::pow(term1 + term2 - term3, -0.5);
       
       E1 += U_i_sq * Ei;
       E2 += Ei;
     }
     
     arma::vec lambda_old = lambda;
     
     double step_l = step[l];
     lambda = lambda - (step_l / N) * lambda % E1 + (step_l / N) * delta * E2;
     
     double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));
     
     // Save if iteration is in out_index
     if (std::find(out_index.begin(), out_index.end(), l + 1) != out_index.end()) {
       lambdalist.insert_cols(lambdalist.n_cols, lambda);
     }
     
     // Optional stopping criterion:
     if (eps < epsilon) break;
   }
   
   return List::create(
     Named("vp") = lambda,
     Named("niter") = niter,
     Named("vplist") = lambdalist
   );
 }


// [[Rcpp::export]]
List update_median_covarianceRcpp(const arma::mat& X,
                                  arma::rowvec m,
                                  arma::rowvec moyennem,
                                  arma::mat V,
                                  arma::mat moyenneV,
                                  int Ninit,
                                  int niterr,
                                  int batch,
                                  double gamman,
                                  double w,
                                  double sslog) {
  
  int d = X.n_cols;
  arma::rowvec gradm = arma::zeros<arma::rowvec>(d);
  arma::rowvec x_l, diff;
  
  int start_index = Ninit - 1 + (niterr - 1) * batch;
  
  // ----- 1. Updating the estimates of the median -----
  if (batch == 1) {
    x_l = X.row(start_index);
    diff = x_l - m;
    gradm = diff / std::sqrt(arma::accu(arma::square(diff)));
  } else {
    for (int l = 0; l < batch; ++l) {
      x_l = X.row(start_index + l);
      diff = x_l - m;
      gradm += diff / std::sqrt(arma::accu(arma::square(diff)));
    }
    gradm /= batch;
  }
  
  m = m + gamman * gradm;
  
  double weight = std::pow(std::log(niterr + 1), w);
  sslog += weight;
  moyennem += (weight / sslog) * (m - moyennem);
  
  // ----- 2. Updating the estimates of the MCM -----
  arma::mat gradV = arma::zeros<arma::mat>(d, d);
  arma::rowvec centered;
  arma::mat outer;
  
  if (batch == 1) {
    x_l = X.row(start_index);
    centered = x_l - moyennem;
    outer = centered.t() * centered;
    arma::mat delta = outer - V;
    double normF = arma::norm(delta, "fro");
    gradV = delta / normF;
  } else {
    for (int l = 0; l < batch; ++l) {
      x_l = X.row(start_index + l);
      centered = x_l - moyennem;
      outer = centered.t() * centered;
      arma::mat delta = outer - V;
      double normF = arma::norm(delta, "fro");
      gradV += delta / normF;
    }
    gradV /= batch;
  }
  
  V = V + gamman * gradV;
  moyenneV += (weight / sslog) * (V - moyenneV);
  
  return List::create(
    Named("m") = m,
    Named("moyennem") = moyennem,
    Named("V") = V,
    Named("moyenneV") = moyenneV,
    Named("sslog") = sslog
  );
}




//' Internal
 //' To normalize vectors
 //'
 
 // [[Rcpp::export]]
 arma::mat normalize_columnsRcpp(const arma::mat& V) {
   arma::mat Vnorm = V;
   arma::rowvec norms = sqrt(sum(square(V), 0));
   Vnorm.each_row() /= norms;
   return Vnorm;
 }




//' Internal
 //' To normalize vectors
 // [[Rcpp::export]]
 arma::mat reconstruct_covarianceRcpp(const arma::mat& VP, const arma::rowvec& lambda) {
   arma::mat A = VP * arma::diagmat(arma::sqrt(lambda.t()));  // lambda.t() : rowvec → vec
   return A * A.t();
 }

//' Internal
 //' Mahalanobis
 //'
 
 // [[Rcpp::export]]
 double mahalanobis_generalizedRcpp(const arma::rowvec& x,
                                    const arma::rowvec& moyennem,
                                    const arma::mat& eigvecs,
                                    const arma::rowvec& lambdaInit) {
   arma::rowvec centered = x - moyennem;
   arma::vec proj = eigvecs.t() * centered.t();  // Projeter x - moyennem sur les vecteurs propres
   arma::vec scaled = arma::square(proj) / lambdaInit.t();  // lambdaInit est un rowvec → transpose en vec
   return arma::accu(scaled);
 }


