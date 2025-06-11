#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List RobbinsMC2_cpp(int mc_sample_size = 10000, 
                          arma::vec vp = arma::vec(), 
                          double epsilon = 1e-8, 
                          double alpha = 0.75, 
                          double c = 1.0, 
                          int w = 2, 
                          Rcpp::IntegerVector samp = Rcpp::IntegerVector::create(), 
                          arma::vec init = arma::vec(), 
                          arma::vec initbarre = arma::vec(), 
                          double cbarre = 0, 
                          double ctilde = 0, 
                          double slog = 1.0) {
  
  int p = vp.n_elem;
  int niter = 0;
  
  arma::vec vp2 = vp;
  arma::vec lambda = vp;
  
  std::vector<arma::vec> lambdalist_vec;
  std::vector<arma::vec> vplist_vec;
  
  arma::mat Y = arma::randn(mc_sample_size, p);
  
  for (int i = 0; i < mc_sample_size; ++i) {
    arma::rowvec Z = Y.row(i);
    arma::vec Z2 = arma::square(Z.t());
    
    double norm_factor = std::sqrt(arma::accu(arma::square(vp - lambda % Z2)));
    if (norm_factor < 1e-10) norm_factor = 1e-10;
    
    arma::vec E1 = Z2 / norm_factor;
    double E2 = 1.0 / norm_factor;
    
    arma::vec lambda_prev = lambda;
    
    lambda -= c * std::pow(i + ctilde, -alpha) * lambda % E1
    - c * std::pow(i + ctilde, -alpha) * vp * E2;
    
    slog += std::log(i + 1 + cbarre) * w;
    vp2 += std::log(i + 1 + cbarre) * w * (1.0 / slog) * (lambda - vp2);
    
    double eps = std::sqrt(arma::accu(arma::square(vp2 - lambda_prev)));
    
    if (std::find(samp.begin(), samp.end(), i + 1) != samp.end()) {
      lambdalist_vec.push_back(lambda);
      vplist_vec.push_back(vp2);
    }
    
    if (eps < epsilon) break;
  }
  
  // Convert std::vector<arma::vec> to arma::mat
  arma::mat lambdalist(p, lambdalist_vec.size());
  arma::mat vplist(p, vplist_vec.size());
  
  for (size_t j = 0; j < lambdalist_vec.size(); ++j) {
    lambdalist.col(j) = lambdalist_vec[j];
    vplist.col(j) = vplist_vec[j];
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vp") = vp2,
    Rcpp::Named("niter") = niter,
    Rcpp::Named("lambdalist") = lambdalist.t(), // transpose to match R expectations
    Rcpp::Named("vplist") = vplist.t(),
    Rcpp::Named("lambda") = lambda
  );
}
