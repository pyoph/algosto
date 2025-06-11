#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List RobbinsMC2_cpp(int mc_sample_size = 10000, 
                    NumericVector vp = NumericVector::create(), 
                    double epsilon = 1e-8, 
                    double alpha = 0.75, 
                    double c = 1.0, 
                    int w = 2, 
                    NumericVector samp = NumericVector::create(), 
                    NumericVector init = NumericVector::create(), 
                    NumericVector initbarre = NumericVector::create(), 
                    double cbarre = 0, 
                    double ctilde = 0, 
                    double slog = 1)   {
  
  int p = vp.size();
  int niter = 0;
  NumericVector vp2 = clone(vp);
  NumericVector lambda = clone(vp);
  //Rcout << "init received from R (size: " << init.size() << "): " << init << "\n";
  
  NumericMatrix lambdalist(0, p); // Stocke les valeurs de lambda
  NumericMatrix vplist(0, p);     // Stocke les valeurs de vp2
  
  std::vector<NumericVector> lambdalist_vec;
  std::vector<NumericVector> vplist_vec;

  // Génération de Y (matrice de valeurs normales aléatoires)
  
    arma::mat Y = arma::randn(mc_sample_size, p);
  
  
  for (int i = 0; i < mc_sample_size; i++) {
    arma::rowvec Z = Y.row(i);
    NumericVector Zvec = as<NumericVector>(wrap(Z));
    
// Rcout<<Zvec;
    
    double norm_factor = sqrt(sum(pow(vp - lambda * pow(Zvec, 2), 2)));
    NumericVector E1 = pow(Zvec, 2) / norm_factor;
    double E2 = 1.0 / norm_factor;
    //Rcout << "Iteration " << i << ":\n";
    //Rcout << "vp: " << vp << "\n";
    //Rcout << "lambda: " << lambda << "\n";
    //Rcout << "Zvec: " << Zvec << "\n";

    // Mise à jour de lambda
    NumericVector lambda_prev = clone(lambda);
    if (norm_factor < 1e-10) norm_factor = 1e-10;  
    lambda = lambda - c * pow(i + ctilde, -alpha) * lambda * E1 + 
      c * pow(i + ctilde, -alpha) * vp * E2;
    //Rcout << "Iteration: " << i << " | norm_factor: " << norm_factor << " | E1: " << E1 << " | E2: " << E2 << "\n";
    // Mise à jour de slog et vp2
    slog += log(i + 1 + cbarre) * w;
    vp2 = vp2 + log(i + 1 + cbarre) * w * (1.0 / slog) * (lambda - vp2);
    
    // Calcul de l'erreur et stockage des valeurs
    double eps = sqrt(sum(pow(vp2 - lambda_prev, 2)));
    
    if (is_true(any(samp == (i + 1)))) { // Vérifie si i+1 est dans samp
      lambdalist_vec.push_back(lambda); // Ajoute à la liste
      vplist_vec.push_back(vp2);
    }
    
    // Arrêt anticipé si convergence
    if (eps < epsilon) break;
  }
  
  return List::create(Named("vp") = vp2, 
                      Named("niter") = niter, 
                      Named("lambdalist") = lambdalist, 
                      Named("vplist") = vplist, 
                      Named("lambda") = lambda);
}
