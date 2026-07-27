#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Internal Robbins-Monro method for the estimation of the eigen values of the variance-covariance matrix
//'
//' Given the eigen values \eqn{\delta} of the median covariation matrix, this function estimates
//' the eigen values \eqn{\lambda} of a variance-covariance matrix using the Robbins-Monro method and Monte Carlo approximation.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List robbinsMCList(const arma::mat &U,
                        const arma::vec &delta,
                        Rcpp::Nullable<arma::vec> init = R_NilValue,
                        Rcpp::Nullable<arma::vec> init_bar = R_NilValue,
                        double gamma = 0.75,
                        double c = 2,
                        double w = 2,
                        double c_bar = 0,
                        double c_tilde = 0,
                        double sumlog = 1,
                        double epsilon = 1e-8,
                        IntegerVector out_index = IntegerVector::create())
{
  int N = U.n_rows;
  int p = U.n_cols;

  arma::vec lambda;
  arma::vec lambda_bar;

  if (init.isNotNull())
    lambda = Rcpp::as<arma::vec>(init);
  else
    lambda = delta;

  if (init_bar.isNotNull())
    lambda_bar = Rcpp::as<arma::vec>(init_bar);
  else
    lambda_bar = delta;

  arma::mat lambdalist(p, 0);
  arma::mat lambdabarlist(p, 0);

  for (int i = 0; i < N; ++i)
  {
    arma::rowvec U_i = U.row(i);
    arma::vec U_i_sq = trans(U_i % U_i); // U_i^2

    arma::vec lambda_Ui2 = lambda % U_i_sq;
    arma::vec delta_minus_lambda_Ui2 = delta - lambda_Ui2;

    double E2 = std::pow(arma::dot(delta_minus_lambda_Ui2, delta_minus_lambda_Ui2) + arma::as_scalar(lambda_Ui2.t() * lambda_Ui2) - arma::accu(lambda_Ui2 % lambda_Ui2),-0.5);

    arma::vec E1 = U_i_sq * E2;
    arma::vec lambda_old = lambda;

    double i_gamma = std::pow(i + 1 + c_tilde, -gamma);
    lambda = lambda - c * i_gamma * lambda % E1 + c * i_gamma * delta * E2;

    sumlog = sumlog + std::pow(std::log(i + 2 + c_bar), w);
    lambda_bar = lambda_bar + std::pow(std::log(i + 2 + c_bar), w) / sumlog * (lambda - lambda_bar);

    double eps = std::sqrt(arma::accu(arma::square(lambda - lambda_old)));

    // If index i is in out_index, save
    if (std::find(out_index.begin(), out_index.end(), i + 1) != out_index.end())
    {
      lambdalist.insert_cols(lambdalist.n_cols, lambda);
      lambdabarlist.insert_cols(lambdabarlist.n_cols, lambda_bar);
    }

    if (eps < epsilon)
      break;
  }

  return List::create(
    Named("lambda_bar") = lambda_bar,
    Named("lambda") = lambda,
    Named("lambda_barlist") = lambdabarlist,
    Named("lambdalist") = lambdalist);
}

//' @keywords internal
// [[Rcpp::export]]
List updateMedianCovarianceRcpp(const arma::mat& X,
                                arma::rowvec m,
                                arma::rowvec averaged_m,
                                arma::mat V,
                                arma::mat averaged_V,
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
  averaged_m += (weight / sslog) * (m - averaged_m);

  // ----- 2. Updating the estimates of the MCM -----
  arma::mat gradV = arma::zeros<arma::mat>(d, d);
  arma::rowvec centered;
  arma::mat outer;

  if (batch == 1) {
    x_l = X.row(start_index);
    centered = x_l - averaged_m;
    outer = centered.t() * centered;
    arma::mat delta = outer - V;
    double normF = arma::norm(delta, "fro");
    gradV = delta / normF;
  } else {
    for (int l = 0; l < batch; ++l) {
      x_l = X.row(start_index + l);
      centered = x_l - averaged_m;
      outer = centered.t() * centered;
      arma::mat delta = outer - V;
      double normF = arma::norm(delta, "fro");
      gradV += delta / normF;
    }
    gradV /= batch;
  }

  V = V + gamman * gradV;
  averaged_V += (weight / sslog) * (V - averaged_V);

  return List::create(
    Named("m") = m,
    Named("averaged_m") = averaged_m,
    Named("V") = V,
    Named("averaged_V") = averaged_V,
    Named("sslog") = sslog
  );
}


//' Internal
//' To normalize vectors
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat normalizeColumnsRcpp(const arma::mat& V) {
  return arma::normalise(V, 2, 0);
}


//' Internal
//' To build covariance
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat buildCovarianceRcpp(const arma::mat& VP, const arma::rowvec& lambda) {
  arma::mat A = VP * arma::diagmat(arma::sqrt(lambda.t()));  // lambda.t() : rowvec → vec
  return A * A.t();
}

//' Internal
//' Mahalanobis
//'
//' @keywords internal
// [[Rcpp::export]]
double mahalanobisGeneralizedRcpp(const arma::rowvec& x,
                                  const arma::rowvec& averaged_m,
                                  const arma::mat& eigvecs,
                                  const arma::rowvec& lambdaInit) {
  arma::rowvec centered = x - averaged_m;
  arma::vec proj = eigvecs.t() * centered.t();  // Projeter x - averaged_m sur les vecteurs propres
  arma::vec scaled = arma::square(proj) / lambdaInit.t();  // lambdaInit est un rowvec → transpose en vec
  return arma::accu(scaled);
}


//' @title Stochastic median clustering
//'
//' @description Gmedian
//'
//' @param X a matrix
//' @param Xtot a matrix
//' @param centers a matrix
//' @param gamma a double
//' @param alpha a double
//'
//' @return A list
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List stoKmed_rcpp(const arma::mat X, const arma::mat Xtot, const arma::mat centers, const double gamma=2, const double alpha = 0.75) // TODO --> passer en internal
{
  // Inputs
  // X : n * p  matrix
  // Xtot: n * p  matrix
  // centers : k * p matrix
  int n = X.n_rows,  p = X.n_cols;
  int k = centers.n_rows; /* Nombre de classes */

  int j_best;
  arma::mat  centersrm = centers;
  arma::mat centersav = centers;
  arma::vec nc(k);
  arma::vec wss(k);

  nc.fill(1); /* initialisation des tailles de classe, pour éviter la division par zéro */
  arma::vec clusterindex(n);

  double poids, best, dd;

  for(int i = 0; i < n; i++) { /* boucle sur les individus */
    best = R_PosInf;
    for(int j = 0; j < k; j++) {  /*calcul du centre le plus proche*/
    dd = arma::norm(X.row(i)-centersrm.row(j))/sqrt(double(p));
      if(dd < best) {
       best = dd;
       j_best = j;
      }
    }
    clusterindex(i) = j_best+1; /* affectation de l'ind. dans la nouvelle classe */ //TODO : je l'écrase ensuite, pas besoin de le remplir ici ?

    /* mise a jour de RM*/
    if (best>0){
     poids = gamma/(pow(double(nc[j_best]),alpha)*sqrt(best));
     centersrm.row(j_best) +=   (X.row(i) - centersrm.row(j_best))*poids;
    }
    nc(j_best)++; /* incrementation de la taille de la classe */
    /* mise a jour de AVeraged*/
    centersav.row(j_best) +=  (centersrm.row(j_best)-centersav.row(j_best))/nc[j_best];
  } /* fin boucle sur individus */

  /* Calcul des distances intra-classes */
  /* initialisation des classes, des tailles de classe et des erreurs */
  wss.fill(0.0);
  nc.fill(0);
  clusterindex.fill(0);

  /* boucle sur tous les n+k individus de la table xtot*/
  for(int i = 0; i < n; i++) {
    best = R_PosInf;
    for(int j = 0; j < k; j++) { /*calcul du centre j le plus proche*/
      dd = arma::norm(Xtot.row(i)-centersav.row(j))/sqrt(double(p));

      if(dd < best) {
        best = dd;
        j_best = j;
      }
    }

    clusterindex(i) = j_best+1; /* affectation de l'ind. dans la nouvelle classe */
    nc(j_best)++;
    wss(j_best) += sqrt(best);
  }

  // Returns ;
  Rcpp::List ret ;
  ret["wss"] = wss ;
  ret["centers"] = centersav ;
  ret["cl"] = clusterindex ;
  ret["nc"] = nc ;
  return Rcpp::wrap(ret);
}
