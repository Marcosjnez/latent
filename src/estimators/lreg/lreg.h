/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Linear regression
 */

class lreg: public estimators {

public:

  arma::vec X;
  arma::vec dX;

  // Provide these in choose_estimator:
  arma::vec Y;
  arma::mat predictors;
  arma::vec res;
  // arma::uvec target_indices;
  // arma::uvec indices;

  void param() {

    X = parameters;

    res = Y - predictors * X;

  }

  void F() {

    f = arma::accu(res % res);

  }

  void G() {

    g = -2 * predictors.t() * res;

  }

  void dG() {

    dX = dparameters;
    dg = 2 * predictors.t() * (predictors * dX);

  }

  void H() {

    // Rcpp::stop("H not available");
    hessian = 2 * predictors.t() * predictors;

  }

  void E() {}

  void M() {}

  void outcomes() {}

};
