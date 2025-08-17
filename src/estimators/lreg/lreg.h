/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 29/04/2025
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

  void param() {

    res = Y - predictors * transparameters;

  }

  void F() {

    f = arma::accu(res % res);

  }

  void G() {

    grad = -2 * predictors.t() * res;

  }

  void dG() {

    dX = dparameters;
    dg = 2 * predictors.t() * (predictors * dX);

  }

  void H() {

    // Rcpp::stop("H not available");
    // hess = 2 * predictors.t() * predictors;

  }

  void E() {}

  void M() {}

  void outcomes() {}

};

lreg* choose_lreg(const Rcpp::List& estimator_setup) {

  lreg* myestimator = new lreg();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat Y = estimator_setup["Y"];
  arma::mat predictors = estimator_setup["predictors"];

  myestimator->indices = indices;
  myestimator->Y = Y;
  myestimator->predictors = predictors;

  return myestimator;

}
