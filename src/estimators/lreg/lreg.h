/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 31/08/2025
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

  void param(arguments_optim& x) {

    res = Y - predictors * transparameters;

  }

  void F(arguments_optim& x) {

    f = arma::accu(res % res);

  }

  void G(arguments_optim& x) {

    grad = -2 * predictors.t() * res;

  }

  void dG(arguments_optim& x) {

    dX = dparameters;
    dg = 2 * predictors.t() * (predictors * dX);

  }

  void H(arguments_optim& x) {

    // Rcpp::stop("H not available");
    // hess = 2 * predictors.t() * predictors;

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void outcomes(arguments_optim& x) {}

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
