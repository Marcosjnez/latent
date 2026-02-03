/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2026
 */

/*
 * Linear regression
 */

class lreg: public estimators {

public:

  arma::mat X, res;
  arma::vec y, beta;

  void param(arguments_optim& x) {

    beta = x.transparameters(indices[0]);
    res = y - X * beta;

  }

  void F(arguments_optim& x) {

    f = arma::accu(res % res);
    x.f += f;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices[0]) += arma::vectorise(-2 * X.t() * res);

  }

  void dG(arguments_optim& x) {

    arma::vec dbeta = x.dtransparameters(indices[0]);
    x.dgrad.elem(indices[0]) += arma::vectorise(2 * X.t() * (X * dbeta));

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] =  f;

  }

};

lreg* choose_lreg(const Rcpp::List& estimator_setup) {

  lreg* myestimator = new lreg();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat y = estimator_setup["y"];
  arma::mat X = estimator_setup["X"];

  myestimator->indices = indices;
  myestimator->y = y;
  myestimator->X = X;

  return myestimator;

}
