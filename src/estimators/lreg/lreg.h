/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 15/02/2026
 */

/*
 * Least-squares
 */

class lreg: public estimators {

public:

  double loss;
  arma::uvec indices;
  arma::vec y, beta;
  arma::mat X, res;

  void param(arguments_optim& x) {

    beta = x.transparameters(indices);
    res = y - X * beta;

  }

  void F(arguments_optim& x) {

    loss = arma::accu(res % res);
    x.f += loss;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices) += arma::vectorise(-2 * X.t() * res);

  }

  void dG(arguments_optim& x) {

    arma::vec dbeta = x.dtransparameters(indices);
    x.dgrad.elem(indices) += arma::vectorise(2 * X.t() * (X * dbeta));

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = loss;

  }

};

lreg* choose_lreg(const Rcpp::List& estimator_setup) {

  lreg* myestimator = new lreg();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat y = estimator_setup["y"];
  arma::mat X = estimator_setup["X"];

  myestimator->indices = indices[0];
  myestimator->y = y;
  myestimator->X = X;

  return myestimator;

}
