/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/04/2026
 */

/*
 * Logarithm of determinant of a matrix
 */

class logdetmat: public estimators {

public:

  int p;
  double f, tr, logdetw;
  arma::mat X, dX, Xinv;
  arma::uvec indices, lower_indices;

  void param(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices), p, p);
    tr = arma::trace(X);

  }

  void F(arguments_optim& x) {

    f = logdetw * (arma::log_det_sympd(X) - p*std::log(tr/p));
    x.f -= f;

  }

  void G(arguments_optim& x) {

    Xinv = arma::inv_sympd(X, arma::inv_opts::allow_approx);
    arma::mat Xtr(p, p, arma::fill::eye);
    Xtr.diag() *= p/tr;
    x.grad.elem(indices) -= logdetw * (Xinv - Xtr);

  }

  void dG(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices), p, p);

    double dtr = arma::trace(dX);
    arma::mat dXinv = -Xinv * dX * Xinv;
    arma::mat dXtr(p, p, arma::fill::eye);
    dXtr.diag() *= -p/(tr*tr)*dtr;
    arma::mat term = dXinv - dXtr;
    x.dgrad.elem(indices) -= logdetw * term;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    doubles[0] =  -f;       // loss   actual model
    doubles[1] =  0.00;     // loglik actual model
    doubles[2] =  0.00;     // loglik independence model
    doubles[3] =  0.00;     // loglik saturated model
    doubles[4] =  -f;       // penalty

  };

};

logdetmat* choose_logdetmat(const Rcpp::List& estimator_setup) {

  logdetmat* myestimator = new logdetmat();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::uvec lower_indices = estimator_setup["lower_indices"];
  double logdetw = estimator_setup["logdetw"];
  int p = estimator_setup["p"];

  arma::mat X(p, p, arma::fill::zeros);

  myestimator->indices = indices[0];
  myestimator->lower_indices = lower_indices;
  myestimator->X = X;
  myestimator->dX = X;
  myestimator->logdetw = logdetw;
  myestimator->p = p;

  return myestimator;

}
