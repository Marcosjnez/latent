/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

/*
 * Logarithm of determinant of a matrix
 */

class logdetmat: public estimators {

public:

  arma::mat X, dX, Xinv;
  arma::uvec lower_indices;
  double tr, logdetw;
  int p;
  double loglik = 0.00;

  void param(arguments_optim& x) {

    X.elem(lower_indices) = x.transparameters(indices[0]);
    X = arma::symmatl(X);
    tr = arma::trace(X);

  }

  void F(arguments_optim& x) {

    f = logdetw * (arma::log_det_sympd(X) - p*std::log(tr/p));
    x.f -= f;

  }

  void G(arguments_optim& x) {

    Xinv = arma::inv_sympd(X, arma::inv_opts::allow_approx);
    arma::mat Xinv2 = 2*Xinv;
    Xinv2.diag() *= 0.5;
    arma::mat Xtr(p, p, arma::fill::eye);
    Xtr.diag() *= p/tr;
    x.grad.elem(indices[0]) -= logdetw * (Xinv2.elem(lower_indices) -
      Xtr.elem(lower_indices));

  }

  void dG(arguments_optim& x) {

    dX.elem(lower_indices) = x.dtransparameters(indices[0]);
    dX = arma::symmatl(dX);

    double dtr = arma::trace(dX);
    arma::mat dXinv = -Xinv * dX * Xinv;
    arma::mat dXtr(p, p, arma::fill::eye);
    dXtr.diag() *= -p/(tr*tr)*dtr;
    dXinv *= 2;
    dXinv.diag() *= 0.5;
    arma::mat term = dXinv - dXtr;
    x.dgrad.elem(indices[0]) -= logdetw * (term.elem(lower_indices));

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] = f;
    doubles[1] = loglik;

  };

};

logdetmat* choose_logdetmat(const Rcpp::List& estimator_setup) {

  logdetmat* myestimator = new logdetmat();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::uvec lower_indices = estimator_setup["lower_indices"];
  double logdetw = estimator_setup["logdetw"];
  int p = estimator_setup["p"];

  arma::mat X(p, p, arma::fill::zeros);
  myestimator->indices = indices;
  myestimator->lower_indices = lower_indices;
  myestimator->X = X;
  myestimator->dX = X;
  myestimator->logdetw = logdetw;
  myestimator->p = p;

  return myestimator;

}
