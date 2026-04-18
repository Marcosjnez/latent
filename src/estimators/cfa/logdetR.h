/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 16/04/2026
 */

/*
 * Logarithm of determinant of a matrix
 */

class logdetR: public estimators {

public:

  int p;
  double f, logdetw;
  arma::mat X, dX, Xinv;
  arma::uvec indices, lower_indices;

  void param(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices), p, p);

  }

  void F(arguments_optim& x) {

    f = logdetw * (arma::log_det_sympd(X) - arma::sum(arma::log(X.diag())));
    x.f -= f;

  }

  void G(arguments_optim& x) {

    Xinv = arma::inv_sympd(X, arma::inv_opts::allow_approx);

    arma::mat Dinv(p, p, arma::fill::zeros);
    Dinv.diag() = arma::pow(X.diag(), -1.0);

    x.grad.elem(indices) -= logdetw * (Xinv - Dinv);

  }

  void dG(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices), p, p);

    arma::mat dXinv = -Xinv * dX * Xinv;

    arma::mat dDinv(p, p, arma::fill::zeros);
    dDinv.diag() = -dX.diag() % arma::pow(X.diag(), -2.0);

    arma::mat term = dXinv - dDinv;

    x.dgrad.elem(indices) -= logdetw * term;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    doubles[0] =  -f;       // loss actual model
    doubles[1] =  0.00;     // loglik actual model
    doubles[2] =  0.00;     // loglik independence model
    doubles[3] =  0.00;     // loglik saturated model
    doubles[4] =  -f;       // penalty

  }

};


logdetR* choose_logdetR(const Rcpp::List& estimator_setup) {

  logdetR* myestimator = new logdetR();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::uvec lower_indices = estimator_setup["lower_indices"];
  double logdetw = estimator_setup["logdetw"];
  int p = estimator_setup["p"];

  arma::mat X(p, p, arma::fill::zeros);

  myestimator->indices = indices[0];
  myestimator->lower_indices = lower_indices;
  myestimator->X = X;
  myestimator->dX = X;
  myestimator->Xinv = X;
  myestimator->logdetw = logdetw;
  myestimator->p = p;

  return myestimator;

}
