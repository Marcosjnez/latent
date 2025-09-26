/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/09/2025
 */

/*
 * Logarithm of determinant of a matrix
 */

class logdetmat: public estimators {

public:

  arma::mat X;
  arma::uvec lower_indices;
  double logdetw;
  int p;

  void param(arguments_optim& x) {

    X.elem(lower_indices) = x.transparameters(indices[0]);
    X = arma::symmatl(X);

  }

  void F(arguments_optim& x) {

    // double val;
    // double sign;
    // bool ok = arma::log_det(val, sign, X);
    // f = val;

    f = logdetw*arma::log_det_sympd(X);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    arma::mat Xinv = arma::inv_sympd(X, arma::inv_opts::allow_approx);
    x.grad.elem(indices[0]) -= logdetw*Xinv.elem(lower_indices);

  }

  void dG(arguments_optim& x) {

    // dg.set_size(x.transparameters.n_elem); dg.zeros();
    // x.dgrad.set_size(x.transparameters.n_elem); x.dgrad.zeros();
    // x.dgrad(indices[0]) += arma::vectorise();

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

    // x.uniquenesses = x.R.diag() - arma::diagvec(x.Rhat);
    // x.Rhat.diag() = x.R.diag();

    doubles.resize(1);
    doubles[0] = f;

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
  myestimator->logdetw = logdetw;
  myestimator->p = p;

  return myestimator;

}
