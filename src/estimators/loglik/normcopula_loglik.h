/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 09/03/2026
 */

/*
 * Gaussian copula logarithm likelihood
 */

double log_normcopula(const arma::rowvec& z, const arma::mat& Rinv,
                      double logdetR) {
  arma::rowvec tmp = z * (Rinv - arma::eye(Rinv.n_rows, Rinv.n_cols));
  return -0.5 * logdetR - 0.5 * arma::as_scalar(tmp * z.t());
}

class normcopula_loglik: public estimators {

public:

  int p, q;
  double N;
  arma::uvec indices;
  arma::mat U;
  arma::mat Z;
  arma::mat R;
  arma::mat Rinv;
  arma::vec log_copula;
  double logdetR;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    R = arma::reshape(pars, q, q);
    R = arma::symmatu(R);
    arma::log_det(logdetR, std::ignore, R);
    Rinv = arma::inv_sympd(R);

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        Z(i,j) = R::qnorm(U(i,j), 0.0, 1.0, 1, 0);
      }
      log_copula(i) = log_normcopula(Z.row(i), Rinv, logdetR);
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_copula);
    x.f += -f/N;
  }

  void outcomes(arguments_optim& x) {
    doubles.resize(1);
    doubles[0] = f;
  }
};

normcopula_loglik* choose_normcopula_loglik(const Rcpp::List& estimator_setup) {

  normcopula_loglik* myestimator = new normcopula_loglik();

  arma::uvec indices = estimator_setup["indices"];
  arma::mat U = estimator_setup["U"];
  double N = estimator_setup["N"];

  int p = U.n_rows;
  int q = U.n_cols;

  arma::mat Z(p,q);
  arma::vec log_copula(p);

  myestimator->indices = indices;
  myestimator->U = U;
  myestimator->Z = Z;
  myestimator->log_copula = log_copula;

  myestimator->N = N;
  myestimator->p = p;
  myestimator->q = q;

  return myestimator;
}
