/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Exponential (rate parametrization) logarithm likelihood
 */

double log_exponential(double x, double lambda) {
  return std::log(lambda) - lambda * x;
}

double dlog_exponential_lambda(double x, double lambda) {
  return 1.0 / lambda - x;
}

class exponential_loglik: public estimators {

public:

  int p, q;
  arma::uvec indices;
  arma::mat X;
  arma::mat lambda;
  arma::mat log_de;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    lambda = arma::reshape(pars.subvec(0, p*q-1), p, q);

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        log_de(i,j) = log_exponential(X(i,j), lambda(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_de);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(p*q);
    grad.zeros();
    int k = 0;

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        double dl_dlambda = dlog_exponential_lambda(X(i,j), lambda(i,j));
        grad(k) = -dl_dlambda;
        k++;
      }
    }

    x.grad.elem(indices) += grad;
  }

  void outcomes(arguments_optim& x) {
    doubles.resize(1);
    doubles[0] = f;
  }
};

exponential_loglik* choose_exponential_loglik(const Rcpp::List& estimator_setup) {
  exponential_loglik* myestimator = new exponential_loglik();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p = X.n_rows;
  int q = X.n_cols;

  arma::mat lambda(p, q);
  arma::mat log_de(p, q);

  myestimator->indices = indices[0];
  myestimator->X = X;
  myestimator->p = p;
  myestimator->q = q;

  myestimator->lambda = lambda;
  myestimator->log_de = log_de;

  return myestimator;
}
