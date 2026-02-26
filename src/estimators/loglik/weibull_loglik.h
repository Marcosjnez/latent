/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Weibull logarithm likelihood
 */

double log_weibull(double x, double k, double lambda) {
  double logx = std::log(x);
  double loglam = std::log(lambda);
  double z = std::pow(x / lambda, k);

  return std::log(k) - k * loglam + (k - 1.0) * logx - z;
}

double dlog_weibull_lambda(double x, double k, double lambda) {
  double z = std::pow(x / lambda, k);

  return (k / lambda) * (z - 1.0);
}

double dlog_weibull_k(double x, double k, double lambda) {
  double log_ratio = std::log(x / lambda);
  double z = std::pow(x / lambda, k);

  return 1.0 / k - std::log(lambda) + std::log(x) - z * log_ratio;
}

class weibull_loglik: public estimators {

public:

  int p, q;
  arma::uvec indices;
  arma::mat X;
  arma::mat k, lambda;
  arma::mat log_dw;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    k      = arma::reshape(pars.subvec(0, p*q-1), p, q);
    lambda = arma::reshape(pars.subvec(p*q, 2*p*q-1), p, q);

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        log_dw(i,j) = log_weibull(X(i,j), k(i,j), lambda(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_dw);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(2*p*q);
    grad.zeros();
    int idx = 0;

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        double dl_dk = dlog_weibull_k(X(i,j), k(i,j), lambda(i,j));
        grad(idx) = -dl_dk;
        idx++;
      }
    }

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        double dl_dlambda = dlog_weibull_lambda(X(i,j), k(i,j), lambda(i,j));
        grad(idx) = -dl_dlambda;
        idx++;
      }
    }

    x.grad.elem(indices) += grad;
  }

  void outcomes(arguments_optim& x) {
    doubles.resize(1);
    doubles[0] = f;
  }
};

weibull_loglik* choose_weibull_loglik(const Rcpp::List& estimator_setup) {
  weibull_loglik* myestimator = new weibull_loglik();

  arma::uvec indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p = X.n_rows;
  int q = X.n_cols;

  arma::mat k(p, q);
  arma::mat lambda(p, q);
  arma::mat log_dw(p, q);

  myestimator->indices = indices;
  myestimator->X = X;
  myestimator->p = p;
  myestimator->q = q;

  myestimator->k = k;
  myestimator->lambda = lambda;
  myestimator->log_dw = log_dw;

  return myestimator;
}
