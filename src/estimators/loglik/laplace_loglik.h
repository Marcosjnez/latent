/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Laplace (double exponential) logarithm likelihood
 */

double log_laplace(double x, double mu, double sigma) {
  return -std::log(2.0 * sigma) - std::fabs(x - mu) / sigma;
}

double dlog_laplace_mu(double x, double mu, double sigma) {
  double diff = x - mu;
  if(diff > 0.0) return  1.0 / sigma;
  if(diff < 0.0) return -1.0 / sigma;
  return 0.0;
}

double dlog_laplace_sigma(double x, double mu, double sigma) {
  return -1.0 / sigma + std::fabs(x - mu) / (sigma * sigma);
}

class laplace_loglik: public estimators {

public:

  int p, q;
  arma::uvec indices;
  arma::mat X;
  arma::mat mu, sigma;
  arma::mat log_dl;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    mu    = arma::reshape(pars.subvec(0, p*q-1), p, q);
    sigma = arma::reshape(pars.subvec(p*q, 2*p*q-1), p, q);

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        log_dl(i,j) = log_laplace(X(i,j), mu(i,j), sigma(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_dl);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(2*p*q);
    grad.zeros();
    int k = 0;

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        grad(k) = -dlog_laplace_mu(X(i,j), mu(i,j), sigma(i,j));
        k++;
      }
    }

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        grad(k) = -dlog_laplace_sigma(X(i,j), mu(i,j), sigma(i,j));
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

laplace_loglik* choose_laplace_loglik(const Rcpp::List& estimator_setup) {

  laplace_loglik* myestimator = new laplace_loglik();

  arma::uvec indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p = X.n_rows;
  int q = X.n_cols;

  arma::mat mu(p, q);
  arma::mat sigma(p, q);
  arma::mat log_dl(p, q);

  myestimator->indices = indices;
  myestimator->X = X;
  myestimator->p = p;
  myestimator->q = q;

  myestimator->mu = mu;
  myestimator->sigma = sigma;
  myestimator->log_dl = log_dl;

  return myestimator;
}
