/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Three-parameter t logarithm likelihood
 */

double log_t3(double x, double mu, double sigma, double nu) {
  double z = (x - mu)/sigma;
  double A = 1.0 + (z*z)/nu;

  return std::lgamma((nu+1.0)/2.0) - std::lgamma(nu/2.0) -
         0.5*std::log(nu*M_PI) - std::log(sigma) - 0.5*(nu+1.0)*std::log(A);
}

double dlog_t3_mu(double x, double mu, double sigma, double nu) {
  double z = (x - mu)/sigma;
  return (nu+1.0)/(nu + z*z) * (z/sigma);
}

double dlog_t3_sigma(double x, double mu, double sigma, double nu) {
  double z = (x - mu)/sigma;
  return -1.0/sigma + (nu+1.0)/(nu + z*z) * (z*z/sigma);
}

double dlog_t3_nu(double x, double mu, double sigma, double nu) {
  double z = (x - mu)/sigma;
  double A = 1.0 + (z*z)/nu;

  return 0.5 * (R::digamma((nu+1.0)/2.0) - R::digamma(nu/2.0) - 1.0/nu) -
         0.5*std::log(A) + 0.5*(nu+1.0) * (z*z)/(nu*(nu+z*z));
}

class t3_loglik: public estimators {

public:

  int p, q;
  arma::uvec indices;
  arma::mat X;
  arma::mat mu, sigma, nu;
  arma::mat log_dt;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    mu    = arma::reshape(pars.subvec(0, p*q-1), p, q);
    sigma = arma::reshape(pars.subvec(p*q, 2*p*q-1), p, q);
    nu    = arma::reshape(pars.subvec(2*p*q, 3*p*q-1), p, q);

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        log_dt(i,j) = log_t3(X(i,j), mu(i,j), sigma(i,j), nu(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_dt);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(3*p*q);
    grad.zeros();
    int k = 0;

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        grad(k) = -dlog_t3_mu(X(i,j), mu(i,j), sigma(i,j), nu(i,j));
        k++;
      }
    }

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        grad(k) =-dlog_t3_sigma(X(i,j), mu(i,j), sigma(i,j), nu(i,j));
        k++;
      }
    }

    for(int i=0;i<p;++i){
      for(int j=0;j<q;++j){
        grad(k) = -dlog_t3_nu(X(i,j), mu(i,j), sigma(i,j), nu(i,j));
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

t3_loglik* choose_t3_loglik(const Rcpp::List& estimator_setup) {
  t3_loglik* myestimator = new t3_loglik();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p = X.n_rows;
  int q = X.n_cols;

  // Allocate parameter matrices
  arma::mat mu(p, q);
  arma::mat sigma(p, q);
  arma::mat nu(p, q);
  arma::mat log_dt(p, q);

  myestimator->indices = indices[0];
  myestimator->X = X;
  myestimator->p = p;
  myestimator->q = q;

  myestimator->mu = mu;
  myestimator->sigma = sigma;
  myestimator->nu = nu;
  myestimator->log_dt = log_dt;

  return myestimator;
}
