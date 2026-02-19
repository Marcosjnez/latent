/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Poisson logarithm likelihood
 */

// const double logSQRT2M_PI = std::log(std::sqrt(2 * M_PI));

double log_poisson(double x, double lambda) {
  return x * std::log(lambda) - lambda - std::lgamma(x + 1.0);
}

double dlog_poisson(double x, double lambda) {
  return x / lambda - 1.0;
}

double ddlog_poisson(double x, double lambda) {
  return -x / (lambda * lambda);
}

class poisson_loglik: public estimators {
  
public:
  
  int p, q;
  double alpha, N;
  arma::uvec indices;
  arma::mat X, lambdas, log_dpois;
  
  void param(arguments_optim& x) {
    lambdas = arma::reshape(x.transparameters(indices), p, q);
    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        log_dpois(i,j) = log_poisson(X(i,j), lambdas(i,j));
      }
    }
  }
  
  void F(arguments_optim& x) {
    f = arma::accu(log_dpois);
    x.f += -f/N;
  }
  
  void G(arguments_optim& x) {
    arma::mat df_dlambda(p, q);
    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        df_dlambda(i,j) = -dlog_poisson(X(i,j), lambdas(i,j));
      }
    }
    x.grad.elem(indices) += arma::vectorise(df_dlambda)/N;
  }
  
  void dG(arguments_optim& x) {
    arma::mat dlambda = reshape(x.dtransparameters(indices), p, q);
    arma::mat ddf_dlambda(p, q);
    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        ddf_dlambda(i,j) = 
          -ddlog_poisson(X(i,j), lambdas(i,j)) * dlambda(i,j);
      }
    }
    x.dgrad.elem(indices) += arma::vectorise(ddf_dlambda)/N;
  }
  
  void outcomes(arguments_optim& x) {
    doubles.resize(1);
    doubles[0] = f;
  };
};

poisson_loglik* choose_poisson_loglik(const Rcpp::List& estimator_setup) {
  
  poisson_loglik* myestimator = new poisson_loglik();
  
  arma::uvec indices = estimator_setup["indices"];
  double alpha = estimator_setup["alpha"];
  double N = estimator_setup["N"];
  arma::mat X = estimator_setup["X"];
  
  int p = X.n_rows;
  int q = X.n_cols;
  
  arma::mat log_dpois(p, q);
  arma::mat lambdas(p, q);
  
  myestimator->indices = indices;
  myestimator->alpha = alpha;
  myestimator->N = N;
  myestimator->X = X;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->log_dpois = log_dpois;
  myestimator->lambdas = lambdas;
  
  return myestimator;
}
