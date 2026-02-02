/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 02/02/2026
 */

/*
 * Standard normal logarithm likelihood
 */

// const double logSQRT2M_PI = std::log(std::sqrt(2 * M_PI));

double log_stdnorm(double x) {
  // loglik of the std normal
  return -0.5*x*x - logSQRT2M_PI;
}

double dlog_stdnorm(double x, double dx = 1.00) {
  // Differential of the loglik of the std normal
  return -x*dx;
}

double ddlog_stdnorm(double dx = 1.00) {
  // Differential of the derivative of the loglik of the std normal
  return -dx;
}

class gaussian_loglik: public estimators {

public:

  arma::mat X, means, sds, Z, log_dnorms, df_DZ;
  double alpha, N;
  int p, q;

  void param(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices[0]), p, q);

    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        Z(i,j) = (X(i,j) - means(i,j)) / sds(i,j);
        log_dnorms(i,j) = log_stdnorm(Z(i,j));
      }
    }

  }

  void F(arguments_optim& x) {

    f = arma::accu(log_dnorms)/N;
    x.f += -f;

  }

  void G(arguments_optim& x) {

    arma::mat df_dZ(p, q);
    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        df_dZ(i, j) = -dlog_stdnorm(Z(i, j), 1.00)/sds(i,j);
      }
    }

    x.grad.elem(indices[0]) += arma::vectorise(df_dZ)/N;

  }

  void dG(arguments_optim& x) {

    arma::mat dZ = reshape(x.dtransparameters(indices[0]), p, q);

    arma::mat ddf_dZ(p, q);
    for(int i=0; i < p; ++i) {
      for(int j=0; j < q; ++j) {
        ddf_dZ(i, j) = dZ(i,j)/sds(i,j)/sds(i,j);
      }
    }

    x.dgrad.elem(indices[0]) += arma::vectorise(ddf_dZ)/N;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] =  f;

  };

};

gaussian_loglik* choose_gaussian_loglik(const Rcpp::List& estimator_setup) {

  gaussian_loglik* myestimator = new gaussian_loglik();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double alpha = estimator_setup["alpha"];
  double N = estimator_setup["N"];
  arma::mat means = estimator_setup["means"];
  arma::mat sds = estimator_setup["sds"];
  int p = means.n_rows;
  int q = means.n_cols;

  arma::mat log_dnorms(p, q);
  arma::mat Z(p, q);

  myestimator->indices = indices;
  myestimator->alpha = alpha;
  myestimator->N = N;
  myestimator->means = means;
  myestimator->sds = sds;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->log_dnorms = log_dnorms;
  myestimator->Z = Z;

  return myestimator;

}
