/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Gamma logarithm likelihood
 */

double log_gamma_dist(double x, double alpha, double beta) {
  return alpha * std::log(beta) - std::lgamma(alpha) +
         (alpha - 1.0) * std::log(x) - beta * x;
}

double dlog_gamma_alpha(double x, double alpha, double beta) {
  return std::log(beta) - R::digamma(alpha) + std::log(x);
}

double dlog_gamma_beta(double x, double alpha, double beta) {
  return alpha / beta - x;
}

class gamma_loglik: public estimators {

public:

  int p_dim, q_dim;
  arma::uvec indices;
  arma::mat X;
  arma::mat alpha, beta;
  arma::mat log_dg;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    arma::mat alpha = arma::reshape(pars.subvec(0, p_dim*q_dim-1),
                                    p_dim, q_dim);
    arma::mat beta = arma::reshape(pars.subvec(p_dim*q_dim,
                                   2*p_dim*q_dim-1), p_dim, q_dim);

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        log_dg(i,j) = log_gamma_dist(X(i,j), alpha(i,j), beta(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_dg);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(2*p_dim*q_dim);
    grad.zeros();
    int k = 0;

    // ---- alpha block ----
    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        double dl_dalpha = dlog_gamma_alpha(X(i,j), alpha(i,j), beta(i,j));
        grad(k) = -dl_dalpha * alpha(i,j);
        k++;
      }
    }

    // ---- beta block ----
    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        double dl_dbeta = dlog_gamma_beta(X(i,j), alpha(i,j), beta(i,j));
        grad(k) = -dl_dbeta * beta(i,j);
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

gamma_loglik* choose_gamma_loglik(const Rcpp::List& estimator_setup) {
  gamma_loglik* myestimator = new gamma_loglik();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p_dim = X.n_rows;
  int q_dim = X.n_cols;

  arma::mat alpha(p_dim, q_dim);
  arma::mat beta(p_dim, q_dim);
  arma::mat log_dg(p_dim, q_dim);

  myestimator->indices = indices[0];
  myestimator->X = X;
  myestimator->p_dim = p_dim;
  myestimator->q_dim = q_dim;

  myestimator->alpha = alpha;
  myestimator->beta = beta;
  myestimator->log_dg = log_dg;

  return myestimator;
}
