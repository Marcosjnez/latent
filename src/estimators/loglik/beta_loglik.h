/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Mean_precision Beta logarithm likelihood
 */

double log_beta_mean(double x, double mu, double phi) {
  double alpha = mu * phi;
  double beta  = (1.0 - mu) * phi;

  return std::lgamma(phi) - std::lgamma(alpha) - std::lgamma(beta) +
    (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log(1.0 - x);
}

double dlog_beta_mu(double x, double mu, double phi) {
  double alpha = mu * phi;
  double beta  = (1.0 - mu) * phi;

  return phi * (R::digamma(beta) - R::digamma(alpha) + std::log(x/(1.0 - x)));
}

double dlog_beta_phi(double x, double mu, double phi) {
  double alpha = mu * phi;
  double beta  = (1.0 - mu) * phi;

  return R::digamma(phi) - mu * R::digamma(alpha) -
    (1.0 - mu) * R::digamma(beta) + mu * std::log(x) +
    (1.0 - mu) * std::log(1.0 - x);
}

class beta_loglik: public estimators {

public:

  int p_dim, q_dim;
  arma::uvec indices;
  arma::mat X;
  arma::mat mu, phi;
  arma::mat log_db;

  void param(arguments_optim& x) {
    arma::vec pars = x.transparameters(indices);
    arma::mat mu = arma::reshape(pars.subvec(0, p_dim*q_dim-1), p_dim, q_dim);
    arma::mat phi = arma::reshape(pars.subvec(p_dim*q_dim, 2*p_dim*q_dim-1),
                                  p_dim, q_dim);

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        log_db(i,j) = log_beta_mean(X(i,j), mu(i,j), phi(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_db);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(2*p_dim*q_dim);
    grad.zeros();
    int k = 0;

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        double dl_dmu = dlog_beta_mu(X(i,j), mu(i,j), phi(i,j));
        grad(k) = -dl_dmu * mu(i,j) * (1.0 - mu(i,j));
        k++;
      }
    }

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        double dl_dphi = dlog_beta_phi(X(i,j), mu(i,j), phi(i,j));
        grad(k) = -dl_dphi * phi(i,j);
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

beta_mean_loglik* choose_beta_mean_loglik(const Rcpp::List& estimator_setup) {
  beta_mean_loglik* myestimator = new beta_mean_loglik();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];

  int p_dim = X.n_rows;
  int q_dim = X.n_cols;

  arma::mat mu(p_dim, q_dim);
  arma::mat phi(p_dim, q_dim);
  arma::mat log_db(p_dim, q_dim);

  myestimator->indices = indices[0];
  myestimator->X = X;
  myestimator->p_dim = p_dim;
  myestimator->q_dim = q_dim;

  myestimator->mu = mu;
  myestimator->phi = phi;
  myestimator->log_db = log_db;

  return myestimator;
}
