/*
 * Author: Vithor R. Franco
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/02/2026
 */

/*
 * Binomial logarithm likelihood
 */

double log_binom(double x, double n, double p) {
  return std::lgamma(n + 1.0) - std::lgamma(x + 1.0) -
         std::lgamma(n - x + 1.0) + x * std::log(p) +
         (n - x) * std::log(1.0 - p);
}

double dlog_binom(double x, double n, double p) {
  return x/p - (n - x)/(1.0 - p);
}

class binomial_loglik: public estimators {

public:

  int p_dim, q_dim;
  arma::uvec indices;
  arma::mat X;
  arma::mat Ntrials;
  arma::mat eta;
  arma::mat prob;
  arma::mat log_db;

  void param(arguments_optim& x) {
    eta = arma::reshape(x.transparameters(indices), p_dim, q_dim);
    prob = 1.0 / (1.0 + arma::exp(-eta));

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        log_db(i,j) =
          log_binom(X(i,j), Ntrials(i,j), prob(i,j));
      }
    }
  }

  void F(arguments_optim& x) {
    f = arma::accu(log_db);
    x.f += -f;
  }

  void G(arguments_optim& x) {
    arma::vec grad(p_dim*q_dim);
    grad.zeros();
    int k = 0;

    for(int i=0;i<p_dim;++i){
      for(int j=0;j<q_dim;++j){
        double dl_dp = dlog_binom(X(i,j), Ntrials(i,j), prob(i,j));
        grad(k) = -dl_dp * prob(i,j) * (1.0 - prob(i,j));
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

binomial_loglik* choose_binomial_loglik(const Rcpp::List& estimator_setup) {
  binomial_loglik* myestimator = new binomial_loglik();

  arma::uvec indices = estimator_setup["indices"];
  arma::mat X = estimator_setup["X"];
  arma::mat Ntrials = estimator_setup["Ntrials"];

  int p_dim = X.n_rows;
  int q_dim = X.n_cols;

  arma::mat eta(p_dim, q_dim);
  arma::mat prob(p_dim, q_dim);
  arma::mat log_db(p_dim, q_dim);

  myestimator->indices = indices;
  myestimator->X = X;
  myestimator->Ntrials = Ntrials;
  myestimator->p_dim = p_dim;
  myestimator->q_dim = q_dim;

  myestimator->eta = eta;
  myestimator->prob = prob;
  myestimator->log_db = log_db;

  return myestimator;
}
