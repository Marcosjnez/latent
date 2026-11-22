/*
 * Author: Vithor R. Franco
 * Modified by: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 21/02/2026
 */

/*
 * Poisson logarithm likelihood
 */

class poisson_loglik: public estimators {

public:

  arma::uvec indices;
  arma::vec counts, lambdas;
  double log_factorials, loss;

  void param(arguments_optim& x) override {
    lambdas = x.transparameters.elem(indices);
  }

  void F(arguments_optim& x) override {
    loss = log_factorials;
    for(arma::uword i = 0; i < counts.n_elem; ++i) {
      const double rate = lambdas(i);
      if(!std::isfinite(rate) || rate < 0.0 || (rate == 0.0 && counts(i) > 0.0)) {
        loss = arma::datum::inf;
        break;
      }
      loss += counts(i) > 0.0 ? rate - counts(i)*std::log(rate) : rate;
    }
    x.f += loss;
  }

  void G(arguments_optim& x) override {
    for(arma::uword i = 0; i < counts.n_elem; ++i) {
      double score = 1.0;
      if(counts(i) > 0.0) score -= counts(i) / lambdas(i);
      x.grad(indices(i)) += score;
    }
  }

  void dG(arguments_optim& x) override {
    for(arma::uword i = 0; i < counts.n_elem; ++i) {
      if(counts(i) > 0.0) {
        x.dgrad(indices(i)) += (counts(i) / lambdas(i)) *
          (x.dtransparameters(indices(i)) / lambdas(i));
      }
    }
  }

  void outcomes(arguments_optim& x) override {
    doubles = {loss, -loss};
    names_doubles = {"loss", "loglik"};
  }
};

poisson_loglik* choose_poisson_loglik(const Rcpp::List& estimator_setup) {
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat observations = estimator_setup["X"];
  if(indices.size() != 1 || indices[0].n_elem != observations.n_elem) {
    Rcpp::stop("poisson_loglik requires one rate per observation.");
  }
  if(!observations.is_finite() || arma::any(arma::vectorise(observations < 0)) ||
     arma::any(arma::vectorise(observations != arma::floor(observations)))) {
    Rcpp::stop("Poisson observations must be non-negative integer counts.");
  }

  poisson_loglik* myestimator = new poisson_loglik();
  myestimator->indices = indices[0];
  myestimator->counts = arma::vectorise(observations);
  myestimator->log_factorials = arma::accu(arma::lgamma(observations + 1.0));
  return myestimator;
}
