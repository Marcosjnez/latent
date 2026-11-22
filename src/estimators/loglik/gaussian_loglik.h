/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Standard normal logarithm likelihood
 */

// const double logSQRT2M_PI = std::log(std::sqrt(2 * M_PI));


class gaussian_loglik: public estimators {

public:

  double N, loss, normalizer;
  arma::uvec indices;
  arma::vec means, inverse_sds, z;

  void param(arguments_optim& x) {
    z = (x.transparameters.elem(indices) - means) % inverse_sds;
  }

  void F(arguments_optim& x) {
    loss = (0.5*arma::dot(z, z) + normalizer)/N;
    x.f += loss;
  }

  void G(arguments_optim& x) {
    x.grad.elem(indices) += (z % inverse_sds)/N;
  }

  void dG(arguments_optim& x) {
    x.dgrad.elem(indices) += (x.dtransparameters.elem(indices) % arma::square(inverse_sds))/N;
  }

  void outcomes(arguments_optim& x) {
    doubles = {loss};
    names_doubles = {"penalty"};
  }


};

gaussian_loglik* choose_gaussian_loglik(const Rcpp::List& estimator_setup) {

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double N = estimator_setup["N"];
  arma::mat means = estimator_setup["means"];
  arma::mat sds = estimator_setup["sds"];
  if(indices.size() != 1 || indices[0].n_elem != means.n_elem ||
     sds.n_rows != means.n_rows || sds.n_cols != means.n_cols ||
     !means.is_finite() || !sds.is_finite() || arma::any(arma::vectorise(sds) <= 0) ||
     !std::isfinite(N) || N <= 0) {
    Rcpp::stop("Gaussian coefficient penalties require compatible means, positive SDs and sample size.");
  }

  gaussian_loglik* myestimator = new gaussian_loglik();
  myestimator->indices = indices[0];
  myestimator->N = N;
  myestimator->means = arma::vectorise(means);
  myestimator->inverse_sds = 1.0/arma::vectorise(sds);
  myestimator->normalizer = means.n_elem*logSQRT2M_PI + arma::accu(arma::log(sds));

  return myestimator;

}
