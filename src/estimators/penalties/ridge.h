/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 15/02/2026
 */

/*
 * Ridge penalization
 */

class ridge: public estimators {

public:

  double lambda, power, N, loss;
  arma::uvec indices;
  arma::vec beta;

  void param(arguments_optim& x) {

    beta = x.transparameters.elem(indices);

  }

  void F(arguments_optim& x) {

    arma::vec beta_power = arma::pow(arma::abs(beta), power);
    loss = (lambda / (N*power)) * arma::accu(beta_power);
    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::vec grad = lambda * (arma::sign(beta) % arma::pow(arma::abs(beta), power - 1.0));
    x.grad.elem(indices) += grad/N;

  }

  void dG(arguments_optim& x) {
    if(power == 1.0) return;
    arma::vec dbeta = x.dtransparameters.elem(indices);

    arma::vec dgrad = lambda * (power - 1.0) *
      (arma::pow(arma::abs(beta), power - 2.0) % dbeta);

    x.dgrad.elem(indices) += dgrad/N;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = loss;
    names_doubles = {"penalty"};

  }

};

ridge* choose_ridge(const Rcpp::List& estimator_setup) {

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double lambda = estimator_setup["lambda"];
  double power = estimator_setup["power"];
  double N = estimator_setup["N"];

  if(indices.size() != 1 || !std::isfinite(lambda) || lambda < 0 ||
     !std::isfinite(power) || power <= 0 || !std::isfinite(N) || N <= 0) {
    Rcpp::stop("Ridge penalties require non-negative lambda and positive power and N.");
  }
  ridge* myestimator = new ridge();
  myestimator->indices = indices[0];
  myestimator->lambda = lambda;
  myestimator->power = power;
  myestimator->N = N;

  return myestimator;

}
