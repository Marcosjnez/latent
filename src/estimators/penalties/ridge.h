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

  double lambda, power, N;
  arma::uvec indices;
  arma::vec beta;

  void param(arguments_optim& x) {

    beta = x.transparameters.elem(indices);

  }

  void F(arguments_optim& x) {

    arma::vec beta_power = arma::pow(arma::abs(beta), power);
    f = -(lambda / power) * arma::accu(beta_power);
    x.f += f/N;

  }

  void G(arguments_optim& x) {

    arma::vec grad = lambda * (arma::sign(beta) % arma::pow(arma::abs(beta), power - 1.0));
    x.grad.elem(indices) += -grad/N;

  }

  void dG(arguments_optim& x) {
    arma::vec dbeta = x.dtransparameters.elem(indices);

    arma::vec dgrad = lambda * (power - 1.0) *
      (arma::pow(arma::abs(beta), power - 2.0) % dbeta);

    x.dgrad.elem(indices) += -dgrad/N;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = f;

  }

};

ridge* choose_ridge(const Rcpp::List& estimator_setup) {

  ridge* myestimator = new ridge();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double lambda = estimator_setup["lambda"];
  double power = estimator_setup["power"];
  double N = estimator_setup["N"];

  myestimator->indices = indices[0];
  myestimator->lambda = lambda;
  myestimator->power = power;
  myestimator->N = N;

  return myestimator;

}
