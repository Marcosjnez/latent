/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Oblimin
 */

class oblimin: public estimators {

public:

  int p, q;
  double loss, gamma;
  arma::uvec indices_lambda;
  arma::mat lambda, dlambda, L2, IgCL2N;

  arma::mat other_factors(const arma::mat& L2) {
    arma::mat result = -L2;
    result.each_col() += arma::sum(L2, 1);
    result.each_row() -= gamma * arma::mean(result, 0);
    return result;
  }

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);

    L2 = lambda % lambda;
    IgCL2N = other_factors(L2);

  }

  void F(arguments_optim& x) {

    loss = arma::accu(L2 % IgCL2N) / 4;
    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda = lambda % IgCL2N;
    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    arma::mat ddf_dlambda = dlambda % IgCL2N +
                            lambda % other_factors(2*dlambda % lambda);

    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  loss;
    doubles[1] =  0.00;

  }

};

oblimin* choose_oblimin(const Rcpp::List& estimator_setup) {

  oblimin* myestimator = new oblimin();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  double gamma = estimator_setup["gamma"];

  myestimator->indices_lambda = indices[0];
  myestimator->p = p;
  myestimator->q = q;
  myestimator->gamma = gamma;

  return myestimator;

}
