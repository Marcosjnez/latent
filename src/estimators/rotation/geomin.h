/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Geomin
 */

class geomin: public estimators {

public:

  int p, q;
  arma::mat lambda, dlambda, L2, LoL2;
  double q2, epsilon;
  arma::vec term;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices[0]), p, q);

    L2 = lambda % lambda;
    L2 += epsilon;
    term = arma::trunc_exp(arma::sum(arma::trunc_log(L2), 1) / q);

  }

  void F(arguments_optim& x) {

    f = arma::accu(term);
    x.f += f;

  }

  void G(arguments_optim& x) {

    LoL2 = lambda / L2;
    arma::mat df_dlambda = LoL2 * q2;
    df_dlambda.each_col() %= term;

    x.grad.elem(indices[0]) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices[0]), p, q);

    arma::mat c1 = (epsilon - lambda % lambda) / (L2 % L2) % dlambda;
    c1.each_col() %= term;
    arma::mat c2 = LoL2;
    arma::vec term2 = q2 * term % arma::sum(LoL2 % dlambda, 1);
    c2.each_col() %= term2;

    arma::mat ddf_dlambda = q2 * (c1 + c2);

    x.dgrad.elem(indices[0]) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

geomin* choose_geomin(const Rcpp::List& estimator_setup) {

  geomin* myestimator = new geomin();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  double epsilon = estimator_setup["epsilon"];

  double q2 = 2/(q + 0.0);

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->q2 = q2;
  myestimator->epsilon = epsilon;

  return myestimator;

}
