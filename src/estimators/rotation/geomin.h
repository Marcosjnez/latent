/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Geomin
 */

class geomin: public estimators {

public:

  int p, q;
  double q2, epsilon;
  arma::uvec indices_lambda;
  arma::vec term;
  arma::mat lambda, dlambda, L2, LoL2;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);

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

    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);

    arma::mat c1 = (epsilon - lambda % lambda) / (L2 % L2) % dlambda;
    c1.each_col() %= term;
    arma::mat c2 = LoL2;
    arma::vec term2 = q2 * term % arma::sum(LoL2 % dlambda, 1);
    c2.each_col() %= term2;

    arma::mat ddf_dlambda = q2 * (c1 + c2);

    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

geomin* choose_geomin(const Rcpp::List& estimator_setup) {

  geomin* myestimator = new geomin();

  arma::uvec indices_lambda = estimator_setup["indices_lambda"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  double epsilon = estimator_setup["epsilon"];

  double q2 = 2/(q + 0.0);

  myestimator->indices_lambda = indices_lambda;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->q2 = q2;
  myestimator->epsilon = epsilon;

  return myestimator;

}
