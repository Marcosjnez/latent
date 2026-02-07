/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Varimax
 */

class varimax: public estimators {

public:

  arma::mat lambda, dlambda, L2, HL2, Hh, hL;
  int p, q;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices[0]), p, q);

    L2 = lambda % lambda;
    HL2 = Hh * L2;

  }

  void F(arguments_optim& x) {

    f = -arma::accu(HL2 % HL2) / 4;
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda = -lambda % HL2;
    x.grad.elem(indices[0]) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices[0]), p, q);
    arma::mat dL2 = 2 * dlambda % lambda;

    arma::mat ddf_dlambda = -dlambda % HL2 - lambda % (Hh * dL2);
    x.dgrad.elem(indices[0]) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

varimax* choose_varimax(const Rcpp::List& estimator_setup) {

  varimax* myestimator = new varimax();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];

  arma::vec v(p, arma::fill::ones);
  arma::mat I(p, p, arma::fill::eye);
  arma::mat Hh = I - v * v.t() / (p + 0.0);

  myestimator->indices = indices;
  myestimator->Hh = Hh;
  myestimator->p = p;
  myestimator->q = q;

  return myestimator;

}
