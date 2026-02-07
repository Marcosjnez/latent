/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * target
 */

class target: public estimators {

public:

  int p, q;
  arma::mat lambda, dlambda, target, weight, f1, weight2;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices[0]), p, q);
    f1 = weight % (lambda - target);

  }

  void F(arguments_optim& x) {

    f = 0.5*arma::accu(f1 % f1);
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda = weight % f1;

    x.grad.elem(indices[0]) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices[0]), p, q);
    arma::mat ddf_dlambda = weight2 % dlambda;

    x.dgrad.elem(indices[0]) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

target* choose_target(const Rcpp::List& estimator_setup) {

  target* myestimator = new target();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat target = estimator_setup["target"];
  arma::mat weight = estimator_setup["weight"];

  arma::mat weight2 = weight % weight;
  int p = target.n_rows;
  int q = target.n_cols;

  myestimator->indices = indices;
  myestimator->target = target;
  myestimator->weight = weight;
  myestimator->weight2 = weight2;
  myestimator->p = p;
  myestimator->q = q;

  return myestimator;

}
