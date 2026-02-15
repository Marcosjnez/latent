/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * target
 */

class target: public estimators {

public:

  int p, q;
  arma::uvec indices_lambda;
  arma::mat lambda, dlambda, target, weight, f1, weight2;

  void param(arguments_optim& x) {

    // Parameterization
    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);
    f1 = weight % (lambda - target);

  }

  void F(arguments_optim& x) {

    // Compute loss
    f = 0.5*arma::accu(f1 % f1);
    x.f += f;

  }

  void G(arguments_optim& x) {

    // Compute gradient
    arma::mat df_dlambda = weight % f1;

    x.grad(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    arma::mat ddf_dlambda = weight2 % dlambda;

    x.dgrad(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[1] =  0.00;

  }

};

target* choose_target(const Rcpp::List& estimator_setup) {

  target* myestimator = new target();

  arma::uvec indices_lambda = estimator_setup["indices_lambda"];
  arma::mat target = estimator_setup["target"];
  arma::mat weight = estimator_setup["weight"];

  arma::mat weight2 = weight % weight;
  int p = target.n_rows;
  int q = target.n_cols;

  myestimator->indices_lambda = indices_lambda;
  myestimator->target = target;
  myestimator->weight = weight;
  myestimator->weight2 = weight2;
  myestimator->p = p;
  myestimator->q = q;

  return myestimator;

}
