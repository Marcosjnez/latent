/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Linear CLF
 */

class lclf: public estimators {

public:

  int p, q;
  double k, epsilon, a, b, f1, f2;
  arma::uvec indices_lambda, lower, larger;
  arma::mat lambda, dlambda, L2, N, Mm, L2N, ML2, hL;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);

    arma::mat absL = arma::abs(lambda);
    lower = arma::find(absL <= epsilon);
    larger = arma::find(absL > epsilon);
    f1 = arma::accu(a + b*absL.elem(lower) % absL.elem(lower));
    f2 = arma::accu(absL.elem(larger));

  }

  void F(arguments_optim& x) {

    f = f1 + f2;
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda(p, q, arma::fill::zeros);
    df_dlambda.elem(lower) = 2*b*lambda.elem(lower);
    df_dlambda.elem(larger) = arma::sign(lambda.elem(larger));

    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    arma::mat ddf_dlambda(p, q, arma::fill::zeros);

    ddf_dlambda.elem(lower) = 2*b*dlambda.elem(lower);
    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

lclf* choose_lclf(const Rcpp::List& estimator_setup) {

  lclf* myestimator = new lclf();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  double epsilon = estimator_setup["epsilon"];

  double b = 1 / (2*epsilon);
  double a = epsilon - b*epsilon*epsilon;

  myestimator->indices_lambda = indices[0];
  myestimator->p = p;
  myestimator->q = q;
  myestimator->epsilon = epsilon;
  myestimator->a = a;
  myestimator->b = b;

  return myestimator;

}
