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
  arma::uvec indices_lambda;
  arma::mat lambda, I_gamma_C, N, dlambda, L2, IgCL2N;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);

    L2 = lambda % lambda;
    IgCL2N = I_gamma_C * L2 * N;

  }

  void F(arguments_optim& x) {

    f = arma::accu(L2 % IgCL2N) / 4;
    x.f += f;
    // Rprintf("f = %.6f\n", f);

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda = lambda % IgCL2N;
    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    arma::mat ddf_dlambda = dlambda % IgCL2N +
                            lambda % (I_gamma_C * (2*dlambda % lambda) * N);

    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

oblimin* choose_oblimin(const Rcpp::List& estimator_setup) {

  oblimin* myestimator = new oblimin();

  arma::uvec indices_lambda = estimator_setup["indices_lambda"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  double gamma = estimator_setup["gamma"];

  arma::mat N(q, q, arma::fill::ones);
  N.diag(0).zeros();
  arma::mat I(p, p, arma::fill::eye), gamma_C(p, p, arma::fill::ones);
  gamma_C *= (gamma/p);
  arma::mat I_gamma_C = (I - gamma_C);
  // I_gamma_C and N must be symmetric

  myestimator->indices_lambda = indices_lambda;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->N = N;
  myestimator->I_gamma_C = I_gamma_C;

  return myestimator;

}
