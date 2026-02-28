/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Crawford-Ferguson family
 */

class cf: public estimators {

public:

  int p, q;
  double k, ff1, ff2;
  arma::uvec indices_lambda;
  arma::mat lambda, dlambda, L2, N, Mm, L2N, ML2, hL;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);

    L2 = lambda % lambda;
    L2N = L2 * N;
    ML2 = Mm * L2;
    ff1 = (1-k) * arma::accu(L2 % L2N) / 4;
    ff2 = k * arma::accu(L2 % ML2) / 4;

  }

  void F(arguments_optim& x) {

    f = ff1 + ff2;
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat f1 = (1-k) * lambda % L2N;
    arma::mat f2 = k * lambda % ML2;
    arma::mat df_dlambda = f1 + f2;

    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);

    arma::mat dL2 = 2 * dlambda % lambda;
    arma::mat df1 = (1-k) * dlambda % L2N + (1-k) * lambda % (dL2 * N);
    arma::mat df2 = k * dlambda % ML2 + k * lambda % (Mm * dL2);
    arma::mat ddf_dlambda = df1 + df2;

    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

cf* choose_cf(const Rcpp::List& estimator_setup) {

 cf* myestimator = new cf();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double k = estimator_setup["k"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];

  arma::mat M(p, p, arma::fill::ones);
  M.diag(0).zeros();
  arma::mat N(q, q, arma::fill::ones);
  N.diag(0).zeros();

  myestimator->indices_lambda = indices[0];
  myestimator->p = p;
  myestimator->q = q;
  myestimator->Mm = M;
  myestimator->N = N;
  myestimator->k = k;

  return myestimator;

}
