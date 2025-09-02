/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Oblimin
 */

class oblimin: public estimators {

public:

  arma::mat lambda, I_gamma_C, N; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, IgCL2N, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param(arguments_optim& x) {

    X = arma::reshape(parameters, q, q);

    if(orth) {

      L = lambda*X;
      Phi = arma::mat(q, q, arma::fill::eye);

    } else {

      Phi = X.t() * X;
      Inv_X = arma::inv(X);
      L = lambda * Inv_X.t();

    }

    L2 = L % L;
    IgCL2N = I_gamma_C * L2 * N;

  }

  void F(arguments_optim& x) {

    f = arma::accu(L2 % IgCL2N) / 4;

  }

  void G(arguments_optim& x) {

    gL = L % IgCL2N;

    if(orth) {

      g = lambda.t() * gL;

    } else {

      g = - Inv_X.t() * gL.t() * L;

    }

  }

  void dG(arguments_optim& x) {

    dX = arma::reshape(dparameters, q, q);
    g = arma::reshape(g, q, q);

    if(orth) {

      dL = lambda * dX;

      dgL = dL % IgCL2N + L % (I_gamma_C * (2*dL % L) * N);

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      dgL = dL % IgCL2N + L % (I_gamma_C * (2*dL % L) * N);

      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void H(arguments_optim& x) {

    // Rcpp::stop("H not available");
    hess.set_size(parameters.n_elem, parameters.n_elem); hess.zeros();

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void outcomes(arguments_optim& x) {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    arma::mat c1 = arma::diagmat(arma::vectorise(IgCL2N));
    arma::mat diagL = arma::diagmat(arma::vectorise(L));
    arma::mat diag2L = 2*diagL;
    arma::mat c2 = diagL * arma::kron(N, I_gamma_C) * diag2L;

    modhessian = c1 + c2;

  }

};

oblimin* choose_oblimin(const Rcpp::List& estimator_setup) {

  oblimin* myestimator = new oblimin();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  double gamma = estimator_setup["gamma"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int p = lambda.n_rows;
  int q = lambda.n_cols;
  arma::mat N(q, q, arma::fill::ones);
  N.diag(0).zeros();
  arma::mat I(p, p, arma::fill::eye), gamma_C(p, p, arma::fill::ones);
  gamma_C *= (gamma/p);
  arma::mat I_gamma_C = (I - gamma_C);

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->q = q;
  myestimator->N = N;
  myestimator->I_gamma_C = I_gamma_C;

  return myestimator;

}
