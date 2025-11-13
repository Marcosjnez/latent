/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Crawford-Ferguson family
 */

class cf: public estimators {

public:

  double k; // Provide this
  arma::mat lambda, N, Mm; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, L2N, ML2, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);
  double ff1, ff2;

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
    L2N = L2 * N;
    ML2 = Mm * L2;
    ff1 = (1-k) * arma::accu(L2 % L2N) / 4;
    ff2 = k * arma::accu(L2 % ML2) / 4;

  }

  void F(arguments_optim& x) {

    f = ff1 + ff2;

  }

  void G(arguments_optim& x) {

    arma::mat f1 = (1-k) * L % L2N;
    arma::mat f2 = k * L % ML2;
    gL = f1 + f2;

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

      arma::mat dL2 = 2 * dL % L;
      arma::mat df1 = (1-k) * dL % L2N + (1-k) * L % (dL2 * N);
      arma::mat df2 = k * dL % ML2 + k * L % (Mm * dL2);
      dgL = df1 + df2;

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      arma::mat dL2 = 2 * dL % L;
      arma::mat df1 = (1-k) * dL % L2N + (1-k) * L % (dL2 * N);
      arma::mat df2 = k * dL % ML2 + k * L % (Mm * dL2);
      dgL = df1 + df2;

      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void H(arguments_optim& x) {

    // Rcpp::stop("H not available");
    hess.set_size(parameters.n_elem, parameters.n_elem); hess.zeros();

  }

  void outcomes(arguments_optim& x) {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    arma::colvec L_vector = arma::vectorise(L);
    int p = lambda.n_rows;
    arma::mat Ip(p, p, arma::fill::eye);
    arma::mat c1 = arma::kron(N.t(), Ip) * arma::diagmat(2*L_vector);
    c1.each_col() %= L_vector;
    arma::mat c2 = arma::diagmat(arma::vectorise(L2 * N));
    arma::mat gf1 = (1-k) * (c1 + c2);

    arma::mat Iq(q, q, arma::fill::eye);
    arma::mat c3 = arma::kron(Iq, Mm) * arma::diagmat(2*L_vector);
    c3.each_col() %= L_vector;
    arma::mat c4 = arma::diagmat(arma::vectorise(Mm * L2));
    arma::mat gf2 = k * (c3 + c4);
    modhessian = gf1 + gf2;

  }

};

cf* choose_cf(const Rcpp::List& estimator_setup) {

 cf* myestimator = new cf();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  double k = estimator_setup["k"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int p = lambda.n_rows;
  int q = lambda.n_cols;
  arma::mat M(p, p, arma::fill::ones);
  M.diag(0).zeros();
  arma::mat N(q, q, arma::fill::ones);
  N.diag(0).zeros();

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->q = q;
  myestimator->Mm = M;
  myestimator->N = N;
  myestimator->k = k;

  return myestimator;

}
