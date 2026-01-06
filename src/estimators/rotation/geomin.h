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

  double q2, epsilon; // Provide epsilon
  arma::vec term;
  arma::mat lambda; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, LoL2, hL;
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

    q2 = 2/(q + 0.0);
    L2 = L % L;
    L2 += epsilon;
    term = arma::trunc_exp(arma::sum(arma::trunc_log(L2), 1) / q);

  }

  void F(arguments_optim& x) {

    f = arma::accu(term);

  }

  void G(arguments_optim& x) {

    LoL2 = L / L2;
    gL = LoL2 * q2;
    gL.each_col() %= term;

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

      arma::mat c1 = (epsilon - L % L) / (L2 % L2) % dL;
      c1.each_col() %= term;
      arma::mat c2 = LoL2;
      arma::vec term2 = q2 * term % arma::sum(LoL2 % dL, 1);
      c2.each_col() %= term2;

      dgL = q2 * (c1 + c2);
      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      arma::mat c1 = (epsilon - L % L) / (L2 % L2) % dL;
      c1.each_col() %= term;
      arma::mat c2 = LoL2;
      arma::vec term2 = q2 * term % arma::sum(LoL2 % dL, 1);
      c2.each_col() %= term2;

      dgL = q2 * (c1 + c2);
      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void outcomes(arguments_optim& x) {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    arma::mat cx = q2 * LoL2;

    arma::mat c1 = q2*(arma::vectorise(L2) - arma::vectorise(2*L % L)) /
      arma::vectorise(L2 % L2);
    arma::mat gcx = arma::diagmat(c1);
    arma::mat c2 = (1/L2) % (2*L) / q;
    c2.each_col() %= term;
    arma::mat gterm = cbind_diag(c2);
    arma::mat v = gterm.t() * cx;
    modhessian = gcx;
    arma::mat term2 = term;
    for(int i=0; i < (q-1); ++i) term2 = arma::join_cols(term2, term);
    modhessian.each_col() %= term2;
    modhessian += kdiag(v);

  }

};

geomin* choose_geomin(const Rcpp::List& estimator_setup) {

  geomin* myestimator = new geomin();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  double epsilon = estimator_setup["epsilon"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int p = lambda.n_rows;
  int q = lambda.n_cols;

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->epsilon = epsilon;

  return myestimator;

}
