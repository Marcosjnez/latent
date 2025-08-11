/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/02/2025
 */

/*
 * Varimax
 */

class varimax: public estimators {

public:

  arma::mat lambda, Hh; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, HL2, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param() {

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
    HL2 = Hh * L2;


  }

  void F() {

    f = -arma::trace(HL2.t() * HL2) / 4;

  }

  void G() {

    gL = -L % HL2;

    if(orth) {

      g = lambda.t() * gL;

    } else {

      g = - Inv_X.t() * gL.t() * L;

    }

  }

  void dG() {

    dX = arma::reshape(dparameters, q, q);
    g = arma::reshape(g, q, q);

    if(orth) {

      dL = lambda * dX;

      arma::mat dL2 = 2 * dL % L;
      dgL = -dL % HL2 - L % (Hh * dL2);

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      arma::mat dL2 = 2 * dL % L;
      dgL = -dL % HL2 - L % (Hh * dL2);

      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void H() {

    // Rcpp::stop("H not available");
    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() {}

  void M() {}

  void outcomes() {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    arma::mat Iq(q, q, arma::fill::eye);
    arma::mat c1 = arma::diagmat(arma::vectorise(HL2));
    arma::colvec L_vector = arma::vectorise(L);
    arma::mat c2 = arma::kron(Iq, Hh) * arma::diagmat(2*L_vector);
    c2.each_col() %= L_vector;
    modhessian = -(c1 + c2);

  }

};

varimax* choose_varimax(const Rcpp::List& estimator_setup) {

  varimax* myestimator = new varimax();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int p = lambda.n_rows;
  int q = lambda.n_cols;

  arma::vec v(p, arma::fill::ones);
  arma::mat I(p, p, arma::fill::eye);
  arma::mat H = I - v * v.t() / (p + 0.0);

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->Hh = H;
  myestimator->q = q;

  return myestimator;

}
