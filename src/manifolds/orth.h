/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/02/2025
 */

// Orthogonal manifold:

class orth:public manifolds {

public:

  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param() {

    X = arma::reshape(parameters, q, q);

  }

  void proj() {

    g.reshape(q, q);
    rg = X * skew(X.t() * g);

  }

  void hess() {

    g.reshape(q, q);
    dg.reshape(q, q);
    dX = arma::reshape(dparameters, q, q);
    arma::mat drg = dg - dX * symm(X.t() * g);
    dH = X * skew(X.t() * drg);

  }

  void retr() {

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    parameters = arma::vectorise(Q);

  }

  void dconstraints() {

  }

  void outcomes() {

  }

};

orth* choose_orth(Rcpp::List manifold_setup) {

  orth* mymanifold = new orth();

  // Provide these:
  arma::uvec indices = manifold_setup["indices"];
  std::size_t q = manifold_setup["q"];
  mymanifold->indices = indices;
  mymanifold->q = q;

  return mymanifold;

}
