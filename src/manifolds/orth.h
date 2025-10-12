/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Orthogonal manifold:

class orth:public manifolds {

public:

  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);
  arma::vec parameters, dir, dparameters;
  arma::mat g, dg;

  void param(arguments_optim& x) {

    parameters = x.parameters(indices[0]);
    X = arma::reshape(parameters, q, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices[0]), q, q);
    x.rg.elem(indices[0]) = arma::vectorise(X * skew(X.t() * g));

  }

  void hess(arguments_optim& x) {

    dg = arma::reshape(x.dg.elem(indices[0]), q, q);
    dparameters = x.dparameters.elem(indices[0]);
    dX = arma::reshape(dparameters, q, q);
    arma::mat drg = dg - dX * symm(X.t() * g);
    x.dH.elem(indices[0]) = arma::vectorise(X * skew(X.t() * drg));

  }

  void retr(arguments_optim& x) {

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    parameters = arma::vectorise(Q);
    x.parameters(indices[0]) = parameters;

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

orth* choose_orth(Rcpp::List manifold_setup) {

  orth* mymanifold = new orth();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->q = q;

  return mymanifold;

}

arma::mat orth(arma::mat X) {

  arma::mat Q, R;
  arma::qr_econ(Q, R, X);

  return Q;

}

arma::mat rorth(int p, int q) {

  arma::mat X(p, q, arma::fill::randn);
  X = orth(X);

  return X;

}
