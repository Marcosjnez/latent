/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Orthogonal manifold:

class orth:public manifolds {

public:

  arma::uvec indices;
  std::size_t p, q;
  arma::mat X, dX, g, dg;

  void param(arguments_optim& x) {

    X = arma::reshape(x.parameters, q, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), q, q);
    x.rg.elem(indices) = arma::vectorise(X * skew(X.t() * g));

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), q, q);
    dg = arma::reshape(x.dg.elem(indices), q, q);
    dX = arma::reshape(x.dparameters.elem(indices), q, q);
    arma::mat drg = dg - dX * symm(X.t() * g);
    x.dH.elem(indices) = arma::vectorise(X * skew(X.t() * drg));

  }

  void retr(arguments_optim& x) {

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    x.parameters(indices) = arma::vectorise(Q);

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

orth* choose_orth(Rcpp::List manifold_setup) {

  orth* mymanifold = new orth();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->p = p;
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
