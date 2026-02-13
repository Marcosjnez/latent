/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Oblique manifold:

class oblq:public manifolds {

public:

  arma::uvec indices;
  std::size_t q;
  arma::mat X, dX, g, dg;

  void param(arguments_optim& x) {

    X = arma::reshape(x.parameters(indices), q, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), q, q);
    x.rg.elem(indices) = arma::vectorise(g - X * arma::diagmat(X.t() * g));

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), q, q);
    dg = arma::reshape(x.dg.elem(indices), q, q);
    dX = arma::reshape(x.dparameters.elem(indices), q, q);
    x.dH.elem(indices) = arma::vectorise(dg - dX * arma::diagmat(X.t() * g) -
      X * arma::diagmat(X.t() * dg));
    // arma::mat drg = dg - dX * arma::diagmat( X.t() * g) - X * arma::diagmat(dX.t() * g) -
    // X * arma::diagmat(X.t() * dg);
    // dH = drg - X * arma::diagmat(X.t() * drg);

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) = arma::vectorise(X * arma::diagmat(1 / sqrt(arma::sum(X % X, 0))));

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

oblq* choose_oblq(Rcpp::List manifold_setup) {

  oblq* mymanifold = new oblq();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->q = q;

  return mymanifold;

}

arma::mat oblq(arma::mat X) {

  X *= arma::diagmat(1 / sqrt(arma::sum(X % X, 0)));

  return X;

}

arma::mat roblq(int p, int q) {

  arma::mat X(p, q, arma::fill::randn);
  X = oblq(X);

  return X;

}
