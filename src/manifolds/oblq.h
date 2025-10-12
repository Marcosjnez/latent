/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Oblique manifold:

class oblq:public manifolds {

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
    x.rg.elem(indices[0]) = arma::vectorise(g - X * arma::diagmat(X.t() * g));

  }

  void hess(arguments_optim& x) {

    dg = arma::reshape(x.dg.elem(indices[0]), q, q);
    dparameters = x.dparameters.elem(indices[0]);
    dX = arma::reshape(dparameters, q, q);
    x.dH.elem(indices[0]) = arma::vectorise(dg - dX * arma::diagmat(X.t() * g) -
      X * arma::diagmat(X.t() * dg));
    // arma::mat drg = dg - dX * arma::diagmat( X.t() * g) - X * arma::diagmat(dX.t() * g) -
    // X * arma::diagmat(X.t() * dg);
    // dH = drg - X * arma::diagmat(X.t() * drg);

  }

  void retr(arguments_optim& x) {

    parameters = arma::vectorise(X * arma::diagmat(1 / sqrt(arma::sum(X % X, 0))));
    x.parameters(indices[0]) = parameters;

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

oblq* choose_oblq(Rcpp::List manifold_setup) {

  oblq* mymanifold = new oblq();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];
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
