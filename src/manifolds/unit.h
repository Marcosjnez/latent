/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Unit manifold (Sphere):

class unit:public manifolds {

public:

  arma::vec X;
  arma::vec dX;
  arma::vec parameters, dir, g, dg, dparameters;

  void param(arguments_optim& x) {

    parameters = x.parameters(indices[0]);
    X = parameters;

  }

  void proj(arguments_optim& x) {

    g = x.g.elem(indices[0]);
    double v = arma::accu(X % g);
    x.rg.elem(indices[0]) = g - X * v;

  }

  void hess(arguments_optim& x) {

    dg = x.dg.elem(indices[0]);
    dparameters = x.dparameters.elem(indices[0]);
    dX = dparameters;
    arma::mat drg = -dX * X.t() * g - X * dX.t() * g;
    // dH = drg - X * X.t() * drg;
    double v2 = arma::accu(X % g);
    arma::vec term = drg - v2 * dX;
    x.dH.elem(indices[0]) = term - X * v2;

  }

  void retr(arguments_optim& x) {

    parameters = X / sqrt(arma::accu(X % X));
    x.parameters(indices[0]) = parameters;

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

unit* choose_unit(Rcpp::List manifold_setup) {

  unit* mymanifold = new unit();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
