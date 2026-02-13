/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Unit manifold (Sphere):

class unit:public manifolds {

public:

  arma::uvec indices;
  arma::vec X, dX, g, dg;

  void param(arguments_optim& x) {

    X = x.parameters(indices);

  }

  void proj(arguments_optim& x) {

    g = x.g.elem(indices);
    double v = arma::accu(X % g);
    x.rg.elem(indices) = g - X * v;

  }

  void hess(arguments_optim& x) {

    g = x.g.elem(indices);
    dg = x.dg.elem(indices);
    dX = x.dparameters.elem(indices);
    arma::mat drg = -dX * X.t() * g - X * dX.t() * g;
    // dH = drg - X * X.t() * drg;
    double v2 = arma::accu(X % g);
    arma::vec term = drg - v2 * dX;
    x.dH.elem(indices) = term - X * v2;

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) = X / sqrt(arma::accu(X % X));

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

unit* choose_unit(Rcpp::List manifold_setup) {

  unit* mymanifold = new unit();

  arma::uvec indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
