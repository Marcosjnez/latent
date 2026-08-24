/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// Euclidean manifold:

class euclidean:public manifolds {

public:

  arma::uvec indices;

  void param(arguments_optim& x) {

  }

  void proj(arguments_optim& x) {

    x.rg.elem(indices) = x.g.elem(indices);

  }

  void hess(arguments_optim& x) {

    x.dH.elem(indices) = x.dg.elem(indices);

  }

  void retr(arguments_optim& x) {

  }

  void tangent_basis(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

euclidean* choose_euclidean(Rcpp::List manifold_setup) {

  euclidean* mymanifold = new euclidean();

  arma::uvec indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
