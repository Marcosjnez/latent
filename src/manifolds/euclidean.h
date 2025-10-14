/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/10/2025
 */

// Euclidean manifold:

class euclidean:public manifolds {

public:

  void param(arguments_optim& x) {

  }

  void proj(arguments_optim& x) {

    x.rg.elem(indices[0]) = x.g.elem(indices[0]);

  }

  void hess(arguments_optim& x) {

    x.dH.elem(indices[0]) = x.dg.elem(indices[0]);

  }

  void retr(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

euclidean* choose_euclidean(Rcpp::List manifold_setup) {

  euclidean* mymanifold = new euclidean();

  std::vector<arma::uvec> indices = manifold_setup["indices"];
  mymanifold->indices = indices;

  return mymanifold;

}
