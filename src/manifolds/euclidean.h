/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Euclidean manifold:

class euclidean:public manifolds {

public:

  arma::vec X;
  arma::vec dX;

  void param() {

    // X = parameters;

  }

  void proj() {

    rg = g;

  }

  void hess() {

    dH = dg;

  }

  void retr() {

    parameters = parameters;

  }

  void dconstraints() {

  }

  void outcomes() {

  }

};

euclidean* choose_euclidean(Rcpp::List manifold_setup) {

  euclidean* mymanifold = new euclidean();

  std::vector<arma::uvec> indices = manifold_setup["indices"];
  mymanifold->indices = indices;

  return mymanifold;

}
