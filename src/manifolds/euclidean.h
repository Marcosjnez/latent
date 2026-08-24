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

    T.zeros(x.nparam, indices.n_elem);

    for(arma::uword i=0L; i < indices.n_elem; ++i) {
      T(indices[i], i) = 1.00;
    }

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

    tangent_basis(x);

    matrices.resize(1);
    matrices[0] = T;
    names_matrices.resize(1);
    names_matrices[0] = "tangent_basis";

  }

};

euclidean* choose_euclidean(Rcpp::List manifold_setup) {

  euclidean* mymanifold = new euclidean();

  arma::uvec indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
