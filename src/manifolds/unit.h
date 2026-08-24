/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
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

    double xg = arma::dot(X, g);
    arma::vec h = dg-xg*dX;

    x.dH.elem(indices) =
      h-X*arma::dot(X, h);

    // dH = drg - X * X.t() * drg;

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) = X / sqrt(arma::accu(X % X));

  }

  void tangent_basis(arguments_optim& x) {

    X = x.parameters(indices);
    arma::mat T_local = tangent_complement(X);

    T.zeros(x.nparam, T_local.n_cols);

    for(arma::uword j=0L; j < T_local.n_cols; ++j) {
      for(arma::uword i=0L; i < indices.n_elem; ++i) {
        T(indices[i], j) = T_local(i, j);
      }
    }

  }

  void dconstraints(arguments_optim& x) {

    if(indices.n_elem == 0L) {
      return;
    }

    X = x.parameters(indices);

    arma::uword ndconstr = x.dconstr.n_cols;
    x.dconstr.resize(x.transparameters.n_elem, ndconstr+1L);

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    for(arma::uword i=0L; i < trans_indices.n_elem; ++i) {
      x.dconstr(trans_indices[i], ndconstr) = 2.00*X[i];
    }

  }

  void outcomes(arguments_optim& x) {

    tangent_basis(x);

    matrices.resize(1);
    matrices[0] = T;
    names_matrices.resize(1);
    names_matrices[0] = "tangent_basis";

  }

};

unit* choose_unit(Rcpp::List manifold_setup) {

  unit* mymanifold = new unit();

  arma::uvec indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
