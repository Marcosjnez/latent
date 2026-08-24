/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// Multinomial manifold (simplex):

class simplex:public manifolds {

public:

  arma::uvec indices;
  arma::vec X, dX, dir, g, dg;

  void param(arguments_optim& x) {

    X = x.parameters(indices);
    dir = x.dir(indices);

  }

  void proj(arguments_optim& x) {

    g = x.g.elem(indices);
    arma::vec xegrad = X % g;
    double v = arma::accu(xegrad);
    x.rg.elem(indices) = xegrad - X * v;

  }

  void hess(arguments_optim& x) {

    g = x.g.elem(indices);
    dg = x.dg.elem(indices);
    dX = x.dparameters.elem(indices);

    double xg = arma::dot(X, g);
    double xdg = arma::dot(X, dg);
    double dxg = arma::dot(dX, g);

    x.dH.elem(indices) =
      X % dg+0.5*(dX % g)-
      X*(xdg+0.5*dxg)-0.5*xg*dX;

    // Previous draft:
    // dX = dparameters;
    // arma::mat drg = -dX * X.t() * g - X * dX.t() * g;
    // // dH = drg - X * X.t() * drg;
    // double v2 = arma::accu(X % g);
    // arma::vec term = drg - v2 * dX;
    // dH = term - X * v2;

  }

  void retr(arguments_optim& x) {

    // Rf_error("OK 37");
    double alpha = arma::accu(dir);

    // for (size_t i = 0; i < dir.n_elem; ++i) {
    //   Rprintf("parameter ");
    //   Rprintf("%u", i);
    //   Rprintf(" = ");
    //   Rprintf("%f", dir[i]);
    //   Rprintf("\n");
    // }
    // Rf_error("OK 40");

    dir = dir - alpha * X;
    // Rf_error("OK 50");

    X = X + x.ss*dir;
    // Rf_error("OK 53");
    arma::vec Y(X.size());
    for(arma::uword j=0L; j < X.size(); ++j) {
      Y[j] = X[j] * exp(x.ss * (dir[j] / X[j]));
    }
    Y = Y / arma::accu(Y);

    double eps = arma::datum::eps;
    // Rf_error("OK 61");
    Y.elem(arma::find(Y < eps)).fill(eps);
    x.parameters(indices) = Y;

  }

  void tangent_basis(arguments_optim& x) {

    const arma::uword n = indices.n_elem;
    const arma::uword dimension = n > 0L ? n-1L : 0L;

    T.zeros(x.nparam, dimension);

    for(arma::uword j=0L; j < dimension; ++j) {

      double scale = std::sqrt(
        static_cast<double>((j+1L)*(j+2L))
      );

      for(arma::uword i=0L; i <= j; ++i) {
        T(indices[i], j) = 1.00/scale;
      }

      T(indices[j+1L], j) =
        -static_cast<double>(j+1L)/scale;

    }

  }

  void dconstraints(arguments_optim& x) {

    if(indices.n_elem == 0L) {
      return;
    }

    // Expand the matrix of constraints derivatives to put in a new column
    // the constraint derivatives of this manifold:
    arma::uword ndconstr = x.dconstr.n_cols;
    x.dconstr.resize(x.transparameters.n_elem, ndconstr+1L);

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    for(arma::uword i=0L; i < trans_indices.n_elem; ++i) {
      x.dconstr(trans_indices[i], ndconstr) = 1.00;
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

simplex* choose_simplex(const Rcpp::List& manifold_setup) {

  simplex* mymanifold = new simplex();

  arma::uvec indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
