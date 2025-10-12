/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/10/2025
 */

// Multinomial manifold (simplex):

class simplex:public manifolds {

public:

  arma::vec X;
  arma::vec dX;
  arma::vec parameters, dir, g, dg;

  void param(arguments_optim& x) {

    parameters = x.parameters(indices[0]);
    X = parameters;
    dir = x.dir(indices[0]);

  }

  void proj(arguments_optim& x) {

    g = x.g.elem(indices[0]);
    arma::vec xegrad = X % g;
    double v = arma::accu(xegrad);
    x.rg.elem(indices[0]) = xegrad - X * v;

  }

  void hess(arguments_optim& x) {

    dg = x.dg.elem(indices[0]);
    x.dH.elem(indices[0]) = x.dg.elem(indices[0]);
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
    for(int j=0; j < X.size(); ++j) {
      Y[j] = X[j] * exp(x.ss * (dir[j] / X[j]));
    }
    Y = Y / arma::accu(Y);

    double eps = arma::datum::eps;
    // Rf_error("OK 61");
    Y.elem(arma::find(Y < eps)).fill(eps);
    parameters = Y;
    x.parameters(indices[0]) = parameters;

  }

  void dconstraints(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

  }

};

simplex* choose_simplex(const Rcpp::List& manifold_setup) {

  simplex* mymanifold = new simplex();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
