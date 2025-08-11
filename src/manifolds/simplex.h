/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Multinomial manifold (simplex):

class simplex:public manifolds {

public:

  arma::vec X;
  arma::vec dX;

  void param() {

    X = parameters;

  }

  void proj() {

    arma::vec xegrad = X % g;
    double v = arma::accu(xegrad);
    rg = xegrad - X * v;

  }

  void hess() {

    // dX = dparameters;
    // arma::mat drg = -dX * X.t() * g - X * dX.t() * g;
    // // dH = drg - X * X.t() * drg;
    // double v2 = arma::accu(X % g);
    // arma::vec term = drg - v2 * dX;
    // dH = term - X * v2;

  }

  void retr() {

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

    X = X + ss*dir;
    // Rf_error("OK 53");
    arma::vec Y(X.size());
    for(int j=0; j < X.size(); ++j) {
      Y[j] = X[j] * exp(ss * (dir[j] / X[j]));
    }
    Y = Y / arma::accu(Y);

    double eps = arma::datum::eps;
    // Rf_error("OK 61");
    Y.elem(arma::find(Y < eps)).fill(eps);
    parameters = Y;

  }

  void dconstraints() {

  }

  void outcomes() {

  }

};

simplex* choose_simplex(const Rcpp::List& manifold_setup) {

  simplex* mymanifold = new simplex();

  // Provide these:
  std::vector<arma::uvec> indices = manifold_setup["indices"];

  mymanifold->indices = indices;

  return mymanifold;

}
