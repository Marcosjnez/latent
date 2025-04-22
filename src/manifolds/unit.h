// Unit manifold (Sphere):

class unit:public manifolds {

public:

  arma::vec X;
  arma::vec dX;

  void param() {

    X = parameters;

  }

  void proj() {

    double v = arma::accu(X % g);
    rg = g - X * v;

  }

  void hess() {

    dX = dparameters;
    arma::mat drg = -dX * X.t() * g - X * dX.t() * g;
    // dH = drg - X * X.t() * drg;
    double v2 = arma::accu(X % g);
    arma::vec term = drg - v2 * dX;
    dH = term - X * v2;

  }

  void retr() {

    parameters = X / sqrt(arma::accu(X % X));

  }

};
