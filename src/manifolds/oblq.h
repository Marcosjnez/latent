// Oblique manifold:

class oblq:public manifolds {

public:

  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param() {

    X = arma::reshape(parameters, q, q);

  }

  void proj() {

    g.reshape(q, q);
    rg = g - X * arma::diagmat(X.t() * g);

  }

  void hess() {

    g.reshape(q, q);
    dg.reshape(q, q);
    dX = arma::reshape(dparameters, q, q);
    dH = dg - dX * arma::diagmat(X.t() * g) -
      X * arma::diagmat(X.t() * dg);
    // arma::mat drg = dg - dX * arma::diagmat( X.t() * g) - X * arma::diagmat(dX.t() * g) -
    // X * arma::diagmat(X.t() * dg);
    // dH = drg - X * arma::diagmat(X.t() * drg);

  }

  void retr() {

    parameters = arma::vectorise(X * arma::diagmat(1 / sqrt(arma::sum(X % X, 0))));

  }

};
