// Orthogonal manifold:

class orth:public manifolds {

public:

  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param() {

    X = arma::reshape(parameters, q, q);

  }

  void proj() {

    g.reshape(q, q);
    rg = X * skew(X.t() * g);

  }

  void hess() {

    g.reshape(q, q);
    dg.reshape(q, q);
    dX = arma::reshape(dparameters, q, q);
    arma::mat drg = dg - dX * symm(X.t() * g);
    dH = X * skew(X.t() * drg);

  }

  void retr() {

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    parameters = arma::vectorise(Q);

  }

};
