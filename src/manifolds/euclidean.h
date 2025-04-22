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

};
