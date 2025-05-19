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

  arma::uvec indices = manifold_setup["indices"];
  mymanifold->indices = indices;

  return mymanifold;

}
