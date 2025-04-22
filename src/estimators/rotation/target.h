/*
 * Target
 */

class target: public estimators {

public:

  arma::mat lambda, Target, Weight, Weight2; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, f1, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

  void param() {

    X = arma::reshape(parameters, q, q);

    if(orth) {

      L = lambda*X;
      Phi = arma::mat(q, q, arma::fill::eye);

    } else {

      Phi = X.t() * X;
      Inv_X = arma::inv(X);
      L = lambda * Inv_X.t();

    }

    f1 = Weight % (L - Target);

  }

  void F() {

    f = 0.5*arma::accu(f1 % f1);

  }

  void G() {

    gL = Weight % f1;

    if(orth) {

      g = lambda.t() * gL;

    } else {

      g = - Inv_X.t() * gL.t() * L;

    }

  }

  void dG() {

    dX = arma::reshape(dparameters, q, q);
    g = arma::reshape(g, q, q);

    if(orth) {

      dL = lambda * dX;

      dgL = Weight2 % dL;

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      dgL = Weight2 % dL;

      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void H() {

    // Rcpp::stop("H not available");
    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() {}

  void M() {}

  void outcomes() {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    arma::mat W2 = Weight % Weight;
    modhessian = arma::diagmat(arma::vectorise(W2));

  }

};
