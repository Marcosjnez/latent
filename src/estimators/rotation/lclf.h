/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Linear CLF
 */

class lclf: public estimators {

public:

  double epsilon; // Provide this
  arma::mat lambda; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);
  arma::uvec lower, larger;
  double a, b, f1, f2;

  void param(arguments_optim& x) {

    X = arma::reshape(parameters, q, q);

    if(orth) {

      L = lambda*X;
      Phi = arma::mat(q, q, arma::fill::eye);

    } else {

      Phi = X.t() * X;
      Inv_X = arma::inv(X);
      L = lambda * Inv_X.t();

    }

    b = 1 / (2*epsilon);
    a = epsilon - b*epsilon*epsilon;

    arma::mat absL = arma::abs(L);
    lower = arma::find(absL <= epsilon);
    larger = arma::find(absL > epsilon);
    f1 = arma::accu(a + b*absL.elem(lower) % absL.elem(lower));
    f2 = arma::accu(absL.elem(larger));

  }

  void F(arguments_optim& x) {

    f = f1 + f2;

  }

  void G(arguments_optim& x) {

    gL.set_size(arma::size(L));
    gL.elem(lower) = 2*b*L.elem(lower);
    gL.elem(larger) = arma::sign(L.elem(larger));

    if(orth) {

      g = lambda.t() * gL;

    } else {

      g = - Inv_X.t() * gL.t() * L;

    }

  }

  void dG(arguments_optim& x) {

    dX = arma::reshape(dparameters, q, q);
    g = arma::reshape(g, q, q);

    dgL.set_size(arma::size(L));
    dgL.zeros();

    if(orth) {

      dL = lambda * dX;

      dgL.elem(lower) = 2*b*dL.elem(lower);

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();

      dgL.elem(lower) = 2*b*dL.elem(lower);

      dg = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;

    }

  }

  void outcomes(arguments_optim& x) {}

};

lclf* choose_lclf(const Rcpp::List& estimator_setup) {

  lclf* myestimator = new lclf();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  double epsilon = estimator_setup["epsilon"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int p = lambda.n_rows;
  int q = lambda.n_cols;

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->epsilon = epsilon;

  return myestimator;

}
