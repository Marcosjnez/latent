/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * Extended Target
 */

class xtarget: public estimators {

public:

  double w;
  arma::mat lambda, Target, Weight, PhiTarget, PhiWeight,
  Weight2, PhiWeight2; // Fixed: Provide these in choose_estimator
  bool orth; // TRUE: orthogonal; FALSE: oblique and poblique
  arma::mat dL, dP, gL, dgL, gP, dgP, Inv_X, L, L2, f1, f2, hL;
  arma::mat Phi;
  arma::mat X = arma::mat(q, q);
  arma::mat dX = arma::mat(q, q);

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

    f1 = Weight % (L - Target);
    f2 = PhiWeight % (Phi - PhiTarget);

  }

  void F(arguments_optim& x) {

    f = 0.5*arma::accu(f1 % f1) + 0.25*w*arma::accu(f2 % f2);

  }

  void G(arguments_optim& x) {

    gL = Weight % f1;

    if(orth) {

      g = lambda.t() * gL;

    } else {

      gP = w * PhiWeight % f2;
      arma::mat g1 = - Inv_X.t() * gL.t() * L;
      arma::mat g2 = X * gP;
      g = g1 + g2;

    }

  }

  void dG(arguments_optim& x) {

    dX = arma::reshape(dparameters, q, q);
    g = arma::reshape(g, q, q);

    if(orth) {

      dL = lambda * dX;

      dgL = Weight2 % dL;

      dg = lambda.t() * dgL;

    } else {

      arma::mat Inv_X_dt = Inv_X * dX;
      dL = - L * Inv_X_dt.t();
      dP = X.t() * dX; dP += dP.t();

      dgL = Weight2 % dL;
      dgP = w * PhiWeight2 % dP;

      arma::mat dg1 = - g * Inv_X_dt.t() -
        (dX * Inv_X).t() * g - (dgL * Inv_X).t() * L;
      arma::mat dg2 = dX * gP + X * dgP;

      dg = dg1 + dg2;

    }

  }

  void H(arguments_optim& x) {

    // Rcpp::stop("H not available");
    hess.set_size(parameters.n_elem, parameters.n_elem); hess.zeros();

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void outcomes(arguments_optim& x) {

    /*
     * Compute the modified hessian (modhessian)
     */

    int nhessian = lambda.n_elem;
    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    hL = arma::diagmat(arma::vectorise(Weight2));

    arma::mat wPW2 = arma::diagmat(arma::vectorise(w*PhiWeight2));
    arma::mat hP = wPW2 + wPW2 * dxt(q, q);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(Phi), -1);
    hP = hP;

    modhessian = hL; // FIX to export hP

  }

};

xtarget* choose_xtarget(const Rcpp::List& estimator_setup) {

  xtarget* myestimator = new xtarget();

  arma::mat lambda = estimator_setup["lambda"];
  bool orth = estimator_setup["orth"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat Target = estimator_setup["Target"];
  arma::mat Weight = estimator_setup["Weight"];
  arma::mat PhiTarget = estimator_setup["PhiTarget"];
  arma::mat PhiWeight = estimator_setup["PhiWeight"];
  double w = estimator_setup["w"];

  arma::mat Weight2 = Weight % Weight;
  arma::mat PhiWeight2 = PhiWeight % PhiWeight;
  int p = lambda.n_rows;
  int q = lambda.n_cols;

  myestimator->lambda = lambda;
  myestimator->orth = orth;
  myestimator->indices = indices;
  myestimator->Target = Target;
  myestimator->Weight = Weight;
  myestimator->Weight2 = Weight2;
  myestimator->PhiTarget = PhiTarget;
  myestimator->PhiWeight = PhiWeight;
  myestimator->PhiWeight2 = PhiWeight2;
  myestimator->w = w;
  myestimator->p = p;
  myestimator->q = q;

  return myestimator;

}
