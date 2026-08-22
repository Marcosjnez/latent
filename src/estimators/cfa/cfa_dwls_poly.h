/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis for ordinal variables
 * (diagonally weighted least-squares)
 *
 * The discrepancy function includes both:
 *   - polychoric correlation residuals
 *   - threshold residuals
 */

class cfa_dwls_poly: public estimators {

public:

  int p;
  double loss, w;
  arma::uvec indices_S, indices_Shat, indices_taus, indices_tauhat,
    diag, lower_diag;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;
  arma::vec taus, tauhat, dtaus, dtauhat, threshold_residuals,
    w_thresholds, W_threshold_residuals;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    taus = x.transparameters(indices_taus);
    tauhat = x.transparameters(indices_tauhat);

    residuals = S-Shat;
    W_residuals = W % residuals;

    threshold_residuals = taus-tauhat;
    W_threshold_residuals = w_thresholds % threshold_residuals;

  }

  void F(arguments_optim& x) {

    loss = w*0.5*(
      arma::accu(residuals % W_residuals) +
      arma::accu(threshold_residuals % W_threshold_residuals)
    );

    x.f += loss;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices_S) += w*arma::vectorise(W_residuals);
    x.grad.elem(indices_Shat) += w*arma::vectorise(-W_residuals);

    x.grad.elem(indices_taus) += w*W_threshold_residuals;
    x.grad.elem(indices_tauhat) += -w*W_threshold_residuals;

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dtaus = x.dtransparameters(indices_taus);
    dtauhat = x.dtransparameters(indices_tauhat);

    arma::vec dthreshold_residuals = dtaus-dtauhat;
    arma::vec W_dthreshold_residuals =
      w_thresholds % dthreshold_residuals;

    x.dgrad.elem(indices_S) +=
      w*arma::vectorise(W % dS-W % dShat);
    x.dgrad.elem(indices_Shat) +=
      w*arma::vectorise(W % dShat-W % dS);

    x.dgrad.elem(indices_taus) += w*W_dthreshold_residuals;
    x.dgrad.elem(indices_tauhat) += -w*W_dthreshold_residuals;

  }

  void outcomes(arguments_optim& x) {

    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();

    arma::mat W_residuals_indep = W % residuals_indep;

    // Thresholds are saturated in the independence model, so their
    // contribution to the independence discrepancy is zero.
    double loss_indep =
      0.5*arma::accu(residuals_indep % W_residuals_indep);

    double loss_sat = 0.00;

    doubles.resize(7);
    doubles[0] =  loss;        // loss actual model
    doubles[1] =  loss_indep;  // loss independence model
    doubles[2] =  loss_sat;    // loss saturated model
    doubles[3] =  0.00;        // loglik actual model
    doubles[4] =  0.00;        // loglik independence model
    doubles[5] =  0.00;        // loglik saturated model
    doubles[6] =  0.00;        // penalty

    vectors.resize(1);
    vectors[0] = threshold_residuals;

    matrices.resize(1);
    matrices[0] = residuals;

  }

};

cfa_dwls_poly* choose_cfa_dwls_poly(const Rcpp::List& estimator_setup) {

  cfa_dwls_poly* myestimator = new cfa_dwls_poly();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  double w = estimator_setup["w"];
  arma::mat W = estimator_setup["W"];
  arma::vec w_thresholds = estimator_setup["w_thresholds"];

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p+1, p*p-1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_tauhat = indices[2];
  myestimator->indices_taus = indices[3];

  myestimator->p = p;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->w_thresholds = w_thresholds;

  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
