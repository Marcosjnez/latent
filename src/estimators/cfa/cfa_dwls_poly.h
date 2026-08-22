/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis for ordinal variables
 * (diagonally weighted least-squares)
 *
 * The discrepancy function compares:
 *   - sample and model-implied polychoric correlations;
 *   - sample and model-implied standardized thresholds.
 */

class cfa_dwls_poly: public estimators {

public:

  int p, ntau;
  double loss, w;
  arma::uvec indices_S, indices_Shat, indices_taus, indices_tauhat;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;
  arma::vec taus, tauhat, dtaus, dtauhat, threshold_residuals,
    w_thresholds, W_threshold_residuals;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);
    tauhat = x.transparameters(indices_tauhat);
    taus = x.transparameters(indices_taus);

    residuals = S-Shat;
    W_residuals = W % residuals;
    threshold_residuals = taus-tauhat;
    W_threshold_residuals = w_thresholds % threshold_residuals;

  }

  void F(arguments_optim& x) {

    loss = w*0.5*(
      arma::accu(residuals % W_residuals)+
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

    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dtauhat = x.dtransparameters(indices_tauhat);
    dtaus = x.dtransparameters(indices_taus);

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

    (void)x;

    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();
    arma::mat W_residuals_indep = W % residuals_indep;

    // Thresholds are saturated in the independence model.
    double loss_indep =
      0.5*arma::accu(residuals_indep % W_residuals_indep);
    double loss_sat = 0.00;

    doubles.resize(7);
    doubles[0] = loss;
    doubles[1] = loss_indep;
    doubles[2] = loss_sat;
    doubles[3] = 0.00;
    doubles[4] = 0.00;
    doubles[5] = 0.00;
    doubles[6] = 0.00;

    vectors.resize(1);
    vectors[0] = threshold_residuals;
    names_vectors.resize(1);
    names_vectors[0] = "threshold residuals";

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

  if(p < 2 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_dwls_poly received invalid p or w values.");
  }

  if(indices.size() != 4L) {
    Rcpp::stop("cfa_dwls_poly requires Shat, S, tauhat, and taus indices.");
  }

  const arma::uword ntau = indices[2].n_elem;

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p) ||
     indices[3].n_elem != ntau || ntau < 1L) {
    Rcpp::stop("The cfa_dwls_poly parameter indices have incompatible dimensions.");
  }

  if(W.n_rows != static_cast<arma::uword>(p) ||
     W.n_cols != static_cast<arma::uword>(p)) {
    Rcpp::stop("cfa_dwls_poly requires W to be a p by p matrix.");
  }

  if(w_thresholds.n_elem != ntau ||
     !w_thresholds.is_finite() || arma::any(w_thresholds < 0.0)) {
    Rcpp::stop("cfa_dwls_poly requires one finite non-negative weight per threshold.");
  }

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_tauhat = indices[2];
  myestimator->indices_taus = indices[3];
  myestimator->p = p;
  myestimator->ntau = static_cast<int>(ntau);
  myestimator->W = W;
  myestimator->w = w;
  myestimator->w_thresholds = w_thresholds;
  myestimator->Shat.zeros(p, p);
  myestimator->S.zeros(p, p);
  myestimator->dShat.zeros(p, p);
  myestimator->dS.zeros(p, p);
  myestimator->tauhat.zeros(ntau);
  myestimator->taus.zeros(ntau);
  myestimator->dtauhat.zeros(ntau);
  myestimator->dtaus.zeros(ntau);
  myestimator->threshold_residuals.zeros(ntau);
  myestimator->W_threshold_residuals.zeros(ntau);

  return myestimator;

}
