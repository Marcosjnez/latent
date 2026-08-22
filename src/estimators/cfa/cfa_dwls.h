/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis (diagonally weighted least-squares)
 *
 * The estimator always receives model-implied and sample covariance matrices,
 * followed by model-implied and sample observed means.
 */

class cfa_dwls: public estimators {

public:

  int p;
  double loss, w;
  arma::uvec indices_S, indices_Shat, indices_meanshat, indices_means;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;
  arma::vec meanshat, means, delta, dmeanshat, dmeans, w_means;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);
    meanshat = x.transparameters(indices_meanshat);
    means = x.transparameters(indices_means);

    residuals = S-Shat;
    W_residuals = W % residuals;
    delta = means-meanshat;

  }

  void F(arguments_optim& x) {

    loss = w*0.5*(
      arma::accu(residuals % W_residuals)+
      arma::accu(w_means % delta % delta)
    );

    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::vec wdelta = w_means % delta;

    x.grad.elem(indices_S) += w*arma::vectorise(W_residuals);
    x.grad.elem(indices_Shat) += w*arma::vectorise(-W_residuals);
    x.grad.elem(indices_meanshat) += -w*wdelta;
    x.grad.elem(indices_means) += w*wdelta;

  }

  void dG(arguments_optim& x) {

    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dmeanshat = x.dtransparameters(indices_meanshat);
    dmeans = x.dtransparameters(indices_means);

    arma::vec ddelta = dmeans-dmeanshat;
    arma::vec wd_delta = w_means % ddelta;

    x.dgrad.elem(indices_S) +=
      w*arma::vectorise(W % dS-W % dShat);
    x.dgrad.elem(indices_Shat) +=
      w*arma::vectorise(W % dShat-W % dS);
    x.dgrad.elem(indices_meanshat) += -w*wd_delta;
    x.dgrad.elem(indices_means) += w*wd_delta;

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();
    arma::mat W_residuals_indep = W % residuals_indep;

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

    matrices.resize(2);
    matrices[0] = residuals;
    matrices[1] = delta;

  }

};

cfa_dwls* choose_cfa_dwls(const Rcpp::List& estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  double w = estimator_setup["w"];
  arma::mat W = estimator_setup["W"];
  arma::vec w_means = estimator_setup["w_means"];

  if(p < 1 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_dwls received invalid p or w values.");
  }

  if(indices.size() != 4L) {
    Rcpp::stop("cfa_dwls requires Shat, S, meanshat, and means indices.");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p) ||
     indices[2].n_elem != static_cast<arma::uword>(p) ||
     indices[3].n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The cfa_dwls parameter indices have incompatible dimensions.");
  }

  if(W.n_rows != static_cast<arma::uword>(p) ||
     W.n_cols != static_cast<arma::uword>(p)) {
    Rcpp::stop("cfa_dwls requires W to be a p by p matrix.");
  }

  if(w_means.n_elem != static_cast<arma::uword>(p) ||
     !w_means.is_finite() || arma::any(w_means < 0.0)) {
    Rcpp::stop("cfa_dwls requires one finite non-negative mean weight per observed variable.");
  }

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_meanshat = indices[2];
  myestimator->indices_means = indices[3];
  myestimator->p = p;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->w_means = w_means;
  myestimator->Shat.zeros(p, p);
  myestimator->S.zeros(p, p);
  myestimator->dShat.zeros(p, p);
  myestimator->dS.zeros(p, p);
  myestimator->meanshat.zeros(p);
  myestimator->means.zeros(p);
  myestimator->delta.zeros(p);
  myestimator->dmeanshat.zeros(p);
  myestimator->dmeans.zeros(p);

  return myestimator;

}
