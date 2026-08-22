/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood discrepancy)
 *
 * General CFA FML estimator. A mean structure is included whenever model and
 * sample mean indices are supplied. If they are omitted, the mean contribution
 * is exactly zero and the estimator reduces to the covariance-only FML model.
 *
 * Rank-deficient sample covariance matrices are permitted for direct FIML
 * missingness-pattern contributions. Terms depending only on a singular fixed
 * sample covariance are omitted while model-dependent terms remain unchanged.
 */

class cfa_fml: public estimators {

public:

  int p, n;
  double loss, w, logdetS, logdetShat, plogpi2;
  bool S_sympd, use_means;
  arma::uvec indices_S, indices_Shat, indices_nu, indices_means;
  arma::mat S, Shat, dShat, dS, Shat_inv, S_inv, gShat, gS, I;
  arma::vec nu, means, delta, dnu, dmeans, gnu, gmeans;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);
    S = 0.5*(S+S.t());

    if(use_means) {
      nu = x.transparameters(indices_nu);
      means = x.transparameters(indices_means);
      delta = means-nu;
    } else {
      delta.zeros();
    }

    if(!Shat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Shat);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      Shat = eigvec*arma::diagmat(d)*eigvec.t();
    }

    Shat_inv = arma::inv_sympd(Shat);
    S_sympd = S.is_sympd();

    if(S_sympd) {
      S_inv = arma::inv_sympd(S);
      logdetS = arma::log_det_sympd(S);
    } else {
      S_inv.zeros(p, p);
      logdetS = 0.0;
    }

  }

  void F(arguments_optim& x) {

    logdetShat = arma::log_det_sympd(Shat);
    double mean_term = arma::as_scalar(delta.t()*Shat_inv*delta);

    if(S_sympd) {
      loss = w*(logdetShat-logdetS+arma::accu(S % Shat_inv)-p+
        mean_term);
    } else {
      loss = w*(logdetShat+arma::accu(S % Shat_inv)+mean_term);
    }

    x.f += loss;

  }

  void G(arguments_optim& x) {

    gShat = Shat_inv-
      Shat_inv*S*Shat_inv-
      Shat_inv*delta*delta.t()*Shat_inv;

    if(S_sympd) {
      gS = Shat_inv-S_inv;
    } else {
      gS = Shat_inv;
    }

    x.grad.elem(indices_S) += w*arma::vectorise(gS);
    x.grad.elem(indices_Shat) += w*arma::vectorise(gShat);

    if(use_means) {
      gnu = -2.0*Shat_inv*delta;
      gmeans = 2.0*Shat_inv*delta;
      x.grad.elem(indices_nu) += w*gnu;
      x.grad.elem(indices_means) += w*gmeans;
    }

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    arma::vec ddelta(p, arma::fill::zeros);

    if(use_means) {
      dnu = x.dtransparameters(indices_nu);
      dmeans = x.dtransparameters(indices_means);
      ddelta = dmeans-dnu;
    }

    arma::mat dD = ddelta*delta.t()+delta*ddelta.t();
    arma::mat dShat_inv = -Shat_inv*dShat*Shat_inv;

    arma::mat dgShat =
      dShat_inv
      -dShat_inv*S*Shat_inv
      -Shat_inv*dS*Shat_inv
      -Shat_inv*S*dShat_inv
      -dShat_inv*delta*delta.t()*Shat_inv
      -Shat_inv*dD*Shat_inv
      -Shat_inv*delta*delta.t()*dShat_inv;

    arma::mat dgS = dShat_inv;

    if(S_sympd) {
      arma::mat dS_inv = -S_inv*dS*S_inv;
      dgS -= dS_inv;
    }

    x.dgrad.elem(indices_S) += w*arma::vectorise(dgS);
    x.dgrad.elem(indices_Shat) += w*arma::vectorise(dgShat);

    if(use_means) {
      arma::vec dgnu = -2.0*(dShat_inv*delta+Shat_inv*ddelta);
      arma::vec dgmeans = 2.0*(dShat_inv*delta+Shat_inv*ddelta);
      x.dgrad.elem(indices_nu) += w*dgnu;
      x.dgrad.elem(indices_means) += w*dgmeans;
    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    double mean_term = arma::as_scalar(delta.t()*Shat_inv*delta);
    double mean_term_indep = arma::as_scalar(delta.t()*delta);
    double loss_indep = S_sympd ?
      w*(-logdetS+arma::trace(S)-p+mean_term_indep) :
      w*(arma::trace(S)+mean_term_indep);
    double loglik = n*0.5*(-plogpi2-
      arma::log_det_sympd(Shat)-arma::accu(S % Shat_inv)-mean_term);
    double loglik_indep = n*0.5*(-plogpi2-arma::trace(S)-
      mean_term_indep);

    doubles.resize(7);
    doubles[0] = loss;
    doubles[1] = loss_indep;
    doubles[2] = S_sympd ? w*mean_term : NA_REAL;
    doubles[3] = loglik;
    doubles[4] = loglik_indep;
    doubles[5] = S_sympd ? n*0.5*(-plogpi2-logdetS-p) : NA_REAL;
    doubles[6] = 0.0;

    matrices.resize(use_means ? 2L : 1L);
    matrices[0] = S-Shat;

    if(use_means) {
      matrices[1] = delta;
    }

  }

};

cfa_fml* choose_cfa_fml(const Rcpp::List& estimator_setup) {

  cfa_fml* myestimator = new cfa_fml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double w = estimator_setup["w"];
  int n = estimator_setup["n"];
  int p = estimator_setup["p"];

  if(p < 1 || n < 1 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_fml received invalid p, n, or w values");
  }

  if(indices.size() != 2L && indices.size() != 4L) {
    Rcpp::stop("cfa_fml requires two or four parameter-index objects");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("The cfa_fml covariance indices have incompatible dimensions");
  }

  bool use_means = false;
  arma::uvec indices_nu;
  arma::uvec indices_means;

  if(indices.size() == 4L) {

    const bool has_nu = indices[2].n_elem > 0L;
    const bool has_means = indices[3].n_elem > 0L;

    if(has_nu != has_means) {
      Rcpp::stop("cfa_fml requires both model and sample means, or neither");
    }

    if(has_nu) {
      if(indices[2].n_elem != static_cast<arma::uword>(p) ||
         indices[3].n_elem != static_cast<arma::uword>(p)) {
        Rcpp::stop("The cfa_fml mean indices have incompatible dimensions");
      }
      use_means = true;
      indices_nu = indices[2];
      indices_means = indices[3];
    }

  }

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_nu = indices_nu;
  myestimator->indices_means = indices_means;
  myestimator->p = p;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->use_means = use_means;
  myestimator->Shat.zeros(p, p);
  myestimator->S.zeros(p, p);
  myestimator->dShat.zeros(p, p);
  myestimator->dS.zeros(p, p);
  myestimator->Shat_inv.zeros(p, p);
  myestimator->S_inv.zeros(p, p);
  myestimator->gShat.zeros(p, p);
  myestimator->gS.zeros(p, p);
  myestimator->I.eye(p, p);
  myestimator->nu.zeros(p);
  myestimator->means.zeros(p);
  myestimator->delta.zeros(p);
  myestimator->dnu.zeros(p);
  myestimator->dmeans.zeros(p);
  myestimator->gnu.zeros(p);
  myestimator->gmeans.zeros(p);
  myestimator->plogpi2 = p*std::log(arma::datum::pi*2.0);
  myestimator->S_sympd = false;

  return myestimator;

}
