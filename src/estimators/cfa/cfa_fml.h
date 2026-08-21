/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 21/08/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood discrepancy)
 *
 * When the sample covariance is positive definite, the ordinary FML discrepancy
 * is used. For rank-deficient missing-pattern covariance matrices, terms that
 * depend only on the fixed sample covariance are omitted. This leaves parameter
 * estimates, model gradients, and log-likelihood contributions unchanged while
 * allowing singleton and low-frequency missingness patterns.
 */

class cfa_fml: public estimators {

public:

  int p, n;
  double loss, w, logdetS, logdetShat, plogpi2;
  bool S_sympd;
  arma::uvec indices_S, indices_Shat;
  arma::mat S, Shat, dShat, dS, Shat_inv, S_inv, gShat, gS, I;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);
    S = 0.5*(S + S.t());

    if(!Shat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Shat);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      Shat = eigvec * arma::diagmat(d) * eigvec.t();
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

    if(S_sympd) {
      loss = w*(logdetShat - logdetS + arma::accu(S % Shat_inv) - p);
    } else {
      loss = w*(logdetShat + arma::accu(S % Shat_inv));
    }

    x.f += loss;

  }

  void G(arguments_optim& x) {

    gShat = Shat_inv * (I - S * Shat_inv);

    if(S_sympd) {
      gS = Shat_inv - S_inv;
    } else {
      gS = Shat_inv;
    }

    x.grad.elem(indices_S) += w*arma::vectorise(gS);
    x.grad.elem(indices_Shat) += w*arma::vectorise(gShat);

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgShat = dShat_inv * (I - S * Shat_inv) -
      Shat_inv * dS * Shat_inv -
      Shat_inv * S * dShat_inv;

    arma::mat dgS = dShat_inv;

    if(S_sympd) {
      arma::mat dS_inv = -S_inv * dS * S_inv;
      dgS -= dS_inv;
    }

    x.dgrad.elem(indices_S) += w*arma::vectorise(dgS);
    x.dgrad.elem(indices_Shat) += w*arma::vectorise(dgShat);

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    double loss_indep = S_sympd ?
      w*(-logdetS + arma::trace(S) - p) :
      w*arma::trace(S);
    double loglik = n*0.5*(-plogpi2 -
      arma::log_det_sympd(Shat) - arma::accu(S % Shat_inv));
    double loglik_indep = n*0.5*(-plogpi2 - arma::trace(S));

    doubles.resize(7);
    doubles[0] = loss;
    doubles[1] = loss_indep;
    doubles[2] = S_sympd ? 0.0 : NA_REAL;
    doubles[3] = loglik;
    doubles[4] = loglik_indep;
    doubles[5] = S_sympd ? n*0.5*(-plogpi2 - logdetS - p) : NA_REAL;
    doubles[6] = 0.0;

    matrices.resize(1);
    matrices[0] = S - Shat;

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

  if(indices.size() != 2L) {
    Rcpp::stop("cfa_fml requires model and sample covariance indices");
  }

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->p = p;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->Shat.zeros(p, p);
  myestimator->S.zeros(p, p);
  myestimator->dShat.zeros(p, p);
  myestimator->dS.zeros(p, p);
  myestimator->Shat_inv.zeros(p, p);
  myestimator->S_inv.zeros(p, p);
  myestimator->gShat.zeros(p, p);
  myestimator->gS.zeros(p, p);
  myestimator->I.eye(p, p);
  myestimator->plogpi2 = p*std::log(arma::datum::pi*2.0);
  myestimator->S_sympd = false;

  return myestimator;

}
