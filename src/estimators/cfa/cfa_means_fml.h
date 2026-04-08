/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 08/04/2026
 */

/*
 * Confirmatory factor analysis with mean structure (maximum-likelihood)
 */

class cfa_means_fml: public estimators {

public:

  int p, n;
  double f, w, logdetS, logdetShat, plogpi2;
  arma::uvec indices_S, indices_Shat, indices_nu, indices_means, lower_diag;
  arma::mat S, Shat, residuals, dShat, dS, Shat_inv, S_inv, gShat, gS, I;
  arma::vec nu, means, delta, dnu, dmeans, gnu, gmeans;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    nu = x.transparameters(indices_nu);
    means = x.transparameters(indices_means);
    delta = means - nu;

    if(!Shat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Shat);
      arma::vec d = arma::clamp(eigval, 0.1, eigval.max());
      Shat = eigvec * arma::diagmat(d) * eigvec.t();
    }

    Shat_inv = arma::inv_sympd(Shat);

  }

  void F(arguments_optim& x) {

    logdetS = arma::log_det_sympd(S);
    logdetShat = arma::log_det_sympd(Shat);

    f = w*(logdetShat - logdetS + arma::accu(S % Shat_inv) - p +
      arma::as_scalar(delta.t() * Shat_inv * delta));
    x.f += f;

  }

  void G(arguments_optim& x) {

    S_inv = arma::inv_sympd(S);
    gS = Shat_inv - S_inv;
    gShat = Shat_inv - Shat_inv * S * Shat_inv - Shat_inv * delta * delta.t() * Shat_inv;
    gnu = -2.0 * Shat_inv * delta;
    gmeans = 2.0 * Shat_inv * delta;

    x.grad.elem(indices_S) += w * arma::vectorise(gS);
    x.grad.elem(indices_Shat) += w * arma::vectorise(gShat);
    x.grad.elem(indices_nu) += w * gnu;
    x.grad.elem(indices_means) += w * gmeans;

  }

  void dG(arguments_optim& x) {

    S_inv = arma::inv_sympd(S);

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dnu = x.dtransparameters(indices_nu);
    dmeans = x.dtransparameters(indices_means);

    arma::vec ddelta = dmeans - dnu;
    arma::mat dD = ddelta * delta.t() + delta * ddelta.t();

    arma::mat dS_inv = -S_inv * dS * S_inv;
    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;

    arma::mat dgS = dShat_inv - dS_inv;
    arma::mat dgShat =
      dShat_inv
      - dShat_inv * S * Shat_inv
    - Shat_inv * dS * Shat_inv
    - Shat_inv * S * dShat_inv
    - dShat_inv * delta * delta.t() * Shat_inv
    - Shat_inv * dD * Shat_inv
    - Shat_inv * delta * delta.t() * dShat_inv;

    arma::vec dgnu = -2.0 * (dShat_inv * delta + Shat_inv * ddelta);
    arma::vec dgmeans = 2.0 * (dShat_inv * delta + Shat_inv * ddelta);

    x.dgrad.elem(indices_S) += w * arma::vectorise(dgS);
    x.dgrad.elem(indices_Shat) += w * arma::vectorise(dgShat);
    x.dgrad.elem(indices_nu) += w * dgnu;
    x.dgrad.elem(indices_means) += w * dgmeans;

  }

  void outcomes(arguments_optim& x) {

    double mean_term = arma::as_scalar(delta.t() * Shat_inv * delta);
    double loglik = n * 0.5 * (-plogpi2 -
                               arma::log_det_sympd(Shat) -
                               arma::accu(S % Shat_inv) -
                               mean_term);
    double loglik_indep = n * 0.5 * (-plogpi2 -
                                     arma::trace(S));
    arma::mat Sinv = arma::inv_sympd(S);
    double loglik_sat = n * 0.5 * (-plogpi2 -
                                   arma::log_det_sympd(S) -
                                   arma::accu(S % Sinv));

    doubles.resize(5);
    doubles[0] = f;
    doubles[1] = loglik;
    doubles[2] = loglik_indep;
    doubles[3] = loglik_sat;
    doubles[4] = 0.00;

    matrices.resize(2);
    matrices[0] = S - Shat;

  };

};

cfa_means_fml* choose_cfa_means_fml(const Rcpp::List& estimator_setup) {

  cfa_means_fml* myestimator = new cfa_means_fml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double w = estimator_setup["w"];
  int n = estimator_setup["n"];
  int p = estimator_setup["p"];

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));
  double plogpi2 = p * std::log(arma::datum::pi * 2);
  arma::mat I(p, p, arma::fill::eye);
  arma::vec nu(p, arma::fill::zeros);
  arma::vec means(p, arma::fill::zeros);

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_nu = indices[2];
  myestimator->indices_means = indices[3];
  myestimator->p = p;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->nu = nu;
  myestimator->means = means;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
