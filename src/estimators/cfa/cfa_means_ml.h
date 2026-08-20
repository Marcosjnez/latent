/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/08/2026
 */

/*
 * Saturated multivariate-normal moments with a mean structure.
 *
 * The sample covariance matrix and sample means are fixed objects supplied in
 * estimator_setup. Only the model covariance matrix and model means are
 * parameters. The optimized criterion is the -2 log-likelihood divided by the
 * total sample size, apart from constants that do not depend on the parameters.
 */

class cfa_means_ml: public estimators {

public:

  int p, n;
  double loss, w, plogpi2;
  arma::uvec indices_Shat, indices_nu;
  arma::mat S, Shat, Shat_inv, dShat, gShat;
  arma::vec means, nu, delta, dnu, gnu;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    nu = x.transparameters(indices_nu);
    delta = means - nu;

    if(!Shat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Shat);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      Shat = eigvec * arma::diagmat(d) * eigvec.t();
    }

    Shat_inv = arma::inv_sympd(Shat);

  }

  void F(arguments_optim& x) {

    loss = w*(arma::log_det_sympd(Shat) +
      arma::accu(S % Shat_inv) +
      arma::as_scalar(delta.t() * Shat_inv * delta));

    x.f += loss;

  }

  void G(arguments_optim& x) {

    gShat = Shat_inv -
      Shat_inv * S * Shat_inv -
      Shat_inv * delta * delta.t() * Shat_inv;

    gnu = -2.0 * Shat_inv * delta;

    x.grad.elem(indices_Shat) += w * arma::vectorise(gShat);
    x.grad.elem(indices_nu) += w * gnu;

  }

  void dG(arguments_optim& x) {

    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dnu = x.dtransparameters(indices_nu);

    arma::vec ddelta = -dnu;
    arma::mat dD = ddelta * delta.t() + delta * ddelta.t();
    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;

    arma::mat dgShat =
      dShat_inv
      - dShat_inv * S * Shat_inv
      - Shat_inv * S * dShat_inv
      - dShat_inv * delta * delta.t() * Shat_inv
      - Shat_inv * dD * Shat_inv
      - Shat_inv * delta * delta.t() * dShat_inv;

    arma::vec dgnu = -2.0 *
      (dShat_inv * delta + Shat_inv * ddelta);

    x.dgrad.elem(indices_Shat) += w * arma::vectorise(dgShat);
    x.dgrad.elem(indices_nu) += w * dgnu;

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    double mean_term = arma::as_scalar(delta.t() * Shat_inv * delta);
    double unweighted_loss = arma::log_det_sympd(Shat) +
      arma::accu(S % Shat_inv) + mean_term;

    double loglik = -0.5*n*(plogpi2 + unweighted_loss);

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = loglik;

    matrices.resize(2);
    matrices[0] = S - Shat;
    matrices[1] = means - nu;

  }

};

cfa_means_ml* choose_cfa_means_ml(const Rcpp::List& estimator_setup) {

  cfa_means_ml* myestimator = new cfa_means_ml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int n = estimator_setup["n"];
  double w = estimator_setup["w"];
  arma::mat S = estimator_setup["S"];
  arma::vec means = estimator_setup["means"];

  if(p < 1 || n < 1 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_means_ml received invalid p, n, or w values");
  }

  if(indices.size() != 2L) {
    Rcpp::stop("cfa_means_ml requires covariance and mean parameter indices");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The cfa_means_ml parameter indices have incompatible dimensions");
  }

  if(S.n_rows != static_cast<arma::uword>(p) ||
     S.n_cols != static_cast<arma::uword>(p)) {
    Rcpp::stop("The fixed sample covariance has incompatible dimensions");
  }

  if(means.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The fixed sample means have incompatible dimensions");
  }

  if(!S.is_finite() || !means.is_finite() ||
     !S.is_symmetric(1e-10)) {
    Rcpp::stop("The fixed sample moments must be finite and symmetric");
  }

  myestimator->indices_Shat = indices[0];
  myestimator->indices_nu = indices[1];
  myestimator->p = p;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->S = S;
  myestimator->means = means;
  myestimator->Shat.zeros(p, p);
  myestimator->Shat_inv.zeros(p, p);
  myestimator->dShat.zeros(p, p);
  myestimator->gShat.zeros(p, p);
  myestimator->nu.zeros(p);
  myestimator->delta.zeros(p);
  myestimator->dnu.zeros(p);
  myestimator->gnu.zeros(p);
  myestimator->plogpi2 = p*std::log(arma::datum::pi*2.0);

  return(myestimator);

}
