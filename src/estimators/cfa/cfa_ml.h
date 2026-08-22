/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis (negative log-likelihood)
 *
 * General CFA likelihood estimator. The estimator can include a mean structure
 * through the model-implied means (nu) and sample means. If the mean vectors are
 * omitted from the transformed-parameter interface, their contribution is zero,
 * reproducing the covariance-only objective.
 */

class cfa_ml: public estimators {

public:

  int p, n;
  double loss, w, plogpi2;
  bool use_means;
  arma::uvec indices_Shat, indices_S, indices_nu, indices_means;
  arma::mat S, Shat, dS, dShat, Shat_inv, gS, gShat;
  arma::vec means, nu, delta, dmeans, dnu, gmeans, gnu;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);

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

  }

  void F(arguments_optim& x) {

    loss = n*0.5*(plogpi2+arma::log_det_sympd(Shat)+
      arma::accu(S % Shat_inv)+
      arma::as_scalar(delta.t()*Shat_inv*delta));

    x.f += loss;

  }

  void G(arguments_optim& x) {

    gS = Shat_inv;
    gShat = Shat_inv-
      Shat_inv*S*Shat_inv-
      Shat_inv*delta*delta.t()*Shat_inv;

    x.grad.elem(indices_S) += n*0.5*arma::vectorise(gS);
    x.grad.elem(indices_Shat) += n*0.5*arma::vectorise(gShat);

    if(use_means) {
      gnu = -2.0*Shat_inv*delta;
      gmeans = 2.0*Shat_inv*delta;
      x.grad.elem(indices_nu) += n*0.5*gnu;
      x.grad.elem(indices_means) += n*0.5*gmeans;
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

    arma::mat dgS = dShat_inv;
    arma::mat dgShat =
      dShat_inv
      -dShat_inv*S*Shat_inv
      -Shat_inv*dS*Shat_inv
      -Shat_inv*S*dShat_inv
      -dShat_inv*delta*delta.t()*Shat_inv
      -Shat_inv*dD*Shat_inv
      -Shat_inv*delta*delta.t()*dShat_inv;

    x.dgrad.elem(indices_S) += n*0.5*arma::vectorise(dgS);
    x.dgrad.elem(indices_Shat) += n*0.5*arma::vectorise(dgShat);

    if(use_means) {
      arma::vec dgnu = -2.0*(dShat_inv*delta+Shat_inv*ddelta);
      arma::vec dgmeans = 2.0*(dShat_inv*delta+Shat_inv*ddelta);
      x.dgrad.elem(indices_nu) += n*0.5*dgnu;
      x.dgrad.elem(indices_means) += n*0.5*dgmeans;
    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    double loglik = -loss;

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = loglik;

    matrices.resize(use_means ? 2L : 1L);
    matrices[0] = S-Shat;

    if(use_means) {
      matrices[1] = means-nu;
    }

  }

};

cfa_ml* choose_cfa_ml(const Rcpp::List& estimator_setup) {

  cfa_ml* myestimator = new cfa_ml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  int n = estimator_setup["n"];
  double w = estimator_setup["w"];

  if(p < 1 || n < 1 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_ml received invalid p, n, or w values");
  }

  if(indices.size() != 2L && indices.size() != 4L) {
    Rcpp::stop("cfa_ml requires two or four parameter-index objects");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("The cfa_ml covariance indices have incompatible dimensions");
  }

  bool use_means = false;
  arma::uvec indices_nu;
  arma::uvec indices_means;

  if(indices.size() == 4L) {

    const bool has_nu = indices[2].n_elem > 0L;
    const bool has_means = indices[3].n_elem > 0L;

    if(has_nu != has_means) {
      Rcpp::stop("cfa_ml requires both model and sample means, or neither");
    }

    if(has_nu) {
      if(indices[2].n_elem != static_cast<arma::uword>(p) ||
         indices[3].n_elem != static_cast<arma::uword>(p)) {
        Rcpp::stop("The cfa_ml mean indices have incompatible dimensions");
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
  myestimator->gS.zeros(p, p);
  myestimator->gShat.zeros(p, p);
  myestimator->nu.zeros(p);
  myestimator->means.zeros(p);
  myestimator->delta.zeros(p);
  myestimator->dnu.zeros(p);
  myestimator->dmeans.zeros(p);
  myestimator->gnu.zeros(p);
  myestimator->gmeans.zeros(p);
  myestimator->plogpi2 = p*std::log(arma::datum::pi*2.0);

  return myestimator;

}
