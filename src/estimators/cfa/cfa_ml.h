/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 21/08/2026
 */

/*
 * Confirmatory factor analysis (negative log-likelihood)
 *
 * The sample covariance matrix and the model-implied covariance matrix use the
 * same transformed-parameter interface as cfa_fml. The estimator evaluates the
 * total negative log-likelihood contribution; therefore, the inverse Hessian
 * is already on the finite-sample information scale.
 */

class cfa_ml: public estimators {

public:

  int p, n;
  double loss, w, plogpi2;
  arma::uvec indices_S, indices_Shat;
  arma::mat S, Shat, dS, dShat, Shat_inv, gS, gShat, I;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    S = arma::reshape(x.transparameters(indices_S), p, p);

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

    loss = n*0.5*(plogpi2 + arma::log_det_sympd(Shat) +
      arma::accu(S % Shat_inv));

    x.f += loss;

  }

  void G(arguments_optim& x) {

    gS = Shat_inv;
    gShat = Shat_inv * (I - S * Shat_inv);

    x.grad.elem(indices_S) += n*0.5*arma::vectorise(gS);
    x.grad.elem(indices_Shat) += n*0.5*arma::vectorise(gShat);

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgS = dShat_inv;
    arma::mat dgShat = dShat_inv * (I - S * Shat_inv) -
      Shat_inv * dS * Shat_inv -
      Shat_inv * S * dShat_inv;

    x.dgrad.elem(indices_S) += n*0.5*arma::vectorise(dgS);
    x.dgrad.elem(indices_Shat) += n*0.5*arma::vectorise(dgShat);

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    double loglik = -loss;

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = loglik;

    matrices.resize(1);
    matrices[0] = S - Shat;

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

  if(indices.size() != 2L) {
    Rcpp::stop("cfa_ml requires model and sample covariance indices");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("The cfa_ml parameter indices have incompatible dimensions");
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
  myestimator->gS.zeros(p, p);
  myestimator->gShat.zeros(p, p);
  myestimator->I.eye(p, p);
  myestimator->plogpi2 = p*std::log(arma::datum::pi*2.0);

  return myestimator;

}
