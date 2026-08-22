/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Confirmatory factor analysis (diagonally weighted least-squares)
 *
 * General CFA DWLS estimator. A mean structure is included whenever model and
 * sample mean indices are supplied. If they are omitted, the mean contribution
 * is exactly zero and the estimator reduces to the covariance-only objective.
 */

class cfa_dwls: public estimators {

public:

  int p;
  double loss, w;
  bool use_means;
  arma::uvec indices_S, indices_Shat, indices_nu, indices_means,
    diag, lower_diag;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;
  arma::vec nu, means, delta, dnu, dmeans, w_means;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);

    residuals = S-Shat;
    W_residuals = W % residuals;

    if(use_means) {
      nu = x.transparameters(indices_nu);
      means = x.transparameters(indices_means);
      delta = means-nu;
    } else {
      delta.zeros();
    }

  }

  void F(arguments_optim& x) {

    loss = w*0.5*(
      arma::accu(residuals % W_residuals)+
      arma::accu(w_means % delta % delta)
    );

    x.f += loss;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices_S) += w*arma::vectorise(W_residuals);
    x.grad.elem(indices_Shat) += w*arma::vectorise(-W_residuals);

    if(use_means) {
      arma::vec wdelta = w_means % delta;
      x.grad.elem(indices_nu) += -w*wdelta;
      x.grad.elem(indices_means) += w*wdelta;
    }

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    x.dgrad.elem(indices_S) +=
      w*arma::vectorise(W % dS-W % dShat);
    x.dgrad.elem(indices_Shat) +=
      w*arma::vectorise(W % dShat-W % dS);

    if(use_means) {
      dnu = x.dtransparameters(indices_nu);
      dmeans = x.dtransparameters(indices_means);
      arma::vec ddelta = dmeans-dnu;
      arma::vec wd_delta = w_means % ddelta;
      x.dgrad.elem(indices_nu) += -w*wd_delta;
      x.dgrad.elem(indices_means) += w*wd_delta;
    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();
    arma::mat W_residuals_indep = W % residuals_indep;

    // Preserve the current reporting convention for the independence loss.
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

    matrices.resize(use_means ? 2L : 1L);
    matrices[0] = residuals;

    if(use_means) {
      matrices[1] = delta;
    }

  }

};

cfa_dwls* choose_cfa_dwls(const Rcpp::List& estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  double w = estimator_setup["w"];
  arma::mat W = estimator_setup["W"];

  if(p < 1 || !std::isfinite(w) || w <= 0.0) {
    Rcpp::stop("cfa_dwls received invalid p or w values");
  }

  if(indices.size() != 2L && indices.size() != 4L) {
    Rcpp::stop("cfa_dwls requires two or four parameter-index objects");
  }

  if(indices[0].n_elem != static_cast<arma::uword>(p*p) ||
     indices[1].n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("The cfa_dwls covariance indices have incompatible dimensions");
  }

  if(W.n_rows != static_cast<arma::uword>(p) ||
     W.n_cols != static_cast<arma::uword>(p)) {
    Rcpp::stop("cfa_dwls requires W to be a p by p matrix");
  }

  bool use_means = false;
  arma::uvec indices_nu;
  arma::uvec indices_means;
  arma::vec w_means(p, arma::fill::zeros);

  if(indices.size() == 4L) {

    const bool has_nu = indices[2].n_elem > 0L;
    const bool has_means = indices[3].n_elem > 0L;

    if(has_nu != has_means) {
      Rcpp::stop("cfa_dwls requires both model and sample means, or neither");
    }

    if(has_nu) {

      if(indices[2].n_elem != static_cast<arma::uword>(p) ||
         indices[3].n_elem != static_cast<arma::uword>(p)) {
        Rcpp::stop("The cfa_dwls mean indices have incompatible dimensions");
      }

      if(!estimator_setup.containsElementNamed("w_means")) {
        Rcpp::stop("cfa_dwls requires w_means when a mean structure is supplied");
      }

      w_means = Rcpp::as<arma::vec>(estimator_setup["w_means"]);

      if(w_means.n_elem != static_cast<arma::uword>(p) ||
         !w_means.is_finite()) {
        Rcpp::stop("cfa_dwls requires one finite mean weight per observed variable");
      }

      use_means = true;
      indices_nu = indices[2];
      indices_means = indices[3];

    }

  }

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p+1, p*p-1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_nu = indices_nu;
  myestimator->indices_means = indices_means;
  myestimator->p = p;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->w_means = w_means;
  myestimator->use_means = use_means;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->nu.zeros(p);
  myestimator->means.zeros(p);
  myestimator->delta.zeros(p);
  myestimator->dnu.zeros(p);
  myestimator->dmeans.zeros(p);
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
