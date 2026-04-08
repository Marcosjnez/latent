/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/04/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_ml: public estimators {

public:

  int p, n;
  double f, w, logdetS, plogpi2;
  arma::uvec indices_S, indices_Shat, diag, lower_diag;
  arma::mat S, Shat, residuals, dS, dShat, Shat_inv, gShat, I;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);

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

    f = w*n*0.5*(plogpi2 +
      arma::log_det_sympd(Shat) +
      arma::accu(S % Shat_inv));
    x.f += f;

  }

  void G(arguments_optim& x) {

    gShat = Shat_inv * (I - S * Shat_inv);

    x.grad.elem(indices_S) += w*n*0.5*arma::vectorise(Shat_inv);
    x.grad.elem(indices_Shat) += w*n*0.5*arma::vectorise(gShat);

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgShat = dShat_inv * (I - S * Shat_inv) - Shat_inv * S * dShat_inv -
      Shat_inv * dS * Shat_inv;

    x.dgrad.elem(indices_S) += w*n*0.5*arma::vectorise(dShat_inv);
    x.dgrad.elem(indices_Shat) += w*n*0.5*arma::vectorise(dgShat);

  }

  void outcomes(arguments_optim& x) {

    double loglik = n*0.5*(-plogpi2 -
                           arma::log_det_sympd(Shat) -
                           arma::accu(S % Shat_inv));
    double loglik_indep = n*0.5*(-plogpi2 -
                                 arma::trace(S));
    arma::mat Rinv = arma::inv_sympd(S);
    double loglik_sat = n*0.5*(-plogpi2 -
                               arma::log_det_sympd(S) -
                               arma::accu(S % Rinv));

    doubles.resize(5);
    doubles[0] =  f;             // loss   actual model
    doubles[1] =  loglik;        // loglik actual model
    doubles[2] =  loglik_indep;  // loglik independence model
    doubles[3] =  loglik_sat;    // loglik saturated model
    doubles[4] =  0.00;          // penalty

    arma::mat W;
    matrices.resize(2);
    matrices[0] = S - Shat;
    matrices[1] = W;

  };

};

cfa_ml* choose_cfa_ml(const Rcpp::List& estimator_setup) {

  cfa_ml* myestimator = new cfa_ml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double w = estimator_setup["w"];
  int p = estimator_setup["p"];
  int n = estimator_setup["n"];

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));
  double plogpi2 = p*std::log(arma::datum::pi*2);
  arma::mat I(p, p, arma::fill::eye);

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->p = p;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
