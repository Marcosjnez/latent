/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/04/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_fml: public estimators {

public:

  int p, n;
  double f, w, logdetS, logdetShat, plogpi2;
  arma::uvec indices_S, indices_Shat, lower_diag;
  arma::mat S, Shat, residuals, dShat, dS, Shat_inv, S_inv, gShat, gS, I;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);

    if(!Shat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Shat);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      Shat = eigvec * arma::diagmat(d) * eigvec.t();
    }

    if(!S.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, S);
      arma::vec d = arma::clamp(eigval, 0.00001, eigval.max());
      S = eigvec * arma::diagmat(d) * eigvec.t();
    }

    Shat_inv = arma::inv_sympd(Shat);

  }

  void F(arguments_optim& x) {

    logdetS = arma::log_det_sympd(S);
    logdetShat = arma::log_det_sympd(Shat);

    f = w*(logdetShat - logdetS + arma::accu(S % Shat_inv) - p);
    x.f += f;

  }

  void G(arguments_optim& x) {

    S_inv = arma::inv_sympd(S);
    gS = Shat_inv - S_inv;
    gShat = Shat_inv * (I - S * Shat_inv);

    x.grad.elem(indices_S) += w*arma::vectorise(gS);
    x.grad.elem(indices_Shat) += w*arma::vectorise(gShat);

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    arma::mat dS_inv = -S_inv * dS * S_inv;
    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgShat = dShat_inv * (I - S * Shat_inv)
      - Shat_inv * dS * Shat_inv
    - Shat_inv * S * dShat_inv;

    x.dgrad.elem(indices_S) += w*arma::vectorise(dShat_inv - dS_inv);
    x.dgrad.elem(indices_Shat) += w*arma::vectorise(dgShat);

  }

  void outcomes(arguments_optim& x) {

    double loss = w*(logdetShat - logdetS + arma::accu(S % Shat_inv) - p);
    double loglik = n*0.5*(-plogpi2 - logdetShat - arma::accu(S % Shat_inv));
    arma::mat I(p, p, arma::fill::eye);
    double loglik_indep = n*0.5*(-plogpi2 - arma::trace(S));
    double loglik_sat = n*0.5*(-plogpi2 - logdetS - p);

    doubles.resize(5);
    doubles[0] =  loss;          // loss   actual model
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

cfa_fml* choose_cfa_fml(const Rcpp::List& estimator_setup) {

  cfa_fml* myestimator = new cfa_fml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double w = estimator_setup["w"];
  int n = estimator_setup["n"];
  int p = estimator_setup["p"];

  arma::mat Shat(p, p, arma::fill::zeros);
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
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
