/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/03/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_fml: public estimators {

public:

  int p, n;
  double w, logdetS, plogpi2;
  arma::uvec indices, lower_diag;
  arma::mat S, Shat, residuals, dShat, Shat_inv, gShat, I;

  void param(arguments_optim& x) {

    Shat = arma::reshape(x.transparameters(indices), p, p);

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

    f = w*(arma::log_det_sympd(Shat) - logdetS + arma::accu(S % Shat_inv) - p);
    x.f += f;

  }

  void G(arguments_optim& x) {

    gShat = Shat_inv * (I - S * Shat_inv);

    x.grad.elem(indices) += w*arma::vectorise(gShat);

  }

  void dG(arguments_optim& x) {

    dShat = arma::reshape(x.dtransparameters(indices), p, p);

    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgShat = dShat_inv * (I - S * Shat_inv) - Shat_inv * S * dShat_inv;

    x.dgrad.elem(indices) += w*arma::vectorise(dgShat);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    double loglik = w*n*0.5*(-plogpi2 -
                       arma::log_det_sympd(Shat) -
                       arma::accu(S % Shat_inv));
    arma::mat I(p, p, arma::fill::eye);
    double loglik_indep = w*n*0.5*(-plogpi2 -
                                   arma::trace(S));
    arma::mat Sinv = arma::inv_sympd(S);
    double loglik_sat = w*n*0.5*(-plogpi2 -
                                 arma::log_det_sympd(S) -
                                 arma::accu(S % Sinv));
    doubles[0] =  f;             // loss   actual model
    doubles[1] =  loglik;        // loglik actual model
    doubles[2] =  w;             // weight scalar
    doubles[3] =  loglik_indep;  // loglik independence model
    doubles[4] =  loglik_sat;    // loglik saturated model

    matrices.resize(2);
    arma::mat W;
    matrices[0] = S - Shat;
    matrices[1] = W;

  };

};

cfa_fml* choose_cfa_fml(const Rcpp::List& estimator_setup) {

  cfa_fml* myestimator = new cfa_fml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat S = estimator_setup["R"];
  double w = estimator_setup["w"];
  int n = estimator_setup["n"];

  int p = S.n_rows;
  double logdetS = arma::log_det_sympd(S);
  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(S));
  double plogpi2 = p*std::log(arma::datum::pi*2);
  arma::mat I(p, p, arma::fill::eye);

  myestimator->indices = indices[0];
  myestimator->p = p;
  myestimator->n = n;
  myestimator->S = S;
  myestimator->w = w;
  myestimator->logdetS = logdetS;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
