/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/11/2025
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_ml: public estimators {

public:

  arma::mat S, Shat, residuals, dShat, Shat_inv, R_Ri, Ri_R_Ri, gShat, I;
  arma::uvec diag, lower_diag;
  double w, logdetS, plogpi2;
  int p, q, n;

  void param(arguments_optim& x) {

    Shat.elem(lower_diag) = x.transparameters(indices[0]);
    Shat = arma::symmatl(Shat);

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

    // f = w*(arma::log_det_sympd(Shat) - logdetR + arma::accu(S % Shat_inv) - p);
    f = w*n*0.5*(plogpi2 +
      arma::log_det_sympd(Shat) +
      arma::accu(S % Shat_inv));
    x.f += f/n;

  }

  void G(arguments_optim& x) {

    // Ri_R_Ri = Shat_inv * S * Shat_inv;
    arma::mat R_Ri = S * Shat_inv;
    gShat = Shat_inv * (I - R_Ri);
    arma::mat temp = 2*gShat;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += w*n*0.5*arma::vectorise(temp(lower_diag))/n;

  }

  void dG(arguments_optim& x) {

    dShat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dShat = arma::symmatl(dShat);

    arma::mat dShat_inv = -Shat_inv * dShat * Shat_inv;
    arma::mat dgShat = 2*(dShat_inv * (I - S * Shat_inv) -
                          Shat_inv * S * dShat_inv);
    dgShat.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += w*n*0.5*arma::vectorise(dgShat(lower_diag))/n;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    // loglik = w*n*0.5*(-plogpi2 -
    //                    arma::log_det_sympd(Shat) -
    //                    arma::accu(S % Shat_inv));
    double loglik_indep = w*n*0.5*(-plogpi2 -
                                   arma::trace(S));
    arma::mat Rinv = arma::inv_sympd(S);
    double loglik_sat = w*n*0.5*(-plogpi2 -
                                 arma::log_det_sympd(S) -
                                 arma::accu(S % Rinv));
    doubles[0] =  f;             // loss   actual model
    doubles[1] = -f;             // loglik actual model
    doubles[2] =  w;
    doubles[3] =  loglik_indep;  // loglik independence model
    doubles[4] =  loglik_sat;    // loglik saturated model

    matrices.resize(2);
    arma::mat W;
    matrices[0] = S - Shat;
    matrices[1] = W;

  };

};

cfa_ml* choose_cfa_ml(const Rcpp::List& estimator_setup) {

  cfa_ml* myestimator = new cfa_ml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat S = estimator_setup["R"];
  double w = estimator_setup["w"];
  int q = estimator_setup["q"];
  int n = estimator_setup["n"];

  int p = S.n_rows;
  double logdetS = arma::log_det_sympd(S);
  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(S));
  double plogpi2 = p*std::log(arma::datum::pi*2);
  arma::mat I(p, p, arma::fill::eye);

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->n = n;
  myestimator->S = S;
  myestimator->w = w;
  myestimator->logdetS = logdetS;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
