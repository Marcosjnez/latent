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

  arma::mat R, Rhat, residuals, dRhat, Rhat_inv, R_Ri, Ri_R_Ri, gRhat, I;
  arma::uvec diag, lower_diag;
  double w, logdetR, plogpi2;
  int p, q, n;

  void param(arguments_optim& x) {

    Rhat.elem(lower_diag) = x.transparameters(indices[0]);
    Rhat = arma::symmatl(Rhat);

    if(!Rhat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Rhat);
      arma::vec d = arma::clamp(eigval, 0.1, eigval.max());
      Rhat = eigvec * arma::diagmat(d) * eigvec.t();
    }
    Rhat_inv = arma::inv_sympd(Rhat);

  }

  void F(arguments_optim& x) {

    // f = w*(arma::log_det_sympd(Rhat) - logdetR + arma::accu(R % Rhat_inv) - p);
    f = w*n*0.5*(plogpi2 +
      arma::log_det_sympd(Rhat) +
      arma::accu(R % Rhat_inv));
    x.f += f;

  }

  void G(arguments_optim& x) {

    // Ri_R_Ri = Rhat_inv * R * Rhat_inv;
    arma::mat R_Ri = R * Rhat_inv;
    gRhat = Rhat_inv * (I - R_Ri);
    arma::mat temp = 2*gRhat;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += w*n*0.5*arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dRhat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dRhat = arma::symmatl(dRhat);

    arma::mat dRhat_inv = -Rhat_inv * dRhat * Rhat_inv;
    arma::mat dgRhat = 2*(dRhat_inv * (I - R * Rhat_inv) -
                          Rhat_inv * R * dRhat_inv);
    dgRhat.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += w*n*0.5*arma::vectorise(dgRhat(lower_diag));

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    // loglik = w*n*0.5*(-plogpi2 -
    //                    arma::log_det_sympd(Rhat) -
    //                    arma::accu(R % Rhat_inv));
    double loglik_indep = w*n*0.5*(-plogpi2 -
                                   arma::trace(R));
    arma::mat Rinv = arma::inv_sympd(R);
    double loglik_sat = w*n*0.5*(-plogpi2 -
                                 arma::log_det_sympd(R) -
                                 arma::accu(R % Rinv));
    doubles[0] =  f;             // loss   actual model
    doubles[1] = -f;             // loglik actual model
    doubles[2] =  w;
    doubles[3] =  loglik_indep;  // loglik independence model
    doubles[4] =  loglik_sat;    // loglik saturated model

    matrices.resize(2);
    arma::mat W;
    matrices[0] = R - Rhat;
    matrices[1] = W;

  };

};

cfa_ml* choose_cfa_ml(const Rcpp::List& estimator_setup) {

  cfa_ml* myestimator = new cfa_ml();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat R = estimator_setup["R"];
  double w = estimator_setup["w"];
  int q = estimator_setup["q"];
  int n = estimator_setup["n"];

  int p = R.n_rows;
  double logdetR = arma::log_det_sympd(R);
  arma::mat Rhat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));
  double plogpi2 = p*std::log(arma::datum::pi*2);
  arma::mat I(p, p, arma::fill::eye);

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->n = n;
  myestimator->R = R;
  myestimator->w = w;
  myestimator->logdetR = logdetR;
  myestimator->Rhat = Rhat;
  myestimator->dRhat = Rhat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;
  myestimator->I = I;

  return myestimator;

}
