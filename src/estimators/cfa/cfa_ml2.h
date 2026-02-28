/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_ml2: public estimators {

public:

  int p, q, n;
  double w, logdetR, plogpi2;
  arma::uvec indices, lower_diag;
  arma::mat R, Rhat, residuals, dRhat, Rhat_inv, Ri_R_Ri, gRhat;

  void param(arguments_optim& x) {

    Rhat.elem(lower_diag) = x.transparameters(indices);
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

    f = w*(arma::log_det_sympd(Rhat) - logdetR + arma::accu(R % Rhat_inv) - p);
    x.f += f;

  }

  void G(arguments_optim& x) {

    Ri_R_Ri = Rhat_inv * R * Rhat_inv;
    gRhat = Rhat_inv - Ri_R_Ri;
    arma::mat temp = 2*gRhat;
    temp.diag() *= 0.5;

    x.grad.elem(indices) += w*arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dRhat.elem(lower_diag) = x.dtransparameters(indices);
    dRhat = arma::symmatl(dRhat);

    arma::mat dRhat_inv = -Rhat_inv * dRhat * Rhat_inv;
    arma::mat dgRhat = 2*(dRhat_inv -
                         (dRhat_inv * R * Rhat_inv +
                          Rhat_inv * R * dRhat_inv));
    dgRhat.diag() *= 0.5;

    x.dgrad.elem(indices) += w*arma::vectorise(dgRhat(lower_diag));

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    double loglik = w*n*0.5*(-plogpi2 -
                       arma::log_det_sympd(Rhat) -
                       arma::accu(R % Rhat_inv));
    arma::mat I(p, p, arma::fill::eye);
    double loglik_indep = w*n*0.5*(-plogpi2 -
                                   arma::trace(R));
    arma::mat Rinv = arma::inv_sympd(R);
    double loglik_sat = w*n*0.5*(-plogpi2 -
                                 arma::log_det_sympd(R) -
                                 arma::accu(R % Rinv));
    doubles[0] =  f;             // loss   actual model
    doubles[1] =  loglik;        // loglik actual model
    doubles[2] =  w;
    doubles[3] =  loglik_indep;  // loglik independence model
    doubles[4] =  loglik_sat;    // loglik saturated model

    matrices.resize(2);
    arma::mat W;
    matrices[0] = R - Rhat;
    matrices[1] = W;

  };

};

cfa_ml2* choose_cfa_ml2(const Rcpp::List& estimator_setup) {

  cfa_ml2* myestimator = new cfa_ml2();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat R = estimator_setup["R"];
  double w = estimator_setup["w"];
  int q = estimator_setup["q"];
  int n = estimator_setup["n"];

  int p = R.n_rows;
  double logdetR = arma::log_det_sympd(R);
  arma::mat Rhat(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));
  double plogpi2 = p*std::log(arma::datum::pi*2);

  myestimator->indices = indices[0];
  myestimator->p = p;
  myestimator->q = q;
  myestimator->n = n;
  myestimator->R = R;
  myestimator->w = w;
  myestimator->logdetR = logdetR;
  myestimator->Rhat = Rhat;
  myestimator->dRhat = Rhat;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;

  return myestimator;

}
