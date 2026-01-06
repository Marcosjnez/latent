/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/12/2025
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_ml_R: public estimators {

public:

  arma::mat R, Rhat, residuals, dR, dRhat, Rhat_inv, Ri_R_Ri, gRhat;
  arma::uvec lower_diag;
  double w, logdetR, plogpi2;
  int p, q, n;

  void param(arguments_optim& x) {

    R.elem(lower_diag) = x.transparameters(indices[0]);
    R = arma::symmatl(R);
    Rhat.elem(lower_diag) = x.transparameters(indices[1]);
    Rhat = arma::symmatl(Rhat);

    if(!Rhat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Rhat);
      arma::vec d = arma::clamp(eigval, 0.1, eigval.max());
      Rhat = eigvec * arma::diagmat(d) * eigvec.t();
    }
    Rhat_inv = arma::inv_sympd(Rhat);
    // logdetR = arma::log_det_sympd(R);

  }

  void F(arguments_optim& x) {

    // f = w*(arma::log_det_sympd(Rhat) - logdetR + arma::accu(R % Rhat_inv) - p);
    f = w*n*0.5*(plogpi2 +
      arma::log_det_sympd(Rhat) +
      arma::accu(R % Rhat_inv));
    x.f += f;

  }

  void G(arguments_optim& x) {

    Ri_R_Ri = Rhat_inv * R * Rhat_inv;
    gRhat = Rhat_inv - Ri_R_Ri;
    arma::mat temp = 2*gRhat;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += w*n*0.5*arma::vectorise(Rhat_inv(lower_diag));
    x.grad.elem(indices[1]) += w*n*0.5*arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dR.elem(lower_diag) = x.dtransparameters(indices[0]);
    dR = arma::symmatl(dR);
    dRhat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dRhat = arma::symmatl(dRhat);

    arma::mat dRhat_inv = -Rhat_inv * dRhat * Rhat_inv;
    arma::mat dgRhat = 2*(dRhat_inv -
                         (dRhat_inv * R * Rhat_inv +
                          Rhat_inv * R * dRhat_inv +
                          Rhat_inv * dR * Rhat_inv));
    dgRhat.diag() *= 0.5;

    arma::mat dgR = 2*(dRhat_inv * R * Rhat_inv +
                       Rhat_inv * dR * Rhat_inv +
                       Rhat_inv * R * dRhat_inv);
    dgR.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += w*n*0.5*arma::vectorise(dgR(lower_diag));
    x.dgrad.elem(indices[1]) += w*n*0.5*arma::vectorise(dgRhat(lower_diag));

  }

  void H(arguments_optim& x) {

    // FIX THIS
    arma::mat hx = -(arma::kron(Rhat_inv, Rhat_inv) -
      (arma::kron(Ri_R_Ri, Rhat_inv) + arma::kron(Rhat_inv, Ri_R_Ri)));
    hx.rows(lower_diag) *= 2;
    hx = arma::diagmat(arma::vectorise(hx.elem(lower_diag)));

    x.hess(indices[0], indices[0]) += w*n*0.5*hx;
    x.hess(indices[1], indices[1]) += w*n*0.5*hx;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    // loglik = w*n*0.5*(-plogpi2 -
    //                    arma::log_det_sympd(Rhat) -
    //                    arma::accu(R % Rhat_inv));
    arma::mat I(p, p, arma::fill::eye);
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

cfa_ml_R* choose_cfa_ml_R(const Rcpp::List& estimator_setup) {

  cfa_ml_R* myestimator = new cfa_ml_R();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  double w = estimator_setup["w"];
  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  int n = estimator_setup["n"];

  arma::mat R(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));
  double plogpi2 = p*std::log(arma::datum::pi*2);

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->n = n;
  myestimator->w = w;
  myestimator->R = R;
  myestimator->dR = R;
  myestimator->Rhat = R;
  myestimator->dRhat = R;
  myestimator->lower_diag = lower_diag;
  myestimator->plogpi2 = plogpi2;

  return myestimator;

}
