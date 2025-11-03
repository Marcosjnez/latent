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

  arma::mat R, Rhat, residuals, dRhat, Rhat_inv, Ri_R_Ri, gRhat;
  arma::uvec lower_diag;
  double w, logdetR, loglik;
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

    f = w*(arma::log_det_sympd(Rhat) - logdetR + arma::accu(R % Rhat_inv) - p);
    x.f += f;

  }

  void G(arguments_optim& x) {

    Ri_R_Ri = Rhat_inv * R * Rhat_inv;
    gRhat = Rhat_inv - Ri_R_Ri;
    arma::mat temp = 2*gRhat;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += w*arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dRhat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dRhat = arma::symmatl(dRhat);

    arma::mat dRhat_inv = -Rhat_inv * dRhat * Rhat_inv;
    arma::mat dgRhat = 2*(dRhat_inv -
                         (dRhat_inv * R * Rhat_inv +
                          Rhat_inv * R * dRhat_inv));
    dgRhat.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += w*arma::vectorise(dgRhat(lower_diag));

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

    arma::mat hx = -(arma::kron(Rhat_inv, Rhat_inv) -
      (arma::kron(Ri_R_Ri, Rhat_inv) + arma::kron(Rhat_inv, Ri_R_Ri)));
    hx.rows(lower_diag) *= 2;
    hx = arma::diagmat(arma::vectorise(hx.elem(lower_diag)));

    x.hess(indices[0], indices[0]) += w*hx;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(3);
    loglik = -w*(0.5*n*p*std::log(2*arma::datum::pi) -
                 0.5*n*arma::log_det_sympd(Rhat) -
                 0.5*n*arma::trace(R*Rhat_inv));
    doubles[0] = f;
    doubles[1] = loglik;
    doubles[2] = w;

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
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->n = n;
  myestimator->R = R;
  myestimator->w = w;
  myestimator->logdetR = logdetR;
  myestimator->Rhat = Rhat;
  myestimator->dRhat = Rhat;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
