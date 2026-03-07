/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/03/2026
 */

/*
 * Confirmatory factor analysis (weighted least-squares)
 */

class cfa_dwls_R: public estimators {

public:

  int p;
  double w;
  arma::uvec indices_S, indices_Shat, diag, lower_diag;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);

    residuals = S - Shat;
    W_residuals = W % residuals;

  }

  void F(arguments_optim& x) {

    f = w*0.5*arma::accu(residuals % W_residuals);
    x.f += f;

  }

  void G(arguments_optim& x) {

    x.grad.elem(indices_S) += w*arma::vectorise(W_residuals);
    x.grad.elem(indices_Shat) += w*arma::vectorise(-W_residuals);

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);

    x.dgrad.elem(indices_S) += w*arma::vectorise(W % dS - W % dShat);
    x.dgrad.elem(indices_Shat) += w*arma::vectorise(W % dShat - W % dS);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(5);
    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();
    arma::mat W_residuals_indep = W % residuals_indep;
    double loss_indep = w*0.5*arma::accu(residuals_indep % W_residuals_indep);
    double loss_sat = 0.00;
    doubles[0] =  f;
    doubles[1] =  0.00;        // loglik actual model
    doubles[2] =  w;
    doubles[3] =  loss_indep;  // loglik independence model
    doubles[4] =  loss_sat;    // loglik saturated model

    matrices.resize(2);
    matrices[0] = residuals;
    matrices[1] = W;

  };

};

cfa_dwls_R* choose_cfa_dwls_R(const Rcpp::List& estimator_setup) {

  cfa_dwls_R* myestimator = new cfa_dwls_R();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  double w = estimator_setup["w"];
  arma::mat W = estimator_setup["W"];

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->p = p;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
