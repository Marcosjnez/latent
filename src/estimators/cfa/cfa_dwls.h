/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/11/2025
 */

/*
 * Confirmatory factor analysis (least-squares, weigthed and unweighted)
 */

class cfa_dwls: public estimators {

public:

  arma::mat S, Shat, dShat, residuals, W, W_residuals;
  arma::uvec diag, lower_diag;
  int p, q;
  double w;

  void param(arguments_optim& x) {

    Shat.elem(lower_diag) = x.transparameters(indices[0]);
    Shat = arma::symmatl(Shat);

    residuals = S - Shat;
    W_residuals = W % residuals;

  }

  void F(arguments_optim& x) {

    f = w*0.5*arma::accu(residuals % W_residuals);
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat temp = -2*W_residuals;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += w*arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dShat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dShat = arma::symmatl(dShat);

    arma::mat dgShat = 2*W % dShat;
    dgShat.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += w*arma::vectorise(dgShat(lower_diag));

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

cfa_dwls* choose_cfa_dwls(const Rcpp::List& estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat S = estimator_setup["R"];
  arma::mat W = estimator_setup["W"];
  double w = estimator_setup["w"];
  int q = estimator_setup["q"];

  int p = S.n_rows;
  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(S));

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->S = S;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
