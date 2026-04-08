/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 08/04/2026
 */

/*
 * Confirmatory factor analysis with mean structure (weighted least-squares)
 */

class cfa_means_dwls: public estimators {

public:

  int p;
  double f, w;
  arma::uvec indices_S, indices_Shat, indices_nu, indices_means, diag, lower_diag;
  arma::mat S, Shat, dShat, dS, residuals, W, W_residuals;
  arma::vec nu, means, delta, dnu, dmeans, w_means;

  void param(arguments_optim& x) {

    S = arma::reshape(x.transparameters(indices_S), p, p);
    Shat = arma::reshape(x.transparameters(indices_Shat), p, p);
    nu = x.transparameters(indices_nu);
    means = x.transparameters(indices_means);

    residuals = S - Shat;
    W_residuals = W % residuals;
    delta = means - nu;

  }

  void F(arguments_optim& x) {

    f = w*0.5*(arma::accu(residuals % W_residuals) + arma::accu(w_means % delta % delta));
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::vec wdelta = w_means % delta;

    x.grad.elem(indices_S) += w*arma::vectorise(W_residuals);
    x.grad.elem(indices_Shat) += w*arma::vectorise(-W_residuals);
    x.grad.elem(indices_nu) += -w*wdelta;
    x.grad.elem(indices_means) += w*wdelta;

  }

  void dG(arguments_optim& x) {

    dS = arma::reshape(x.dtransparameters(indices_S), p, p);
    dShat = arma::reshape(x.dtransparameters(indices_Shat), p, p);
    dnu = x.dtransparameters(indices_nu);
    dmeans = x.dtransparameters(indices_means);

    arma::vec ddelta = dmeans - dnu;
    arma::vec wd_delta = w_means % ddelta;

    x.dgrad.elem(indices_S) += w*arma::vectorise(W % dS - W % dShat);
    x.dgrad.elem(indices_Shat) += w*arma::vectorise(W % dShat - W % dS);
    x.dgrad.elem(indices_nu) += -w*wd_delta;
    x.dgrad.elem(indices_means) += w*wd_delta;

  }

  void outcomes(arguments_optim& x) {

    arma::mat residuals_indep = S;
    residuals_indep.diag().zeros();
    arma::mat W_residuals_indep = W % residuals_indep;
    double loss_indep = w*0.5*arma::accu(residuals_indep % W_residuals_indep);
    double loss_sat = 0.00;

    doubles.resize(5);
    doubles[0] =  f;
    doubles[1] =  0.00;
    doubles[2] =  loss_indep;
    doubles[3] =  loss_sat;
    doubles[4] =  0.00;

    matrices.resize(2);
    matrices[0] = residuals;
    matrices[1] = arma::reshape(means - nu, p, 1);

  };

};

cfa_means_dwls* choose_cfa_means_dwls(const Rcpp::List& estimator_setup) {

  cfa_means_dwls* myestimator = new cfa_means_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  int p = estimator_setup["p"];
  double w = estimator_setup["w"];
  arma::mat W = estimator_setup["W"];
  arma::vec w_means = estimator_setup["w_means"];

  arma::mat Shat(p, p, arma::fill::zeros);
  arma::uvec diag = arma::regspace<arma::uvec>(0, p + 1, p*p - 1);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(Shat));
  arma::vec nu(p, arma::fill::zeros);
  arma::vec means(p, arma::fill::zeros);

  myestimator->indices_Shat = indices[0];
  myestimator->indices_S = indices[1];
  myestimator->indices_nu = indices[2];
  myestimator->indices_means = indices[3];
  myestimator->p = p;
  myestimator->W = W;
  myestimator->w = w;
  myestimator->w_means = w_means;
  myestimator->Shat = Shat;
  myestimator->dShat = Shat;
  myestimator->nu = nu;
  myestimator->means = means;
  myestimator->diag = diag;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
