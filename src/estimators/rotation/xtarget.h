/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 15/02/2026
 */

/*
 * Extended target
 */

class xtarget: public estimators {

public:

  int p, q;
  double w;
  arma::uvec indices_lambda, indices_psi, lower_psi;
  arma::mat lambda, dlambda, target, weight, psi, dpsi, psitarget, psiweight,
  weight2, psiweight2, f1, f2;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);
    psi.elem(lower_psi) = x.transparameters(indices_psi);
    psi = arma::symmatl(psi);

    f1 = weight % (lambda - target);
    f2 = psiweight % (psi - psitarget);

  }

  void F(arguments_optim& x) {

    f = 0.5*arma::accu(f1 % f1) + 0.5*w*arma::accu(f2 % f2);
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat df_dlambda = weight % f1;
    arma::mat df_dpsi = 2*w * psiweight % f2;
    df_dpsi.diag() *= 0.5;

    x.grad.elem(indices_lambda) += arma::vectorise(df_dlambda);
    x.grad.elem(indices_psi) += df_dpsi.elem(lower_psi);

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    dpsi.elem(lower_psi) = x.dtransparameters(indices_psi);
    dpsi = arma::symmatl(dpsi);

    arma::mat ddf_dlambda = weight2 % dlambda;
    arma::mat ddf_dpsi = 2*w * psiweight2 % dpsi;
    ddf_dpsi.diag() *= 0.5;

    x.dgrad.elem(indices_lambda) += arma::vectorise(ddf_dlambda);
    x.dgrad.elem(indices_psi) += ddf_dpsi.elem(lower_psi);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(2);
    doubles[0] =  f;
    doubles[0] =  0.00;

  }

};

xtarget* choose_xtarget(const Rcpp::List& estimator_setup) {

  xtarget* myestimator = new xtarget();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat target = estimator_setup["target"];
  arma::mat weight = estimator_setup["weight"];
  arma::mat psitarget = estimator_setup["psitarget"];
  arma::mat psiweight = estimator_setup["psiweight"];
  double w = estimator_setup["w"];

  int p = target.n_rows;
  int q = target.n_cols;
  arma::mat psi(q, q, arma::fill::zeros);
  arma::uvec lower_psi = arma::trimatl_ind(arma::size(psi));
  arma::mat weight2 = weight % weight;
  arma::mat psiweight2 = psiweight % psiweight;

  myestimator->indices_lambda = indices[0];
  myestimator->indices_psi = indices[1];
  myestimator->psi = psi;
  myestimator->dpsi = psi;
  myestimator->target = target;
  myestimator->weight = weight;
  myestimator->weight2 = weight2;
  myestimator->psitarget = psitarget;
  myestimator->psiweight = psiweight;
  myestimator->psiweight2 = psiweight2;
  myestimator->w = w;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->lower_psi = lower_psi;

  return myestimator;

}
