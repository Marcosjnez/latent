/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 30/10/2025
 */

/*
 * Confirmatory factor analysis (least-squares, weigthed and unweighted)
 */

class cfa_dwls: public estimators {

public:

  arma::mat R, Rhat, dRhat, residuals, W, W_residuals;
  arma::uvec lower_diag;
  int p, q;

  void param(arguments_optim& x) {

    Rhat.elem(lower_diag) = x.transparameters(indices[0]);
    Rhat = arma::symmatl(Rhat);

    residuals = R - Rhat;
    W_residuals = W % residuals;

  }

  void F(arguments_optim& x) {

    f = 0.5*arma::accu(residuals % W_residuals);
    x.f += f;

  }

  void G(arguments_optim& x) {

    arma::mat temp = -2*W_residuals;
    temp.diag() *= 0.5;

    x.grad.elem(indices[0]) += arma::vectorise(temp(lower_diag));

  }

  void dG(arguments_optim& x) {

    dRhat.elem(lower_diag) = x.dtransparameters(indices[0]);
    dRhat = arma::symmatl(dRhat);

    arma::mat dgRhat = 2*W % dRhat;
    dgRhat.diag() *= 0.5;

    x.dgrad.elem(indices[0]) += arma::vectorise(dgRhat(lower_diag));

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

    arma::mat hx = 2*W;
    hx.diag() *= 0.5;
    hx = arma::diagmat(arma::vectorise(hx.elem(lower_diag)));

    x.hess(indices[0], indices[0]) += hx;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = f;

    matrices.resize(2);
    matrices[0] = residuals;
    matrices[1] = W;

  };

};

cfa_dwls* choose_cfa_dwls(const Rcpp::List& estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat R = estimator_setup["R"];
  arma::mat W = estimator_setup["W"];
  int q = estimator_setup["q"];

  int p = R.n_rows;
  arma::mat Rhat(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(R));

  myestimator->indices = indices;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->R = R;
  myestimator->W = W;
  myestimator->Rhat = Rhat;
  myestimator->dRhat = Rhat;
  myestimator->lower_diag = lower_diag;

  return myestimator;

}
