/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 23/08/2026
 */

// X.t()*Y transformation:

class XtY: public transformations {

public:

  int p, q, r;
  arma::uvec indices_X, indices_Y, indices_XtY;
  arma::mat X, Y, dX, dY, dXtY, grad_out, grad_in_X, grad_in_Y;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), p, r);
    arma::mat XtY = X.t()*Y;
    x.transparameters(indices_XtY) = arma::vectorise(XtY);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_XtY), q, r);

    grad_in_X = Y*grad_out.t();
    grad_in_Y = X*grad_out;

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_Y) += arma::vectorise(grad_in_Y);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters(indices_Y), p, r);
    dXtY = dX.t()*Y+X.t()*dY;

    x.dtransparameters(indices_XtY) = arma::vectorise(dXtY);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_XtY), q, r);

    arma::mat dgrad_in_X = dY*grad_out.t()+Y*dgrad_out.t();
    arma::mat dgrad_in_Y = dX*grad_out+X*dgrad_out;

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);

  }

  void jacobian(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), p, r);

    jacob.zeros(q*r, p*q+p*r);

    arma::mat E_X(p, q, arma::fill::zeros);
    arma::mat E_Y(p, r, arma::fill::zeros);

    for(int k = 0; k < p*q; ++k) {
      E_X.zeros();
      E_X[k] = 1.00;
      jacob.col(k) = arma::vectorise(E_X.t()*Y);
    }

    for(int k = 0; k < p*r; ++k) {
      E_Y.zeros();
      E_Y[k] = 1.00;
      jacob.col(p*q+k) = arma::vectorise(X.t()*E_Y);
    }

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

XtY* choose_XtY(const Rcpp::List& trans_setup) {

  XtY* mytrans = new XtY();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];
  int r = trans_setup["r"];

  if(p < 1 || q < 1 || r < 1) {
    Rcpp::stop("XtY requires positive p, q, and r dimensions.");
  }

  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rcpp::stop("XtY requires X and Y inputs and one output.");
  }

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_Y = indices_in[1];
  arma::uvec indices_XtY = indices_out[0];

  if(indices_X.n_elem != static_cast<arma::uword>(p*q) ||
     indices_Y.n_elem != static_cast<arma::uword>(p*r) ||
     indices_XtY.n_elem != static_cast<arma::uword>(q*r)) {
    Rcpp::stop("The XtY parameter indices have incompatible dimensions.");
  }

  mytrans->indices_in = arma::join_cols(indices_X, indices_Y);
  mytrans->indices_out = indices_XtY;
  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_XtY = indices_XtY;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->r = r;

  return mytrans;

}
