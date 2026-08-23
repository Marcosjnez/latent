/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 23/08/2026
 */

// X.t()*Y*X symmetric congruence transformation:

class XtYX: public transformations {

public:

  int p, q;
  arma::uvec indices_X, indices_Y, indices_XtYX;
  arma::mat X, Y, dX, dY, dXtYX, grad_out, grad_in_X, grad_in_Y;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), p, p);
    arma::mat XtYX = X.t()*Y*X;
    x.transparameters(indices_XtYX) = arma::vectorise(XtYX);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_XtYX), q, q);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    grad_in_X = Y*X*grad_out.t()+Y.t()*X*grad_out;
    grad_in_Y = X*grad_out*X.t();

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_Y) += arma::vectorise(grad_in_Y);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters(indices_Y), p, p);

    dXtYX = dX.t()*Y*X+X.t()*dY*X+X.t()*Y*dX;
    x.dtransparameters(indices_XtYX) = arma::vectorise(dXtYX);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_XtYX), q, q);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    arma::mat dgrad_in_X =
      dY*X*grad_out.t()+Y*dX*grad_out.t()+Y*X*dgrad_out.t()+
      dY.t()*X*grad_out+Y.t()*dX*grad_out+Y.t()*X*dgrad_out;

    arma::mat dgrad_in_Y =
      dX*grad_out*X.t()+X*dgrad_out*X.t()+X*grad_out*dX.t();

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);

  }

  void jacobian(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), p, p);

    jacob.zeros(q*q, p*q+p*p);

    arma::mat E_X(p, q, arma::fill::zeros);
    arma::mat E_Y(p, p, arma::fill::zeros);

    for(int k = 0; k < p*q; ++k) {
      E_X.zeros();
      E_X[k] = 1.00;
      jacob.col(k) = arma::vectorise(E_X.t()*Y*X+X.t()*Y*E_X);
    }

    for(int k = 0; k < p*p; ++k) {
      E_Y.zeros();
      E_Y[k] = 1.00;
      jacob.col(p*q+k) = arma::vectorise(X.t()*E_Y*X);
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

XtYX* choose_XtYX(const Rcpp::List& trans_setup) {

  XtYX* mytrans = new XtYX();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  if(p < 1 || q < 1) {
    Rcpp::stop("XtYX requires positive p and q dimensions.");
  }

  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rcpp::stop("XtYX requires X and Y inputs and one output.");
  }

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_Y = indices_in[1];
  arma::uvec indices_XtYX = indices_out[0];

  if(indices_X.n_elem != static_cast<arma::uword>(p*q) ||
     indices_Y.n_elem != static_cast<arma::uword>(p*p) ||
     indices_XtYX.n_elem != static_cast<arma::uword>(q*q)) {
    Rcpp::stop("The XtYX parameter indices have incompatible dimensions.");
  }

  mytrans->indices_in = arma::join_cols(indices_X, indices_Y);
  mytrans->indices_out = indices_XtYX;
  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_XtYX = indices_XtYX;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
