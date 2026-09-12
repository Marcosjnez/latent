/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

// X*Y*X.t() symmetric congruence transformation:
// X is p x q, Y is q x q, and the output is p x p.

class XYXt: public transformations {

public:

  int p, q;
  arma::uvec indices_X, indices_Y, indices_XYXt;
  arma::mat X, Y, dX, dY, dXYXt, grad_out, grad_in_X, grad_in_Y;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), q, q);
    arma::mat XYXt = X*Y*X.t();
    x.transparameters(indices_XYXt) = arma::vectorise(XYXt);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_XYXt), p, p);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    grad_in_X = grad_out*X*Y.t()+grad_out.t()*X*Y;
    grad_in_Y = X.t()*grad_out*X;

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_Y) += arma::vectorise(grad_in_Y);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters(indices_Y), q, q);

    dXYXt = dX*Y*X.t()+X*dY*X.t()+X*Y*dX.t();
    x.dtransparameters(indices_XYXt) = arma::vectorise(dXYXt);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_XYXt), p, p);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    arma::mat dgrad_in_X =
      dgrad_out*X*Y.t()+grad_out*dX*Y.t()+grad_out*X*dY.t()+
      dgrad_out.t()*X*Y+grad_out.t()*dX*Y+grad_out.t()*X*dY;

    arma::mat dgrad_in_Y =
      dX.t()*grad_out*X+X.t()*dgrad_out*X+X.t()*grad_out*dX;

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);

  }

  void jacobian(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), q, q);

    jacob.zeros(p*p, p*q+q*q);

    arma::mat E_X(p, q, arma::fill::zeros);
    arma::mat E_Y(q, q, arma::fill::zeros);

    for(int k = 0; k < p*q; ++k) {
      E_X.zeros();
      E_X[k] = 1.00;
      jacob.col(k) = arma::vectorise(E_X*Y*X.t()+X*Y*E_X.t());
    }

    for(int k = 0; k < q*q; ++k) {
      E_Y.zeros();
      E_Y[k] = 1.00;
      jacob.col(p*q+k) = arma::vectorise(X*E_Y*X.t());
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

XYXt* choose_XYXt(const Rcpp::List& trans_setup) {

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  if(p < 1 || q < 1) {
    Rcpp::stop("XYXt requires positive p and q dimensions.");
  }

  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rcpp::stop("XYXt requires X and Y inputs and one output.");
  }

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_Y = indices_in[1];
  arma::uvec indices_XYXt = indices_out[0];

  if(indices_X.n_elem != static_cast<arma::uword>(p*q) ||
     indices_Y.n_elem != static_cast<arma::uword>(q*q) ||
     indices_XYXt.n_elem != static_cast<arma::uword>(p*p)) {
    Rcpp::stop("The XYXt parameter indices have incompatible dimensions.");
  }

  XYXt* mytrans = new XYXt();

  mytrans->indices_in = arma::join_cols(indices_X, indices_Y);
  mytrans->indices_out = indices_XYXt;
  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_XYXt = indices_XYXt;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
