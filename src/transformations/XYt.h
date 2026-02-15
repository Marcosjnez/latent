/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

// X*Y.t() transformation:

class XYt:public transformations {

public:

  arma::uvec indices_X, indices_Y, indices_in, indices_out;
  int p, q;
  arma::mat X, Y, dX, dY, dXYt, grad_out, grad_in_X, grad_in_Y, jacob;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), q, q);
    arma::mat XYt = X * Y.t();
    x.transparameters(indices_out) = arma::vectorise(XYt);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out), p, q);

    grad_in_X = grad_out * Y;
    grad_in_Y = grad_out.t() * X;

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_Y) += arma::vectorise(grad_in_Y);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters(indices_Y), q, q);
    dXYt = X * dY.t() + dX * Y.t();
    x.dtransparameters(indices_out) = arma::vectorise(dXYt);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_out), p, q);

    arma::mat dgrad_in_X = dgrad_out * Y + grad_out * dY;
    arma::mat dgrad_in_Y = dgrad_out.t() * X + grad_out.t() * dX;

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);

  }

  void jacobian(arguments_optim& x) {

    arma::mat I_p = arma::eye(p, p);
    arma::mat I_q = arma::eye(q, q);

    arma::mat Kqq = dxt(q, q);

    arma::mat Jx = arma::kron(Y, I_p);
    arma::mat Jy = arma::kron(I_q, X) * Kqq;

    jacob = arma::join_rows(Jx, Jy);

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(indices_X, indices_Y);
    x.vcov(indices_out, indices_out) = jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    vectors.resize(1);

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

XYt* choose_XYt(const Rcpp::List& trans_setup) {

  XYt* mytrans = new XYt();

  arma::uvec indices_X = trans_setup["indices_X"];
  arma::uvec indices_Y = trans_setup["indices_Y"];
  arma::uvec indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
