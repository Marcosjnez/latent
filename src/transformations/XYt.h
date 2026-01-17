/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

// X*Y.t() transformation:

class XYt:public transformations {

public:

  int p, q;
  arma::mat X, Y, dX, dY, dXYt, grad_out, grad_in_X, grad_in_Y;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in[0]), p, q);
    Y = arma::reshape(x.transparameters(indices_in[1]), q, q);
    arma::mat XYt = X * Y.t();
    x.transparameters(indices_out[0]) = arma::vectorise(XYt);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out[0]), p, q);

    grad_in_X = grad_out * Y;
    grad_in_Y = X.t() * grad_out;

    x.grad(indices_in[0]) += arma::vectorise(grad_in_X);
    x.grad(indices_in[1]) += arma::vectorise(grad_in_Y);

    // arma::vec v = x.grad(indices_in[0]);
    // for (arma::uword i = 0; i < v.n_elem; ++i) {
    //   Rprintf("%.6f%s", v(i), (i + 1 < v.n_elem) ? " " : "\n"); // space-separated, then newline
    // }

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_in[0]), p, q);
    dY = arma::reshape(x.dtransparameters(indices_in[1]), q, q);
    dXYt = X * dY.t() + dX * Y.t();
    x.dtransparameters(indices_out[0]) = arma::vectorise(dXYt);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_out[0]), p, q);
    arma::mat dgrad_in_X = dgrad_out * Y + grad_out * dY;
    arma::mat dgrad_in_Y = dX.t() * grad_out + X.t() * dgrad_out;

    x.dgrad.elem(indices_in[0]) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_in[1]) += arma::vectorise(dgrad_in_Y);

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

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
