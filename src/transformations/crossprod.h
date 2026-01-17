/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/10/2025
 */

// Crossproduct transformation:

class crossprod:public transformations {

public:

  int p;
  arma::mat X, grad_out, dX, dXtX, Dp;
  arma::uvec lower_diag;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in[0]), p, p);
    arma::mat XtX = X.t() * X;
    x.transparameters(indices_out[0]) = arma::vectorise(XtX.elem(lower_diag));

  }

  void update_grad(arguments_optim& x) {

    grad_out.elem(lower_diag) = x.grad(indices_out[0]);
    grad_out = arma::symmatl(grad_out);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal
    x.grad(indices_in[0]) += arma::vectorise(2*X * grad_out);

    // arma::vec v = x.grad(indices_in[0]);
    // for (arma::uword i = 0; i < v.n_elem; ++i) {
    //   Rprintf("%.6f%s", v(i), (i + 1 < v.n_elem) ? " " : "\n"); // space-separated, then newline
    // }

    // arma::mat I(p, p, arma::fill::eye);
    // jacob = 2*Dp.t() * arma::kron(I, X.t());
    // x.grad(indices_in[0]) += jacob.t() * x.grad(indices_out[0]);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_in[0]), p, p);
    dXtX = X.t() * dX + dX.t() * X;
    x.dtransparameters(indices_out[0]) = arma::vectorise(dXtX.elem(lower_diag));

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out(p, p, arma::fill::zeros);
    dgrad_out.elem(lower_diag) = x.dgrad.elem(indices_out[0]);
    dgrad_out = arma::symmatl(dgrad_out);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    x.dgrad.elem(indices_in[0]) += arma::vectorise(2*dX * grad_out +
                                                   2*X * dgrad_out);

  }

  void jacobian(arguments_optim& x) {

    arma::mat I(p, p, arma::fill::eye);
    jacob = 2*Dp.t() * arma::kron(I, X.t());

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    arma::vec chisq_p(p, arma::fill::value(p));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

  }

};

crossprod* choose_crossprod(const Rcpp::List& trans_setup) {

  crossprod* mytrans = new crossprod();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];

  arma::mat X(p, p);
  arma::mat grad_out(p, p, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));
  arma::mat Dp = duplication(p);

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->X = X;
  mytrans->grad_out = grad_out;
  mytrans->lower_diag = lower_diag;
  mytrans->Dp = Dp;

  return mytrans;

}
