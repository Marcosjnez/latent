/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/09/2026
 */

// Matrix-vector transformation:
// X is p x q, y has length q, and the output has length p.

class matrix_vector:public transformations {

public:

  int p, q;
  arma::uvec indices_X, indices_y, indices_Xy;
  arma::mat X, dX, grad_in_X;
  arma::vec y, dy, dXy, grad_out, grad_in_y;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    y = x.transparameters(indices_y);
    x.transparameters(indices_Xy) = X*y;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad(indices_Xy);
    grad_in_X = grad_out*y.t();
    grad_in_y = X.t()*grad_out;

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_y) += grad_in_y;

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dy = x.dtransparameters(indices_y);
    dXy = dX*y+X*dy;
    x.dtransparameters(indices_Xy) = dXy;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec dgrad_out = x.dgrad(indices_Xy);
    arma::mat dgrad_in_X = dgrad_out*y.t()+grad_out*dy.t();
    arma::vec dgrad_in_y = dX.t()*grad_out+X.t()*dgrad_out;

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_y) += dgrad_in_y;

  }

  void jacobian(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    y = x.transparameters(indices_y);

    arma::mat I_p = arma::eye(p, p);
    arma::mat Jx = arma::kron(y.t(), I_p);
    arma::mat Jy = X;

    jacob = arma::join_rows(Jx, Jy);

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

matrix_vector* choose_matrix_vector(const Rcpp::List& trans_setup) {

  matrix_vector* mytrans = new matrix_vector();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  if(p < 1 || q < 1) {
    Rf_error("matrix_vector requires positive p and q dimensions.");
  }

  if(indices_in.size() != 2L || indices_out.size() != 1L) {
    Rf_error("matrix_vector requires X and y inputs and one output.");
  }

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_y = indices_in[1];
  arma::uvec indices_Xy = indices_out[0];

  if(indices_X.n_elem != static_cast<arma::uword>(p*q) ||
     indices_y.n_elem != static_cast<arma::uword>(q) ||
     indices_Xy.n_elem != static_cast<arma::uword>(p)) {
    Rf_error("The matrix_vector parameter indices have incompatible dimensions.");
  }

  mytrans->indices_in = arma::join_cols(indices_X, indices_y);
  mytrans->indices_out = indices_Xy;
  mytrans->indices_X = indices_X;
  mytrans->indices_y = indices_y;
  mytrans->indices_Xy = indices_Xy;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
