/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/09/2026
 */

// A = inv(I-X) transformation:

class inverse_of_I_minusX:public transformations {

public:

  int p;
  arma::mat X, I, A, dX, dA, grad_out, grad_in;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in), p, p);

    if(!X.is_finite()) {
      Rf_error("inverse_of_I_minusX requires a finite input matrix.");
    }

    bool ok = arma::inv(A, I-X);

    if(!ok || !A.is_finite()) {
      Rf_error("inverse_of_I_minusX: I-X is singular or its inverse is nonfinite.");
    }

    x.transparameters(indices_out) = arma::vectorise(A);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out), p, p);
    grad_in = A.t()*grad_out*A.t();
    x.grad(indices_in) += arma::vectorise(grad_in);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_in), p, p);
    dA = A*dX*A;
    x.dtransparameters(indices_out) = arma::vectorise(dA);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_out), p, p);

    arma::mat At = A.t();
    arma::mat dAt = dA.t();
    arma::mat dgrad_in = dAt*grad_out*At+
                         At*dgrad_out*At+
                         At*grad_out*dAt;

    x.dgrad(indices_in) += arma::vectorise(dgrad_in);

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::kron(A.t(), A);

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

inverse_of_I_minusX* choose_inverse_of_I_minusX(const Rcpp::List& trans_setup) {

  inverse_of_I_minusX* mytrans = new inverse_of_I_minusX();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];

  if(p < 1) {
    Rf_error("inverse_of_I_minusX requires a positive p dimension.");
  }

  if(indices_in.size() != 1L || indices_out.size() != 1L) {
    Rf_error("inverse_of_I_minusX requires one input and one output.");
  }

  const arma::uword n2 = static_cast<arma::uword>(p)*p;

  if(indices_in[0].n_elem != n2 || indices_out[0].n_elem != n2) {
    Rf_error("The inverse_of_I_minusX parameter indices have incompatible dimensions.");
  }

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];
  mytrans->p = p;
  mytrans->I = arma::mat(p, p, arma::fill::eye);

  return mytrans;

}
