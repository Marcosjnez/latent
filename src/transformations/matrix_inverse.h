/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

// Matrix inverse transformation:

class matrix_inverse:public transformations {

public:

  arma::uvec indices_in, indices_out;
  int p;
  arma::mat X, Xinv, dX, dXinv, grad_out, grad_in, jacob;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in), p, p);
    Xinv = arma::inv(X);
    x.transparameters(indices_out) = arma::vectorise(Xinv);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out), p, p);
    grad_in = -Xinv.t() * grad_out * Xinv.t();
    x.grad(indices_in) += arma::vectorise(grad_in);

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_in), p, p);
    dXinv = -Xinv * dX * Xinv;
    x.dtransparameters(indices_out) = arma::vectorise(dXinv);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad(indices_out), p, p);

    arma::mat A  = Xinv.t();
    arma::mat dA = dXinv.t();
    arma::mat dgrad_in = -(dA * grad_out * A +
                           A  * dgrad_out * A +
                           A  * grad_out * dA);

    x.dgrad(indices_in) += arma::vectorise(dgrad_in);

  }

  void jacobian(arguments_optim& x) {

    jacob = -arma::kron(Xinv.t(), Xinv);

  }

  void update_vcov(arguments_optim& x) {

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

matrix_inverse* choose_matrix_inverse(const Rcpp::List& trans_setup) {

  matrix_inverse* mytrans = new matrix_inverse();

  arma::uvec indices_in = trans_setup["indices_in"];
  arma::uvec indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;

  return mytrans;

}
