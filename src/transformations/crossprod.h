/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 25/08/2026
 */

// Crossproduct transformation:

class crossprod:public transformations {

public:

  int p;
  arma::uvec unique_indices_out, lower_diag;
  arma::mat X, grad_out, dX, dXtX;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in), p, p);
    arma::mat XtX = X.t() * X;
    x.transparameters(indices_out) = arma::vectorise(XtX);

  }

  void update_grad(arguments_optim& x) {

    grad_out = matrix_output_gradient(x.grad, p);

    x.grad(indices_in) += arma::vectorise(X*(grad_out+grad_out.t()));

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_in), p, p);
    dXtX = X.t() * dX + dX.t() * X;
    x.dtransparameters(indices_out) = arma::vectorise(dXtX);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = matrix_output_gradient(x.dgrad, p);

    x.dgrad.elem(indices_in) += arma::vectorise(
      dX*(grad_out+grad_out.t())+X*(dgrad_out+dgrad_out.t()));

  }

  void jacobian(arguments_optim& x) {

    jacob.zeros(p*p, p*p);
    arma::mat E(p, p, arma::fill::zeros);

    for(int k=0; k < p*p; ++k) {

      E.zeros();
      E[k] = 1.00;
      jacob.col(k) =
        arma::vectorise(E.t()*X+X.t()*E);

    }

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

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
  arma::uvec unique_indices_out = arma::unique(indices_out[0]);

  mytrans->indices_in = indices_in[0];
  mytrans->set_matrix_output(indices_out[0]);
  mytrans->p = p;
  mytrans->X = X;
  mytrans->grad_out = grad_out;
  mytrans->unique_indices_out = unique_indices_out;
  mytrans->lower_diag = lower_diag;

  return mytrans;

}
