/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2025
 */

// Identity (no transformation):

class identity:public transformations {

public:

  // arma::vec X;
  // arma::uvec vector_indices;

  void transform(arguments_optim& x) {

    // Rf_error("id 18");
    // X = parameters;
    // transparameters = X;
    // Rf_error("id 21");

  }

  void update_grad(arguments_optim& x) {

    // Rprintf("vector_indices:\n");
    // for (arma::uword i = 0; i < vector_indices.n_elem; ++i) {
    //   Rprintf("%u ", vector_indices[i]);
    // }
    // Rprintf("\n\n");
    //
    // Rprintf("x.grad:\n");
    // for (arma::uword j = 0; j < grad.n_elem; ++j) {
    //   Rprintf("%g ", grad[j]);
    // }
    // Rprintf("\n\n");
    //
    // Rf_error("30");

    // jacob.set_size(parameters.n_elem, parameters.n_elem);
    // jacob.eye();
    // jacob = jacob.cols(vector_indices);
    // g = jacob.t() * grad;
    // grad = grad;

  }

  void dparam(arguments_optim& x) {

  }

  void update_dgrad(arguments_optim& x) {

  }

  void jacobian(arguments_optim& x) {

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

  }

};

identity* choose_identity(const Rcpp::List& trans_setup) {

  identity* mytrans = new identity();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
