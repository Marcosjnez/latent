/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 29/04/2025
 */

// Identity (no transformation):

class identity:public transformations {

public:

  arma::vec X;
  arma::uvec vector_indices;

  void transform() {

    // Rf_error("id 18");
    X(vector_indices) = parameters;
    transparameters = X;
    // Rf_error("id 20");

  }

  void jacobian() {

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

    g = grad(vector_indices);

  }

  void d2jacobian() {

  }

  void dconstraints() {

  }

  void outcomes() {

  }

};

identity* choose_identity(Rcpp::List trans_setup) {

  identity* mytrans = new identity();

  arma::uvec indices = trans_setup["indices"];
  arma::uvec target_indices = trans_setup["target_indices"];
  arma::uvec vector_indices = trans_setup["vector_indices"];
  arma::vec X = trans_setup["X"];

  mytrans->indices = indices;
  mytrans->target_indices = target_indices;
  mytrans->vector_indices = vector_indices;
  mytrans->X = X;

  return mytrans;

}
