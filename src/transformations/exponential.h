/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/06/2025
 */

// Exponential transformation:

class exponential:public transformations {

public:

  arma::vec X;
  arma::uvec vector_indices;

  void transform() {

    X(vector_indices) = parameters;
    transparameters = 0.05 + arma::trunc_exp(X);

  }

  void jacobian() {

    // jacob = arma::diagmat(transparameters);
    // jacob = jacob.cols(vector_indices);
    // g = jacob * grad;
    g = grad % transparameters;
    g = g(vector_indices);

  }

  void d2jacobian() {

  }

  void dconstraints() {

  }

  void outcomes() {

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

exponential* choose_exponential(const Rcpp::List& trans_setup) {

  exponential* mytrans = new exponential();

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
