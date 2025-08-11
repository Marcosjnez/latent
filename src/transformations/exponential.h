/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Exponential transformation:

class exponential:public transformations {

public:

  void transform() {

    transparameters = arma::trunc_exp(transparameters);

  }

  void jacobian() {

    // jacob = arma::diagmat(transparameters);
    // jacob = jacob.cols(vector_indices);
    // g = jacob * grad;
    grad = grad % transparameters;

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

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
