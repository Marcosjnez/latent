/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 29/04/2025
 */

// Exponential transformation:

class expo:public transformations {

public:

  arma::vec X;
  arma::uvec vector_indices;

  void transform() {

    X(vector_indices) = parameters;
    transparameters = arma::trunc_exp(X);

  }

  void jacobian() {

    jacob = transparameters;
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

expo* choose_expo(Rcpp::List trans_setup) {

  expo* mytrans = new expo();

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
