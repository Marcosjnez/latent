/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 29/04/2025
 */

// Softmax transformation:

class softmaxt:public transformations {

public:

  arma::vec X;
  arma::uvec vector_indices;

  void transform() {

    X(vector_indices) = parameters;
    transparameters = softmax(X, 1.00);

  }

  void jacobian() {

    jacob = -transparameters * transparameters.t();
    jacob.diag() = transparameters % (1-transparameters);
    g = jacob * grad;
    g = g(vector_indices);
    g.elem( arma::find_nonfinite(g) ).zeros();

  }

  void d2jacobian() {



  }

  void dconstraints() {

    arma::vec v = arma::vec(transparameters.n_elem).fill(1.00);
    dconstr = v;

  }

  void outcomes() {

    dconstraints();
    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = arma::conv_to<arma::vec>::from(target_indices);

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

softmaxt* choose_softmaxt(Rcpp::List trans_setup) {

  softmaxt* mytrans = new softmaxt();

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
