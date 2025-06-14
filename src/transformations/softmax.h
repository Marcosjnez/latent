/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/06/2025
 */

// Softmax transformation:

class softmax:public transformations {

public:

  arma::vec X;
  arma::uvec vector_indices;

  void transform() {

    X(vector_indices) = parameters;
    transparameters = soft(X, 1.00);

  }

  void jacobian() {

    jacob = -transparameters * transparameters.t();
    jacob.diag() = transparameters % (1-transparameters);
    jacob = jacob.cols(vector_indices);
    g = jacob.t() * grad;

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

softmax* choose_softmax(const Rcpp::List& trans_setup) {

  softmax* mytrans = new softmax();

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
