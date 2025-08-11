/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Softmax transformation:

class softmax:public transformations {

public:

  void transform() {

    transparameters = soft(transparameters, 1.00);

  }

  void jacobian() {

    jacob = -transparameters * transparameters.t();
    jacob.diag() = transparameters % (1-transparameters);
    // jacob = jacob.cols(indices_in);
    grad = jacob.t() * grad;

    // for (size_t i = 0; i < grad.n_elem; ++i) {
    //   Rprintf("grad[%zu] = %f\n", i, grad[i]);
    // }

  }

  void d2jacobian() {



  }

  void dconstraints() {

    arma::vec v = arma::vec(transparameters.n_elem).fill(1.00);
    dconstr = v;

  }

  void outcomes() {

    dconstraints();
    vectors.resize(1);
    vectors[0] = dconstr;
    // vectors[1] = arma::conv_to<arma::vec>::from(indices_out);

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

softmax* choose_softmax(const Rcpp::List& trans_setup) {

  softmax* mytrans = new softmax();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
