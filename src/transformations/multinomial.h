/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/08/2025
 */

// Logarithm multinomial probability transformation:

class multinomial:public transformations {

public:

  arma::vec small_peta;
  arma::vec small_logpeta;
  arma::vec logpeta;
  arma::uvec peta_indices;

  void transform() {

    small_peta = transparameters; // JxIxK
    peta_indices = indices_in[1]; // SxJxIxK indices
    small_logpeta = arma::trunc_log(small_peta); // JxIxK
    transparameters = small_logpeta.elem(peta_indices); // SxJxI
    // logpeta = small_logpeta.elem(peta_indices); // SxJxIxK
    // transparameters = logpeta.elem(dummy_indices); // SxJxI

  }

  void jacobian() {

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();

    arma::vec dloglik = grad;
    // arma::uvec final_indices = peta_indices(dummy_indices);
    grad.resize(indices_in[0].n_elem); grad.zeros();
    // grad.elem(final_indices) += dloglik;
    grad.elem(peta_indices) += dloglik;

  }

  void d2jacobian() {

  }

  void dconstraints() {

  }

  void outcomes() {

    matrices.resize(1);
    matrices[0] = jacob;

    cubes.resize(1);
    cubes[0] = jacob2;

  }

};

multinomial* choose_multinomial(const Rcpp::List& trans_setup) {

  multinomial* mytrans = new multinomial();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
