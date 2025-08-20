/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/08/2025
 */

// Logarithm multinomial probability transformation:

class multinomial:public transformations {

public:

  arma::vec small_peta;
  arma::vec small_logpeta;
  arma::vec logpeta;
  arma::uvec peta_indices;
  arma::vec dloglik;

  void transform() {

    small_peta = transparameters; // JxIxK
    peta_indices = indices_in[1]; // SxJxIxK indices
    small_logpeta = arma::trunc_log(small_peta); // JxIxK
    transparameters = small_logpeta.elem(peta_indices); // SxJxI

  }

  void jacobian() {

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();

    dloglik = grad_in / small_peta(peta_indices);
    grad_out.resize(indices_in[0].n_elem); grad_out.zeros();
    grad_out.elem(peta_indices) += dloglik;

  }

  void d2jacobian() {

    jacob.set_size(transparameters.n_elem, indices_in[0].n_elem);
    jacob.zeros();
    jacob.cols(peta_indices) += arma::diagmat(1/small_peta(peta_indices));

    // sum_djacob.set_size(indices_in[0].n_elem, indices_in[0].n_elem);
    // sum_djacob.zeros();
    arma::vec djacob(indices_in[0].n_elem);
    djacob.elem(peta_indices) += (1/small_peta(peta_indices))
      % (1/small_peta(peta_indices)) % grad_in;
    sum_djacob = arma::diagmat(-djacob);

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
