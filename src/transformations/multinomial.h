/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2025
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

  void update_grad() {

    // jacob.set_size(transparameters.n_elem, parameters.n_elem);
    // jacob.zeros();

    dloglik = grad_in / small_peta(peta_indices);
    grad_out.resize(indices_in[0].n_elem); grad_out.zeros();
    grad_out.elem(peta_indices) += dloglik;

  }

  void update_hess() {

    const arma::uword n_out = transparameters.n_elem;
    const arma::uword n_in  = indices_in[0].n_elem;

    jacob.set_size(n_out, n_in);
    jacob.zeros();
    sum_djacob.set_size(n_in, n_in);
    sum_djacob.zeros();

    // jacob.cols(peta_indices) += arma::diagmat(1/small_peta(peta_indices));

    // Vectorised version:
    arma::uvec rows = arma::regspace<arma::uvec>(0, n_out - 1);
    arma::uvec lin_peta = rows + peta_indices * n_out;
    arma::vec dpeta = 1/small_peta(peta_indices);
    jacob.elem(lin_peta) += dpeta;

    // arma::vec djacob(n_in);
    // djacob.elem(peta_indices) += (1/small_peta(peta_indices))
    //   % (1/small_peta(peta_indices)) % grad_in;
    // sum_djacob = arma::diagmat(-djacob);

    // Vectorised version:
    lin_peta = peta_indices + peta_indices * n_in;
    sum_djacob.elem(lin_peta) -= dpeta % dpeta % grad_in;

    // hess_out = jacob.t() * hess_in * jacob + sum_djacob;

  }

  void dconstraints() {

    constraints = false;

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
