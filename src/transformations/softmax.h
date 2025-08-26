/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2025
 */

// Softmax transformation:

arma::mat softmax_chain_second_term(const arma::vec& p,
                                    const arma::vec& grad_in,
                                    const arma::mat& J) {
  arma::vec w = grad_in % p;                                    // elementwise
  double sumw = arma::accu(w);

  arma::mat A = arma::diagmat(w) - w * p.t() - p * w.t() + sumw * (p * p.t());

  return A - sumw * J;  // equals Σ_i grad_in(i) * H_i
}

class softmax:public transformations {

public:

  arma::vec theta;

  void transform() {

    theta = transparameters;
    transparameters = soft(transparameters, 1.00);

  }

  void update_grad() {

    jacob = arma::diagmat(transparameters) - transparameters * transparameters.t();
    grad_out = jacob.t() * grad_in;

  }

  void update_hess() {

    sum_djacob = softmax_chain_second_term(transparameters, grad_in, jacob);

    // hess_out = jacob.t() * hess_in * jacob + sum_djacob;

  }

  void dconstraints() {

    constraints = true;
    dconstr.set_size(indices_out[0].n_elem);
    dconstr.ones();

  }

  void outcomes() {

    dconstraints();
    int p = transparameters.n_elem;
    arma::vec chisq_p(p, arma::fill::value(p));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

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
