/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/08/2025
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

  void jacobian() {

    jacob = arma::diagmat(transparameters) - transparameters * transparameters.t();
    grad_out = jacob.t() * grad_in;

  }

  void d2jacobian() {

    jacob = arma::diagmat(transparameters) - transparameters * transparameters.t();
    sum_djacob = softmax_chain_second_term(transparameters, grad_in, jacob);

    // hess_out = jacob.t() * hess_in * jacob + sum_djacob;

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

arma::cube softmax_hessian(const arma::vec& theta) {
  arma::vec p = arma::normalise(arma::exp(theta), 1); // softmax
  const arma::uword n = p.n_elem;

  arma::mat J = arma::diagmat(p) - p * p.t(); // Jacobian
  arma::cube H(n, n, n, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec a = -p;
    a(i) += 1.0; // a = e_i - p
    H.slice(i) = p(i) * (a * a.t() - J);
  }
  return H;
}

