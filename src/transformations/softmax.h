/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
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

  arma::vec theta, probs;

  void transform(arguments_optim& x) {

    theta = x.transparameters(indices_in[0]);
    probs = soft(theta, 1.00);
    x.transparameters.elem(indices_out[0]) = probs;

  }

  void update_grad(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();
    // Fill the gradient:
    x.grad(indices_in[0]) += arma::vectorise(jacob.t() * x.grad(indices_out[0]));

  }

  void update_dparam(arguments_optim& x) {

  }

  void update_dgrad(arguments_optim& x) {

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();

  }

  void update_hess(arguments_optim& x) {

    sum_djacob = softmax_chain_second_term(probs, x.grad(indices_out[0]), jacob);

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = true;
    dconstr.set_size(indices_out[0].n_elem);
    dconstr.ones();

  }

  void outcomes(arguments_optim& x) {

    dconstraints(x);
    int p = indices_out[0].n_elem;
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
