/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2025
 */

// Exponential transformation:

class exponential:public transformations {

public:

  arma::vec trans;

  void transform(arguments_optim& x) {

    trans = arma::trunc_exp(x.transparameters(indices_in[0]));
    x.transparameters.elem(indices_out[0]) = trans;

  }

  void update_grad(arguments_optim& x) {

    // jacob = arma::diagmat(trans);
    // g = jacob.t() * grad;
    // Fill the gradient:
    x.grad(indices_in[0]) += x.grad(indices_out[0]) % trans;

  }

  void update_dparam(arguments_optim& x) {

  }

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    jacob = arma::diagmat(trans);
    sum_djacob = arma::diagmat(trans % x.grad(indices_out[0]));

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov(arguments_optim& x) {

    jacob = arma::diagmat(trans);

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void M(arguments_optim& x) {

    // x.transparameters(indices_in[0]) = arma::trunc_log(x.transparameters(indices_out[0]));

  }

  void outcomes(arguments_optim& x) {

    int p = indices_out[0].n_elem;
    arma::vec chisq_p(p, arma::fill::value(1.00));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

  }

};

exponential* choose_exponential(const Rcpp::List& trans_setup) {

  exponential* mytrans = new exponential();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
