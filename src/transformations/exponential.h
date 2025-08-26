/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2025
 */

// Exponential transformation:

class exponential:public transformations {

public:

  void transform() {

    transparameters = arma::trunc_exp(transparameters);

  }

  void update_grad() {

    // jacob = arma::diagmat(transparameters);
    // g = jacob.t() * grad;
    grad_out = grad_in % transparameters;

  }

  void update_hess() {

    jacob = arma::diagmat(transparameters);
    sum_djacob = arma::diagmat(transparameters % grad_in);

    // hess_out = jacob.t() * hess_in * jacob + sum_djacob;

  }

  void dconstraints() {

    constraints = false;

  }

  void outcomes() {

    int p = transparameters.n_elem;
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
