/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

// Exponential transformation:

class exponential:public transformations {

public:

  arma::uvec indices_in, indices_out;
  arma::vec trans, dtrans;
  arma::mat jacob;

  void transform(arguments_optim& x) {

    trans = arma::trunc_exp(x.transparameters(indices_in));
    x.transparameters.elem(indices_out) = trans;

  }

  void update_grad(arguments_optim& x) {

    // jacob = arma::diagmat(trans);
    // g = jacob.t() * grad;
    x.grad(indices_in) += x.grad(indices_out) % trans;

  }

  void dtransform(arguments_optim& x) {

    dtrans = x.dtransparameters(indices_in) % trans;
    x.dtransparameters(indices_out) = dtrans;

  }

  void update_dgrad(arguments_optim& x) {

    x.dgrad.elem(indices_in) += x.dgrad.elem(indices_out) % trans +
                                   x.grad(indices_out) % dtrans;

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::diagmat(trans);

  }

  void update_vcov(arguments_optim& x) {

    x.vcov(indices_out, indices_out) = jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    int p = indices_out.n_elem;
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

  arma::uvec indices_in = trans_setup["indices_in"];
  arma::uvec indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
