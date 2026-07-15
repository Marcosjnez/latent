/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/03/2026
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

    // x.vcov(indices_out, indices_out) =
    //   jacob * x.vcov(indices_in, indices_in) * jacob.t();

    arma::mat vcov_in(indices_in.n_elem, indices_in.n_elem);

    for(arma::uword j = 0L; j < indices_in.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_in.n_elem; ++i) {
        vcov_in(i, j) = x.vcov(indices_in[i], indices_in[j]);
      }
    }

    arma::mat vcov_out = jacob * vcov_in * jacob.t();

    for(arma::uword j = 0L; j < indices_out.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {
        x.vcov(indices_out[i], indices_out[j]) = vcov_out(i, j);
      }
    }

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

exponential* choose_exponential(const Rcpp::List& trans_setup) {

  exponential* mytrans = new exponential();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];

  return mytrans;

}
