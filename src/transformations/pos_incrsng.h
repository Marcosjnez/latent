/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 23/06/2026
 */

// Positive increasing transformation (lower_triangular_matrix_of_1s * exp(x)):

class pos_incrsng: public transformations {

public:

  int n_in, n_out;
  arma::uvec indices_in, indices_out;
  arma::vec pos_values, pos_incrsng_values;
  arma::vec dpos_values, dpos_incrsng_values;
  arma::vec grad_out, dgrad_out;
  arma::mat lower_trng, jacob;

  void transform(arguments_optim& x) {

    pos_values = arma::trunc_exp(x.transparameters.elem(indices_in));
    pos_incrsng_values = lower_trng * pos_values;

    x.transparameters.elem(indices_out) = pos_incrsng_values;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad.elem(indices_out);

    x.grad.elem(indices_in) +=
      pos_values % (lower_trng.t() * grad_out);

  }

  void dtransform(arguments_optim& x) {

    dpos_values = pos_values % x.dtransparameters.elem(indices_in);
    dpos_incrsng_values = lower_trng * dpos_values;

    x.dtransparameters.elem(indices_out) = dpos_incrsng_values;

  }

  void update_dgrad(arguments_optim& x) {

    dgrad_out = x.dgrad.elem(indices_out);

    x.dgrad.elem(indices_in) +=
      dpos_values % (lower_trng.t() * grad_out) +
      pos_values % (lower_trng.t() * dgrad_out);

  }

  void jacobian(arguments_optim& x) {

    jacob = lower_trng * arma::diagmat(pos_values);

  }

  void update_vcov(arguments_optim& x) {

    x.vcov(indices_out, indices_out) =
      jacob * x.vcov(indices_in, indices_in) * jacob.t();

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

pos_incrsng* choose_pos_incrsng(const Rcpp::List& trans_setup) {

  pos_incrsng* mytrans = new pos_incrsng();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  int n = indices_in[0].n_elem;
  arma::mat lower_trng = arma::trimatl(arma::ones<arma::mat>(n, n));

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];
  mytrans->lower_trng = lower_trng;
  mytrans->n_in = n;
  mytrans->n_out = indices_out[0].n_elem;

  return mytrans;

}

