/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 25/05/2026
 */

// Square root transformation:

class sqrt_vector: public transformations {

public:

  int n_in, n_out;
  arma::uvec indices_in, indices_out;
  arma::vec input, output, dinput, doutput, grad_out, dgrad_out;
  arma::mat jacob;

  void transform(arguments_optim& x) {

    input = x.transparameters.elem(indices_in);
    output = arma::sqrt(input);

    x.transparameters.elem(indices_out) = output;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad.elem(indices_out);

    x.grad.elem(indices_in) += grad_out / (2.0 * output);

  }

  void dtransform(arguments_optim& x) {

    dinput = x.dtransparameters.elem(indices_in);
    doutput = dinput / (2.0 * output);

    x.dtransparameters.elem(indices_out) = doutput;

  }

  void update_dgrad(arguments_optim& x) {

    dgrad_out = x.dgrad.elem(indices_out);

    arma::vec a = 1.0 / (2.0 * output);
    arma::vec da = -dinput / (4.0 * output % output % output);

    x.dgrad.elem(indices_in) += dgrad_out % a + grad_out % da;

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    for (int i = 0; i < n_out; ++i) {
      jacob(i, i) = 1.0 / (2.0 * output(i));
    }

  }

  void update_vcov(arguments_optim& x) {

    x.vcov(indices_out, indices_out) =
      jacob * x.vcov(indices_in, indices_in) * jacob.t();

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

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

sqrt_vector* choose_sqrt_vector(const Rcpp::List& trans_setup) {

  sqrt_vector* mytrans = new sqrt_vector();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  arma::uvec indices_input = indices_in[0];
  arma::uvec indices_output = indices_out[0];

  int n_in = indices_input.n_elem;
  int n_out = indices_output.n_elem;

  mytrans->indices_in = indices_input;
  mytrans->indices_out = indices_output;
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;

  return mytrans;

}
