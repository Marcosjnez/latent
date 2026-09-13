/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/09/2026
 */

// Square transformation:

class square_vector: public transformations {

public:

  int n_in, n_out;
  arma::vec input, output, dinput, doutput, grad_out, dgrad_out;

  void transform(arguments_optim& x) {

    input = x.transparameters.elem(indices_in);
    output = input % input;

    x.transparameters.elem(indices_out) = output;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad.elem(indices_out);

    x.grad.elem(indices_in) += 2.0 * input % grad_out;

  }

  void dtransform(arguments_optim& x) {

    dinput = x.dtransparameters.elem(indices_in);
    doutput = 2.0 * input % dinput;

    x.dtransparameters.elem(indices_out) = doutput;

  }

  void update_dgrad(arguments_optim& x) {

    dgrad_out = x.dgrad.elem(indices_out);

    x.dgrad.elem(indices_in) +=
      2.0 * input % dgrad_out +
      2.0 * dinput % grad_out;

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    for(int i = 0; i < n_out; ++i) {
      jacob(i, i) = 2.0 * input(i);
    }

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

square_vector* choose_square_vector(const Rcpp::List& trans_setup) {

  square_vector* mytrans = new square_vector();

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
