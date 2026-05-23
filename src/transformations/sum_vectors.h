/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 16/05/2026
 */

class sum_vectors: public transformations {

public:

  int V, L, n_in, n_out;
  std::vector<arma::uvec> indices_vectors;
  arma::uvec indices_in, indices_out;
  arma::vec output, doutput, grad_out, dgrad_out;
  arma::mat jacob;

  void transform(arguments_optim& x) {

    output.zeros();

    for (int v = 0; v < V; ++v) {
      output += x.transparameters.elem(indices_vectors[v]);
    }

    x.transparameters.elem(indices_out) = output;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad.elem(indices_out);

    for (int v = 0; v < V; ++v) {
      x.grad.elem(indices_vectors[v]) += grad_out;
    }

  }

  void dtransform(arguments_optim& x) {

    doutput.zeros();

    for (int v = 0; v < V; ++v) {
      doutput += x.dtransparameters.elem(indices_vectors[v]);
    }

    x.dtransparameters.elem(indices_out) = doutput;

  }

  void update_dgrad(arguments_optim& x) {

    dgrad_out = x.dgrad.elem(indices_out);

    for (int v = 0; v < V; ++v) {
      x.dgrad.elem(indices_vectors[v]) += dgrad_out;
    }

  }

  void jacobian(arguments_optim& x) {

    jacob.set_size(n_out, n_in);
    jacob.zeros();

    int offset = 0;
    for (int v = 0; v < V; ++v) {
      for (int l = 0; l < L; ++l) {
        jacob(l, offset + l) = 1.0;
      }
      offset += L;
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

    vectors.resize(1);

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

sum_vectors* choose_sum_vectors(const Rcpp::List& trans_setup) {

  sum_vectors* mytrans = new sum_vectors();

  std::vector<arma::uvec> indices_in_list = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out_list = trans_setup["indices_out"];

  int V = indices_in_list.size();
  int L = indices_in_list[0].n_elem;
  int n_out = indices_out_list[0].n_elem;
  int n_in = 0;

  arma::uvec indices_in = indices_in_list[0];
  n_in += indices_in_list[0].n_elem;

  for (int v = 1; v < V; ++v) {
    indices_in = arma::join_cols(indices_in, indices_in_list[v]);
    n_in += indices_in_list[v].n_elem;
  }

  arma::vec output(L, arma::fill::zeros);
  arma::vec doutput(L, arma::fill::zeros);

  mytrans->indices_vectors = indices_in_list;
  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out_list[0];
  mytrans->V = V;
  mytrans->L = L;
  mytrans->n_in = n_in;
  mytrans->n_out = n_out;
  mytrans->output = output;
  mytrans->doutput = doutput;

  return mytrans;

}
