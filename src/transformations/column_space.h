/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 08/09/2025
 */

// Column space transformation:

class column_space:public transformations {

public:

  arma::mat X;
  arma::vec coefficients, linear_preds;

  void transform(arguments_optim& x) {

    coefficients = x.transparameters(indices_in[0]);
    linear_preds = X * coefficients;
    x.transparameters.elem(indices_out[0]) = linear_preds;

  }

  void update_grad(arguments_optim& x) {

    jacob = X;
    // grad_out = jacob.t() * grad_in; grad_out.zeros();
    // Fill the gradient:
    x.grad(indices_in[0]) += arma::vectorise(jacob.t() * x.grad(indices_out[0]));

  }

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    sum_djacob = arma::diagvec(X);

    // hess_in = jacob.t() * hess_out * jacob + sum_djacob;

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = true;
    dconstr.set_size(indices_out[0].n_elem);
    dconstr.ones();

  }

  void M(arguments_optim& x) {

    // // New number of subjects in each class:
    // arma::vec freqs_i = arma::sum(x.freqs, 0).t();
    //
    // // New proportion of subjects in each class:
    // x.transparameters(indices_out[0]) = freqs_i / arma::accu(freqs_i);
    // x.transparameters(indices_in[0]) = arma::trunc_log(x.transparameters(indices_out[0]));


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

column_space* choose_column_space(const Rcpp::List& trans_setup) {

  column_space* mytrans = new column_space();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::mat X = trans_setup["X"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->X = X;

  return mytrans;

}
