/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

// Logarithmic transformation:

class logarithm:public transformations {

public:

  arma::uvec indices_in, indices_out;
  arma::vec trans, dtrans;
  arma::mat jacob;

  void transform(arguments_optim& x) {

    trans = arma::trunc_log(x.transparameters(indices_in));
    x.transparameters.elem(indices_out) = trans;

  }

  void update_grad(arguments_optim& x) {

    x.grad(indices_in) += x.grad(indices_out) / x.transparameters(indices_in);

  }

  void dtransform(arguments_optim& x) {

    dtrans = x.dtransparameters(indices_in) / x.transparameters(indices_in);
    x.dtransparameters(indices_out) = dtrans;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec x_in   = x.transparameters(indices_in);
    arma::vec g_out  = x.grad(indices_out);
    arma::vec dg_out = x.dgrad(indices_out);

    arma::vec dg_in_add = (dg_out - g_out % dtrans) / x_in;

    x.dgrad(indices_in) += dg_in_add;

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::diagmat(1/x.transparameters(indices_in));

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

logarithm* choose_logarithm(const Rcpp::List& trans_setup) {

  logarithm* mytrans = new logarithm();

  // arma::uvec indices_in = trans_setup["indices_in"];
  // arma::uvec indices_out = trans_setup["indices_out"];
  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];

  return mytrans;

}
