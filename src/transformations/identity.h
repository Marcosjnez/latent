/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2025
 */

// Identity (no transformation):

class identity:public transformations {

public:

  void transform(arguments_optim& x) {

    x.transparameters(indices_out[0]) = arma::vectorise(x.transparameters(indices_in[0]));

  }

  void update_grad(arguments_optim& x) {

    x.grad(indices_in[0]) += arma::vectorise(x.grad(indices_out[0]));

  }

  void dtransform(arguments_optim& x) {

    x.dtransparameters(indices_out[0]) = arma::vectorise(x.dtransparameters(indices_in[0]));

  }

  void update_dgrad(arguments_optim& x) {

    x.dgrad.elem(indices_in[0]) += arma::vectorise(x.dgrad(indices_out[0]));

  }

  void jacobian(arguments_optim& x) {

    arma::uvec v = indices_in[0];
    int p = v.n_elem;
    arma::mat I(p, p, arma::fill::eye);
    jacob = I;

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);

  }

};

identity* choose_identity(const Rcpp::List& trans_setup) {

  identity* mytrans = new identity();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
