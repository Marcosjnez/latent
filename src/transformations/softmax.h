/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

// Softmax transformation:

class softmax:public transformations {

public:

  arma::vec theta, probs, Jdx;

  void transform(arguments_optim& x) {

    theta = x.transparameters(indices_in[0]);
    probs = soft(theta, 1.00);
    x.transparameters.elem(indices_out[0]) = probs;

  }

  void update_grad(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();
    x.grad(indices_in[0]) += arma::vectorise(jacob.t() * x.grad(indices_out[0]));

  }

  void dtransform(arguments_optim& x) {

    arma::vec dtrans = x.dtransparameters(indices_in[0]);
    jacob = arma::diagmat(probs) - probs * probs.t();
    Jdx = jacob * dtrans;
    x.dtransparameters(indices_out[0]) = Jdx;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec g_out  = x.grad(indices_out[0]);
    arma::vec dg_out = x.dgrad(indices_out[0]);

    arma::vec dp = Jdx;

    double gp  = arma::dot(g_out, probs);
    double gdp = arma::dot(g_out, dp);

    arma::vec dJg = g_out % dp - gp * dp - gdp * probs;

    arma::vec dg_theta = dJg + jacob * dg_out;

    x.dgrad(indices_in[0]) += dg_theta;

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = true;
    dconstr.set_size(indices_out[0].n_elem);
    dconstr.ones();

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

softmax* choose_softmax(const Rcpp::List& trans_setup) {

  softmax* mytrans = new softmax();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;

  return mytrans;

}
