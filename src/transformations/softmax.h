/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// Softmax transformation:

class softmax:public transformations {

public:

  bool constraints;
  arma::vec theta, probs, Jdx;

  void transform(arguments_optim& x) {

    theta = x.transparameters(indices_in);
    probs = soft(theta, 1.00);
    x.transparameters.elem(indices_out) = probs;

  }

  void update_grad(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();
    x.grad(indices_in) += arma::vectorise(jacob.t() * x.grad(indices_out));

  }

  void dtransform(arguments_optim& x) {

    arma::vec dtrans = x.dtransparameters(indices_in);
    jacob = arma::diagmat(probs) - probs * probs.t();
    Jdx = jacob * dtrans;
    x.dtransparameters(indices_out) = Jdx;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec g_out  = x.grad(indices_out);
    arma::vec dg_out = x.dgrad(indices_out);

    arma::vec dp = Jdx;

    double gp  = arma::dot(g_out, probs);
    double gdp = arma::dot(g_out, dp);

    arma::vec dJg = g_out % dp - gp * dp - gdp * probs;

    arma::vec dg_theta = dJg + jacob * dg_out;

    x.dgrad(indices_in) += dg_theta;

  }

  void jacobian(arguments_optim& x) {

    jacob = arma::diagmat(probs) - probs * probs.t();

  }

  void dconstraints(arguments_optim& x) {

    arma::sp_mat first_derivatives(
      x.transparameters.n_elem, 1L
    );

    for(arma::uword i=0L; i < indices_out.n_elem; ++i) {
      first_derivatives(indices_out[i], 0L) = 1.00;
    }

    std::vector<arma::sp_mat> second_derivatives(
      1L,
      arma::sp_mat(x.transparameters.n_elem,
                   x.transparameters.n_elem)
    );

    append_constraint_derivatives(
      x,
      first_derivatives,
      second_derivatives
    );

  }

  void outcomes(arguments_optim& x) {

    arma::vec dconst(indices_out.n_elem); dconst.ones();

    vectors.resize(1);
    vectors[0] = dconst;
    names_vectors.resize(1);
    names_vectors[0] = "constraints_deriv";

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

softmax* choose_softmax(const Rcpp::List& trans_setup) {

  softmax* mytrans = new softmax();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];

  return mytrans;

}
