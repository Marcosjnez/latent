/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2025
 */

// Crossproduct transformation:

class crossprod:public transformations {

public:

  int p, q;
  arma::mat X, gradmat;

  void transform(arguments_optim& x) {

    // std::memcpy(X.memptr(), x.transparameters(indices_in[0]).memptr(),
    //             indices_in[0].n_elem * sizeof(double));
    // arma::mat XtX = arma::symmatu(X.t() * X);
    // x.transparameters(indices_out[0]) = arma::vectorise(XtX);

  }

  void update_grad(arguments_optim& x) {

    // gradmat.set_size(p, q);
    // std::memcpy(gradmat.memptr(), grad_in.memptr(), grad_in.n_elem * sizeof(double));
    // grad_out = arma::vectorise(2*X * gradmat); grad_out.zeros();
    // x.grad(indices_in[0]) += arma::vectorise(2*X * gradmat);

  }

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    // jacob = arma::diagmat(arma::vectorise(2*X));
    // sum_djacob = 2*gradmat;

    // jacob = arma::join_cols(
    //   arma::kron(arma::eye(p, p), X.t()),
    //   arma::kron(arma::eye(p, p), X.t()) * arma::kron(arma::eye(p, p), arma::eye(p, p))
    // );


  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void M(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

    arma::vec chisq_p(p, arma::fill::value(p));

    vectors.resize(2);
    vectors[0] = dconstr;
    vectors[1] = chisq_p;

    matrices.resize(2);
    matrices[0] = jacob;
    matrices[1] = sum_djacob;

  }

};

crossprod* choose_crossprod(const Rcpp::List& trans_setup) {

  crossprod* mytrans = new crossprod();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::mat X(p, q);

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->X = X;

  return mytrans;

}
