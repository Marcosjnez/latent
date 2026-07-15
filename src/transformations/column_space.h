/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/03/2026
 */

// Column space transformation:

class column_space:public transformations {

public:

  arma::uvec indices_in, indices_out;
  arma::mat X, coefs, dcoefs, linear_preds, jacob;

  void transform(arguments_optim& x) {

    arma::vec values = x.transparameters(indices_in);
    std::memcpy(coefs.memptr(), values.memptr(), sizeof(double)*values.n_elem);
    linear_preds = X * coefs;
    x.transparameters.elem(indices_out) = arma::vectorise(linear_preds);

  }

  void update_grad(arguments_optim& x) {

    arma::mat I(coefs.n_cols, coefs.n_cols, arma::fill::eye);
    jacob = arma::kron(I, X);
    x.grad(indices_in) += arma::vectorise(jacob.t() * x.grad(indices_out));

  }

  void dtransform(arguments_optim& x) {

    arma::vec dvalues = x.dtransparameters(indices_in);
    std::memcpy(dcoefs.memptr(), dvalues.memptr(), sizeof(double)*dvalues.n_elem);
    arma::mat dlinear_preds = X * dcoefs;
    x.dtransparameters(indices_out) = arma::vectorise(dlinear_preds);

  }

  void update_dgrad(arguments_optim& x) {

    x.dgrad(indices_in) += arma::vectorise(jacob.t() * x.dgrad(indices_out));

  }

  void jacobian(arguments_optim& x) {

    arma::mat I(coefs.n_cols, coefs.n_cols, arma::fill::eye);
    jacob = arma::kron(I, X);

  }

  void update_vcov(arguments_optim& x) {

    // x.vcov(indices_out, indices_out) =
    //   jacob * x.vcov(indices_in, indices_in) * jacob.t();

    arma::mat vcov_in(indices_in.n_elem, indices_in.n_elem);

    for(arma::uword j = 0L; j < indices_in.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_in.n_elem; ++i) {
        vcov_in(i, j) = x.vcov(indices_in[i], indices_in[j]);
      }
    }

    arma::mat vcov_out = jacob * vcov_in * jacob.t();

    for(arma::uword j = 0L; j < indices_out.n_elem; ++j) {
      for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {
        x.vcov(indices_out[i], indices_out[j]) = vcov_out(i, j);
      }
    }

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

column_space* choose_column_space(const Rcpp::List& trans_setup) {

  column_space* mytrans = new column_space();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  arma::mat X = trans_setup["X"];

  int q = indices_in[0].n_elem / X.n_cols;
  arma::mat coefs(X.n_cols, q);

  mytrans->indices_in = indices_in[0];
  mytrans->indices_out = indices_out[0];
  mytrans->X = X;
  mytrans->coefs = coefs;
  mytrans->dcoefs = coefs;

  return mytrans;

}
