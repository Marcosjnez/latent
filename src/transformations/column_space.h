/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/08/2026
 */

// Column space transformation:

class column_space:public transformations {

public:

  arma::mat X, coefs, dcoefs, linear_preds;

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

  void dconstraints(arguments_optim& x) {

    // When X is not full-rank, there are linear dependencies in the columns.
    // If r = rank(X), then there are n-r constraints in each row of X.
    // A basis for their derivatives is null(X.t()) because
    // null((X.t()))

    // Columns of N span the orthogonal complement of the column space of X:
    arma::mat N = arma::null(X.t());

    if(N.n_cols == 0L) return;

    arma::uword n = X.n_rows;
    arma::uword q = coefs.n_cols;

    // Number of existing constraint columns:
    arma::uword ndconstr = x.dconstr.n_cols;

    // Each column of the output has N.n_cols constraints:
    arma::uword nconstraints = N.n_cols * q;

    x.dconstr.resize(x.transparameters.n_elem,
                     ndconstr + nconstraints);

    for(arma::uword j = 0L; j < q; ++j) {
      for(arma::uword k = 0L; k < N.n_cols; ++k) {
        arma::uword column = ndconstr + j*N.n_cols + k;
        for(arma::uword i = 0L; i < n; ++i) {
          arma::uword output_index = j*n + i;
          x.dconstr(indices_out[output_index], column) = N(i, k);
        }
      }
    }

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
