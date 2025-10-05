/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/09/2025
 */

// Crossproduct transformation:

class crossprod:public transformations {

public:

  int p, q;
  arma::mat X, grad_out;
  arma::uvec lower_diag;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_in[0]), p, q);
    // arma::mat XtX = arma::symmatu(X.t() * X);
    arma::mat XtX = X.t() * X;
    x.transparameters(indices_out[0]) = arma::vectorise(XtX.elem(lower_diag));

  }

  void update_grad(arguments_optim& x) {

    // grad_out = arma::reshape(x.grad(indices_out[0]), p, q);
    // grad_out *= 0.50; // Do not double-count the symmetric part
    // grad_out.diag() *= 2; // Restore the diagonal
    grad_out.elem(lower_diag) = x.grad(indices_out[0]);
    grad_out = arma::symmatl(grad_out);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal
    x.grad(indices_in[0]) += arma::vectorise(2*X * grad_out);

  }

  void update_dgrad(arguments_optim& x) {

  }

  void update_hess(arguments_optim& x) {

    jacob = arma::diagmat(arma::vectorise(2*X));
    jacob = jacob.rows(lower_diag);
    sum_djacob = arma::diagmat(2*arma::vectorise(grad_out));

  }

  void update_vcov(arguments_optim& x) {

    jacob = arma::diagmat(arma::vectorise(2*X));
    jacob = jacob.rows(lower_diag);

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
  arma::mat grad_out(q, q, arma::fill::zeros);
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->X = X;
  mytrans->grad_out = grad_out;
  mytrans->lower_diag = lower_diag;

  return mytrans;

}
