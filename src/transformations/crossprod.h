/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/07/2025
 */

// Crossproduct transformation:

class crossprod:public transformations {

public:

  int p, q;
  arma::mat X;

  void transform() {

    std::memcpy(X.memptr(), transparameters.memptr(),
                transparameters.n_elem * sizeof(double));
    arma::mat XtX = arma::symmatu(X.t() * X);
    transparameters = arma::vectorise(XtX);

  }

  void jacobian() {

    // jacob = arma::diagmat(transparameters);
    // jacob = jacob.cols(vector_indices);
    // g = jacob * grad;
    arma::mat gradmat(p, q);
    std::memcpy(gradmat.memptr(), grad.memptr(), grad.n_elem * sizeof(double));
    grad = arma::vectorise(2*X * gradmat);
    // grad = 2*grad % arma::vectorise(X.t());
    // grad = 2*grad % arma::vectorise(X);
    // grad = arma::vectorise((X + X.t()) * gradmat);
    // grad = arma::vectorise(X * gradmat + X.t() * gradmat);

  }

  void d2jacobian() {

  }

  void dconstraints() {

  }

  void outcomes() {

    matrices.resize(1);
    matrices[0] = jacob;

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
