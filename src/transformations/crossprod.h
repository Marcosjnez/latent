/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 18/08/2025
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
    arma::mat gradmat(p, q);
    std::memcpy(gradmat.memptr(), grad_in.memptr(), grad_in.n_elem * sizeof(double));
    grad_out = arma::vectorise(2*X * gradmat);

  }

  void d2jacobian() {

    jacob = arma::diagmat(transparameters);

    sum_djacob.set_size(indices_in[0].n_elem, indices_in[0].n_elem);
    sum_djacob.zeros();

    hess_out.set_size(indices_in[0].n_elem, indices_in[0].n_elem);
    hess_out.zeros();

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
