/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

// X*Y*z transformation:
// X is p x q, Y is q x r, z has length r, and the output has length p.
// Derivative updates use the current transform() caches; update_dgrad()
// also requires update_grad() and dtransform() at the same parameter point.

class XYz: public transformations {

public:

  int p, q, r;
  arma::uvec indices_X, indices_Y, indices_z, indices_XYz;
  arma::mat X, Y, dX, dY, grad_in_X, grad_in_Y;
  arma::vec z, dz, Yz, dYz, dXYz, grad_out, Xt_grad_out, grad_in_z;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), q, r);
    z = x.transparameters(indices_z);

    Yz = Y*z;
    x.transparameters(indices_XYz) = X*Yz;

  }

  void update_grad(arguments_optim& x) {

    grad_out = x.grad(indices_XYz);
    Xt_grad_out = X.t()*grad_out;

    grad_in_X = grad_out*Yz.t();
    grad_in_Y = Xt_grad_out*z.t();
    grad_in_z = Y.t()*Xt_grad_out;

    x.grad(indices_X) += arma::vectorise(grad_in_X);
    x.grad(indices_Y) += arma::vectorise(grad_in_Y);
    x.grad(indices_z) += grad_in_z;

  }

  void dtransform(arguments_optim& x) {

    dX = arma::reshape(x.dtransparameters(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters(indices_Y), q, r);
    dz = x.dtransparameters(indices_z);

    dYz = dY*z+Y*dz;
    dXYz = dX*Yz+X*dYz;
    x.dtransparameters(indices_XYz) = dXYz;

  }

  void update_dgrad(arguments_optim& x) {

    arma::vec dgrad_out = x.dgrad(indices_XYz);
    arma::vec dXt_grad_out = dX.t()*grad_out+X.t()*dgrad_out;

    arma::mat dgrad_in_X = dgrad_out*Yz.t()+grad_out*dYz.t();
    arma::mat dgrad_in_Y = dXt_grad_out*z.t()+Xt_grad_out*dz.t();
    arma::vec dgrad_in_z = dY.t()*Xt_grad_out+Y.t()*dXt_grad_out;

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);
    x.dgrad.elem(indices_z) += dgrad_in_z;

  }

  void jacobian(arguments_optim& x) {

    X = arma::reshape(x.transparameters(indices_X), p, q);
    Y = arma::reshape(x.transparameters(indices_Y), q, r);
    z = x.transparameters(indices_z);
    Yz = Y*z;

    // Column order: vec(X), vec(Y), z.
    arma::mat I_p = arma::eye(p, p);
    arma::mat Jx = arma::kron(Yz.t(), I_p);
    arma::mat Jy = arma::kron(z.t(), X);
    arma::mat Jz = X*Y;

    jacob = arma::join_rows(arma::join_rows(Jx, Jy), Jz);

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

XYz* choose_XYz(const Rcpp::List& trans_setup) {

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];
  int r = trans_setup["r"];

  if(p < 1 || q < 1 || r < 1) {
    Rcpp::stop("XYz requires positive p, q, and r dimensions.");
  }

  if(indices_in.size() != 3L || indices_out.size() != 1L) {
    Rcpp::stop("XYz requires X, Y, and z inputs and one output.");
  }

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_Y = indices_in[1];
  arma::uvec indices_z = indices_in[2];
  arma::uvec indices_XYz = indices_out[0];

  if(indices_X.n_elem != static_cast<arma::uword>(p)*q ||
     indices_Y.n_elem != static_cast<arma::uword>(q)*r ||
     indices_z.n_elem != static_cast<arma::uword>(r) ||
     indices_XYz.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The XYz parameter indices have incompatible dimensions.");
  }

  XYz* mytrans = new XYz();

  mytrans->indices_in = arma::join_cols(arma::join_cols(indices_X, indices_Y),
                                      indices_z);
  mytrans->indices_out = indices_XYz;
  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_z = indices_z;
  mytrans->indices_XYz = indices_XYz;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->r = r;

  return mytrans;

}
