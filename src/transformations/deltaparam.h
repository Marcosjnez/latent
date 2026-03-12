/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/03/2026
 */

// 1 - diag(X * Y * X.t()) transformation:

// This is the delta parameterization in factor analysis
// latent scores have a fixed variance of 1 for identification
// Error variances become 1 - diag(lambda * psi * lambda.t())

class deltaparam: public transformations {

public:

  int p, q;
  arma::uvec indices_X, indices_Y, indices_xyxt, indices_in, indices_out;
  arma::mat X, Y, dX, dY, grad_in_X, grad_in_Y, jacob;
  arma::vec xyxt, dxyxt, grad_out;

  void transform(arguments_optim& x) {

    X = arma::reshape(x.transparameters.elem(indices_X), p, q);
    Y = arma::reshape(x.transparameters.elem(indices_Y), q, q);
    xyxt = 1.0 - arma::diagvec(X * Y * X.t());

    x.transparameters.elem(indices_xyxt) = xyxt;

  }

  void update_grad(arguments_optim& x) {

    X = arma::reshape(x.transparameters.elem(indices_X), p, q);
    Y = arma::reshape(x.transparameters.elem(indices_Y), q, q);
    grad_out = x.grad.elem(indices_xyxt);

    arma::mat Dg = arma::diagmat(grad_out);

    grad_in_X = -Dg * X * (Y + Y.t());
    grad_in_Y = -X.t() * Dg * X;

    x.grad.elem(indices_X) += arma::vectorise(grad_in_X);
    x.grad.elem(indices_Y) += arma::vectorise(grad_in_Y);

  }

  void dtransform(arguments_optim& x) {

    X = arma::reshape(x.transparameters.elem(indices_X), p, q);
    Y = arma::reshape(x.transparameters.elem(indices_Y), q, q);
    dX = arma::reshape(x.dtransparameters.elem(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters.elem(indices_Y), q, q);

    arma::mat dXYXt = X * dY * X.t() + dX * Y * X.t() + X * Y * dX.t();
    dxyxt = -arma::diagvec(dXYXt);

    x.dtransparameters.elem(indices_xyxt) = dxyxt;

  }

  void update_dgrad(arguments_optim& x) {

    X = arma::reshape(x.transparameters.elem(indices_X), p, q);
    Y = arma::reshape(x.transparameters.elem(indices_Y), q, q);
    dX = arma::reshape(x.dtransparameters.elem(indices_X), p, q);
    dY = arma::reshape(x.dtransparameters.elem(indices_Y), q, q);
    grad_out = x.grad.elem(indices_xyxt);
    arma::vec dgrad_out = x.dgrad.elem(indices_xyxt);

    arma::mat Dg = arma::diagmat(grad_out);
    arma::mat Ddg = arma::diagmat(dgrad_out);

    arma::mat dgrad_in_X = -Ddg * X * (Y + Y.t()) - Dg * dX * (Y + Y.t()) - Dg * X * (dY + dY.t());
    arma::mat dgrad_in_Y = -(dX.t() * Dg * X + X.t() * Ddg * X + X.t() * Dg * dX);

    x.dgrad.elem(indices_X) += arma::vectorise(dgrad_in_X);
    x.dgrad.elem(indices_Y) += arma::vectorise(dgrad_in_Y);

  }

  void jacobian(arguments_optim& x) {

    // CHECK THE JACOBIAN

    arma::mat S = Y + Y.t();
    arma::mat Jx(p, p * q, arma::fill::zeros);
    arma::mat Jy(p, q * q, arma::fill::zeros);

    for (int i = 0; i < p; ++i) {
      arma::rowvec rx = -X.row(i) * S;
      for (int j = 0; j < q; ++j) {
        Jx(i, i + j * p) = rx(j);
      }
      Jy.row(i) = -arma::vectorise(X.row(i).t() * X.row(i)).t();
    }

    jacob = arma::join_rows(Jx, Jy);

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(indices_X, indices_Y);
    indices_out = indices_xyxt;
    x.vcov(indices_out, indices_out) = jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    vectors.resize(1);

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

deltaparam* choose_deltaparam(const Rcpp::List& trans_setup) {

  deltaparam* mytrans = new deltaparam();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::uvec indices_X = indices_in[0];
  arma::uvec indices_Y = indices_in[1];
  arma::uvec indices_xyxt = indices_out[0];

  mytrans->indices_X = indices_X;
  mytrans->indices_Y = indices_Y;
  mytrans->indices_xyxt = indices_xyxt;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
