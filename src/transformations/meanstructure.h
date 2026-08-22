/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 22/08/2026
 */

/*
 * Model-implied observed mean transformation:
 *
 * meanshat = nu + lambda*alpha
 *
 * nu       : p observed-variable intercepts
 * lambda   : p by q factor-loading matrix
 * alpha    : q latent-factor means
 * meanshat : p model-implied observed means
 */

class meanstructure: public transformations {

public:

  int p, q;
  arma::uvec indices_nu, indices_lambda, indices_alpha;
  arma::vec nu, alpha, meanshat, dnu, dalpha, dmeanshat, grad_out;
  arma::mat lambda, dlambda;

  void transform(arguments_optim& x) {

    nu = x.transparameters.elem(indices_nu);
    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    alpha = x.transparameters.elem(indices_alpha);

    meanshat = nu+lambda*alpha;
    x.transparameters.elem(indices_out) = meanshat;

  }

  void update_grad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    alpha = x.transparameters.elem(indices_alpha);
    grad_out = x.grad.elem(indices_out);

    arma::vec grad_nu = grad_out;
    arma::mat grad_lambda = grad_out*alpha.t();
    arma::vec grad_alpha = lambda.t()*grad_out;

    x.grad.elem(indices_nu) += grad_nu;
    x.grad.elem(indices_lambda) += arma::vectorise(grad_lambda);
    x.grad.elem(indices_alpha) += grad_alpha;

  }

  void dtransform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    alpha = x.transparameters.elem(indices_alpha);
    dnu = x.dtransparameters.elem(indices_nu);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dalpha = x.dtransparameters.elem(indices_alpha);

    dmeanshat = dnu+dlambda*alpha+lambda*dalpha;
    x.dtransparameters.elem(indices_out) = dmeanshat;

  }

  void update_dgrad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    alpha = x.transparameters.elem(indices_alpha);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dalpha = x.dtransparameters.elem(indices_alpha);
    grad_out = x.grad.elem(indices_out);
    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    arma::vec dgrad_nu = dgrad_out;
    arma::mat dgrad_lambda =
      dgrad_out*alpha.t()+grad_out*dalpha.t();
    arma::vec dgrad_alpha =
      dlambda.t()*grad_out+lambda.t()*dgrad_out;

    x.dgrad.elem(indices_nu) += dgrad_nu;
    x.dgrad.elem(indices_lambda) += arma::vectorise(dgrad_lambda);
    x.dgrad.elem(indices_alpha) += dgrad_alpha;

  }

  void jacobian(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    alpha = x.transparameters.elem(indices_alpha);

    arma::mat J_nu = arma::eye(p, p);
    arma::mat J_lambda = arma::kron(alpha.t(), arma::eye(p, p));
    arma::mat J_alpha = lambda;

    jacob = arma::join_rows(J_nu,
                            arma::join_rows(J_lambda, J_alpha));

  }

  void outcomes(arguments_optim& x) {

    (void)x;

    vectors.resize(1);
    vectors[0] = meanshat;
    names_vectors.resize(1);
    names_vectors[0] = "meanshat";

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

meanstructure* choose_meanstructure(const Rcpp::List& trans_setup) {

  meanstructure* mytrans = new meanstructure();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  if(p < 1 || q < 1) {
    Rcpp::stop("meanstructure requires positive p and q dimensions.");
  }

  if(indices_in.size() != 3L || indices_out.size() != 1L) {
    Rcpp::stop("meanstructure requires nu, lambda, and alpha inputs and one meanshat output.");
  }

  arma::uvec indices_nu = indices_in[0];
  arma::uvec indices_lambda = indices_in[1];
  arma::uvec indices_alpha = indices_in[2];
  arma::uvec indices_meanshat = indices_out[0];

  if(indices_nu.n_elem != static_cast<arma::uword>(p) ||
     indices_lambda.n_elem != static_cast<arma::uword>(p*q) ||
     indices_alpha.n_elem != static_cast<arma::uword>(q) ||
     indices_meanshat.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The meanstructure parameter indices have incompatible dimensions.");
  }

  mytrans->indices_nu = indices_nu;
  mytrans->indices_lambda = indices_lambda;
  mytrans->indices_alpha = indices_alpha;
  mytrans->indices_in = arma::join_cols(
    indices_nu,
    arma::join_cols(indices_lambda, indices_alpha)
  );
  mytrans->indices_out = indices_meanshat;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
