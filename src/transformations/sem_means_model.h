/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

/*
 * Model-implied observed mean transformation:
 *
 * meanshat = nu + lambda*mu
 *
 * nu       : p observed-variable intercepts
 * lambda   : p by q factor-loading matrix
 * mu    : q latent-factor means
 * meanshat : p model-implied observed means
 */

class sem_means_model: public transformations {

public:

  int p, q;
  arma::uvec indices_nu, indices_lambda, indices_mu;
  arma::vec nu, mu, meanshat, dnu, dmu, dmeanshat, grad_out;
  arma::mat lambda, dlambda;

  void transform(arguments_optim& x) {

    nu = x.transparameters.elem(indices_nu);
    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    mu = x.transparameters.elem(indices_mu);

    meanshat = nu + lambda * mu;
    x.transparameters.elem(indices_out) = meanshat;

  }

  void update_grad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    mu = x.transparameters.elem(indices_mu);
    grad_out = x.grad.elem(indices_out);

    arma::vec grad_nu = grad_out;
    arma::mat grad_lambda = grad_out*mu.t();
    arma::vec grad_mu = lambda.t()*grad_out;

    x.grad.elem(indices_nu) += grad_nu;
    x.grad.elem(indices_lambda) += arma::vectorise(grad_lambda);
    x.grad.elem(indices_mu) += grad_mu;

  }

  void dtransform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    mu = x.transparameters.elem(indices_mu);
    dnu = x.dtransparameters.elem(indices_nu);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dmu = x.dtransparameters.elem(indices_mu);

    dmeanshat = dnu+dlambda*mu+lambda*dmu;
    x.dtransparameters.elem(indices_out) = dmeanshat;

  }

  void update_dgrad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    mu = x.transparameters.elem(indices_mu);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dmu = x.dtransparameters.elem(indices_mu);
    grad_out = x.grad.elem(indices_out);
    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    arma::vec dgrad_nu = dgrad_out;
    arma::mat dgrad_lambda =
      dgrad_out*mu.t()+grad_out*dmu.t();
    arma::vec dgrad_mu =
      dlambda.t()*grad_out+lambda.t()*dgrad_out;

    x.dgrad.elem(indices_nu) += dgrad_nu;
    x.dgrad.elem(indices_lambda) += arma::vectorise(dgrad_lambda);
    x.dgrad.elem(indices_mu) += dgrad_mu;

  }

  void jacobian(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    mu = x.transparameters.elem(indices_mu);

    arma::mat J_nu = arma::eye(p, p);
    arma::mat J_lambda = arma::kron(mu.t(), arma::eye(p, p));
    arma::mat J_mu = lambda;

    jacob = arma::join_rows(J_nu,
                            arma::join_rows(J_lambda, J_mu));

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

sem_means_model* choose_sem_means_model(const Rcpp::List& trans_setup) {

  sem_means_model* mytrans = new sem_means_model();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  if(p < 1 || q < 1) {
    Rcpp::stop("sem_means_model requires positive p and q dimensions.");
  }

  if(indices_in.size() != 3L || indices_out.size() != 1L) {
    Rcpp::stop("sem_means_model requires nu, lambda, and mu inputs and one meanshat output.");
  }

  arma::uvec indices_nu = indices_in[0];
  arma::uvec indices_lambda = indices_in[1];
  arma::uvec indices_mu = indices_in[2];
  arma::uvec indices_meanshat = indices_out[0];

  if(indices_nu.n_elem != static_cast<arma::uword>(p) ||
     indices_lambda.n_elem != static_cast<arma::uword>(p*q) ||
     indices_mu.n_elem != static_cast<arma::uword>(q) ||
     indices_meanshat.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The sem_means_model parameter indices have incompatible dimensions.");
  }

  mytrans->indices_nu = indices_nu;
  mytrans->indices_lambda = indices_lambda;
  mytrans->indices_mu = indices_mu;
  mytrans->indices_in = arma::join_cols(
    indices_nu,
    arma::join_cols(indices_lambda, indices_mu)
  );
  mytrans->indices_out = indices_meanshat;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
