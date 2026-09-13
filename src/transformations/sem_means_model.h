/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

/*
 * Model-implied observed mean transformation:
 *
 * meanshat = nu + lambda*latent_means
 *
 * nu       : p observed-variable intercepts
 * lambda   : p by q factor-loading matrix
 * latent_means    : q latent-factor means
 * meanshat : p model-implied observed means
 */

class sem_means_model: public transformations {

public:

  int p, q;
  arma::uvec indices_nu, indices_lambda, indices_latent_means;
  arma::vec nu, latent_means, meanshat, dnu, dlatent_means, dmeanshat, grad_out;
  arma::mat lambda, dlambda;

  void transform(arguments_optim& x) {

    nu = x.transparameters.elem(indices_nu);
    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    latent_means = x.transparameters.elem(indices_latent_means);

    meanshat = nu + lambda * latent_means;
    x.transparameters.elem(indices_out) = meanshat;

  }

  void update_grad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    latent_means = x.transparameters.elem(indices_latent_means);
    grad_out = x.grad.elem(indices_out);

    arma::vec grad_nu = grad_out;
    arma::mat grad_lambda = grad_out*latent_means.t();
    arma::vec grad_latent_means = lambda.t()*grad_out;

    x.grad.elem(indices_nu) += grad_nu;
    x.grad.elem(indices_lambda) += arma::vectorise(grad_lambda);
    x.grad.elem(indices_latent_means) += grad_latent_means;

  }

  void dtransform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    latent_means = x.transparameters.elem(indices_latent_means);
    dnu = x.dtransparameters.elem(indices_nu);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dlatent_means = x.dtransparameters.elem(indices_latent_means);

    dmeanshat = dnu+dlambda*latent_means+lambda*dlatent_means;
    x.dtransparameters.elem(indices_out) = dmeanshat;

  }

  void update_dgrad(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    latent_means = x.transparameters.elem(indices_latent_means);
    dlambda = arma::reshape(x.dtransparameters.elem(indices_lambda), p, q);
    dlatent_means = x.dtransparameters.elem(indices_latent_means);
    grad_out = x.grad.elem(indices_out);
    arma::vec dgrad_out = x.dgrad.elem(indices_out);

    arma::vec dgrad_nu = dgrad_out;
    arma::mat dgrad_lambda =
      dgrad_out*latent_means.t()+grad_out*dlatent_means.t();
    arma::vec dgrad_latent_means =
      dlambda.t()*grad_out+lambda.t()*dgrad_out;

    x.dgrad.elem(indices_nu) += dgrad_nu;
    x.dgrad.elem(indices_lambda) += arma::vectorise(dgrad_lambda);
    x.dgrad.elem(indices_latent_means) += dgrad_latent_means;

  }

  void jacobian(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters.elem(indices_lambda), p, q);
    latent_means = x.transparameters.elem(indices_latent_means);

    arma::mat J_nu = arma::eye(p, p);
    arma::mat J_lambda = arma::kron(latent_means.t(), arma::eye(p, p));
    arma::mat J_latent_means = lambda;

    jacob = arma::join_rows(J_nu,
                            arma::join_rows(J_lambda, J_latent_means));

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
    Rcpp::stop("sem_means_model requires nu, lambda, and latent_means inputs and one meanshat output.");
  }

  arma::uvec indices_nu = indices_in[0];
  arma::uvec indices_lambda = indices_in[1];
  arma::uvec indices_latent_means = indices_in[2];
  arma::uvec indices_meanshat = indices_out[0];

  if(indices_nu.n_elem != static_cast<arma::uword>(p) ||
     indices_lambda.n_elem != static_cast<arma::uword>(p*q) ||
     indices_latent_means.n_elem != static_cast<arma::uword>(q) ||
     indices_meanshat.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("The sem_means_model parameter indices have incompatible dimensions.");
  }

  mytrans->indices_nu = indices_nu;
  mytrans->indices_lambda = indices_lambda;
  mytrans->indices_latent_means = indices_latent_means;
  mytrans->indices_in = arma::join_cols(
    indices_nu,
    arma::join_cols(indices_lambda, indices_latent_means)
  );
  mytrans->indices_out = indices_meanshat;
  mytrans->p = p;
  mytrans->q = q;

  return mytrans;

}
