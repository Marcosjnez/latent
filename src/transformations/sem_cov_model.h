/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/09/2026
 */

// sem_cov_model transformation:

class sem_cov_model:public transformations {

public:

  bool constraints;
  int p, q;
  arma::uvec indices_lambda, indices_latent_cov, indices_theta,
  diag_latent_cov, diag_theta, lower_latent_cov, lower_theta, lower_diag;
  arma::mat R, Shat, lambda, latent_cov, theta, lambda_latent_cov, glambda,
  glatent_cov, gtheta,
  dlambda, dglambda, dlatent_cov, dShat, dglatent_cov, dtheta, dgtheta, grad_out;

  void transform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);
    latent_cov = arma::reshape(x.transparameters(indices_latent_cov), q, q);
    theta = arma::reshape(x.transparameters(indices_theta), p, p);

    Shat = lambda * latent_cov * lambda.t() + theta;
    x.transparameters(indices_out) = arma::vectorise(Shat);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out), p, p);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    lambda_latent_cov = lambda * latent_cov;
    glambda = 2*grad_out * lambda_latent_cov;
    glatent_cov = lambda.t() * grad_out * lambda;
    gtheta = grad_out;

    x.grad.elem(indices_lambda) += arma::vectorise(glambda);
    x.grad.elem(indices_latent_cov) += arma::vectorise(glatent_cov);
    x.grad.elem(indices_theta) += arma::vectorise(gtheta);

  }

  void dtransform(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    dlatent_cov = arma::reshape(x.dtransparameters(indices_latent_cov), q, q);
    dtheta = arma::reshape(x.dtransparameters(indices_theta), p, p);

    dShat = dlambda * latent_cov * lambda.t() +
            lambda * latent_cov * dlambda.t() +
            lambda * dlatent_cov * lambda.t() +
            dtheta;

    x.dtransparameters(indices_out) = arma::vectorise(dShat);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad.elem(indices_out), p, p);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    arma::mat dB = dlambda * latent_cov + lambda * dlatent_cov;
    dglambda = 2*(dgrad_out * lambda_latent_cov + grad_out * dB);

    // dglatent_cov:
    dglatent_cov = dlambda.t() * grad_out * lambda +
      lambda.t() * dgrad_out * lambda +
      lambda.t() * grad_out * dlambda;

    // dgtheta:
    dgtheta = dgrad_out;

    x.dgrad.elem(indices_lambda) += arma::vectorise(dglambda);
    x.dgrad.elem(indices_latent_cov) += arma::vectorise(dglatent_cov);
    x.dgrad.elem(indices_theta) += arma::vectorise(dgtheta);

  }

  void jacobian(arguments_optim& x) {

    lambda_latent_cov = lambda * latent_cov;

    arma::mat I_p = arma::eye(p, p);
    arma::mat comm = dxt(p, q); // Commutation matrix

    arma::mat J_lambda = arma::kron(lambda_latent_cov, I_p) +
      arma::kron(I_p, lambda_latent_cov) * comm;
    arma::mat J_latent_cov    = arma::kron(lambda, lambda);
    arma::mat J_theta  = arma::eye(p*p, p*p);

    arma::mat Dp = duplication(p, true, false);
    jacob = Dp.t() * arma::join_rows(J_lambda, J_latent_cov, J_theta);

  }

  void outcomes(arguments_optim& x) {

    matrices.resize(1);
    matrices[0] = jacob;
    names_matrices.resize(1);
    names_matrices[0] = "jacobian";

  }

};

sem_cov_model* choose_sem_cov_model(const Rcpp::List& trans_setup) {

  sem_cov_model* mytrans = new sem_cov_model();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::uvec indices_lambda = indices_in[0];
  arma::uvec indices_latent_cov = indices_in[1];
  arma::uvec indices_theta = indices_in[2];
  arma::mat grad_out(p, p, arma::fill::zeros);

  arma::mat lambda(p, q, arma::fill::zeros);
  arma::mat latent_cov(q, q, arma::fill::zeros);
  arma::mat theta(p, p, arma::fill::zeros);

  arma::uvec lower_latent_cov = arma::trimatl_ind(arma::size(latent_cov));
  arma::uvec lower_theta = arma::trimatl_ind(arma::size(theta));
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));

  arma::uvec diag_latent_cov   = arma::regspace<arma::uvec>(0, q - 1) * q
  + arma::regspace<arma::uvec>(0, q - 1);
  arma::uvec diag_theta = arma::regspace<arma::uvec>(0, p - 1) * p
  + arma::regspace<arma::uvec>(0, p - 1);


  mytrans->indices_lambda = indices_lambda;
  mytrans->indices_latent_cov = indices_latent_cov;
  mytrans->indices_theta = indices_theta;
  mytrans->indices_in = arma::join_cols(indices_lambda, indices_latent_cov,
                                        indices_theta);
  mytrans->indices_out = indices_out[0];
  mytrans->p = p;
  mytrans->q = q;
  mytrans->grad_out = grad_out;
  mytrans->lambda = lambda;
  mytrans->latent_cov = latent_cov;
  mytrans->theta = theta;
  mytrans->dlatent_cov = latent_cov;
  mytrans->dtheta = theta;
  mytrans->lower_latent_cov = lower_latent_cov;
  mytrans->lower_theta = lower_theta;
  mytrans->lower_diag = lower_diag;
  mytrans->diag_latent_cov = diag_latent_cov;
  mytrans->diag_theta = diag_theta;

  return mytrans;

}
