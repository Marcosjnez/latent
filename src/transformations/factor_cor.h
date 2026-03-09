/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/03/2026
 */

// factor_cor transformation:

class factor_cor:public transformations {

public:

  bool constraints;
  int p, q;
  arma::uvec indices_lambda, indices_psi, indices_theta, indices_in,
  indices_out, diag_psi, diag_theta, lower_psi, lower_theta, lower_diag;
  arma::mat R, Rhat, lambda, psi, theta, lambda_psi, glambda, gpsi, gtheta,
  dlambda, dglambda, dpsi, dRhat, dgpsi, dtheta, dgtheta, grad_out, jacob;

  void transform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);
    psi = arma::reshape(x.transparameters(indices_psi), q, q);
    theta = arma::reshape(x.transparameters(indices_theta), p, p);

    Rhat = lambda * psi * lambda.t() + theta;
    x.transparameters(indices_out) = arma::vectorise(Rhat);

  }

  void update_grad(arguments_optim& x) {

    grad_out = arma::reshape(x.grad(indices_out), p, p);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    lambda_psi = lambda * psi;
    glambda = 2*grad_out * lambda_psi;
    gpsi = lambda.t() * grad_out * lambda;
    gtheta = grad_out;

    x.grad.elem(indices_lambda) += arma::vectorise(glambda);
    x.grad.elem(indices_psi) += arma::vectorise(gpsi);
    x.grad.elem(indices_theta) += arma::vectorise(gtheta);

    // arma::mat I_p = arma::eye(p, p);
    // arma::mat comm = dxt(p, q);
    //
    // arma::mat J_lambda = arma::kron(lambda_psi, I_p) +
    //   arma::kron(I_p, lambda_psi) * comm;
    // arma::mat J_psi    = arma::kron(lambda, lambda);
    // arma::mat J_theta  = arma::eye(p*p, p*p);
    //
    // arma::mat Dp = duplication(p, true, false);
    // jacob = Dp.t() * arma::join_rows(J_lambda, J_psi, J_theta);
    //
    // arma::uvec all_idx =
    //   arma::join_cols(
    //     arma::join_cols(indices_lambda, indices_psi),
    //     indices_theta
    //   );
    // x.grad.elem(all_idx) += jacob.t() * x.grad(indices_out);

  }

  void dtransform(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    dpsi = arma::reshape(x.dtransparameters(indices_psi), q, q);
    dtheta = arma::reshape(x.dtransparameters(indices_theta), p, p);

    dRhat = dlambda * psi * lambda.t() +
            lambda * psi * dlambda.t() +
            lambda * dpsi * lambda.t() +
            dtheta;

    x.dtransparameters(indices_out) = arma::vectorise(dRhat);

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out = arma::reshape(x.dgrad.elem(indices_out), p, p);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    arma::mat dB = dlambda * psi + lambda * dpsi;
    dglambda = 2*(dgrad_out * lambda_psi + grad_out * dB);

    // dgpsi:
    dgpsi = dlambda.t() * grad_out * lambda +
      lambda.t() * dgrad_out * lambda +
      lambda.t() * grad_out * dlambda;

    // dgtheta:
    dgtheta = dgrad_out;

    x.dgrad.elem(indices_lambda) += arma::vectorise(dglambda);
    x.dgrad.elem(indices_psi) += arma::vectorise(dgpsi);
    x.dgrad.elem(indices_theta) += arma::vectorise(dgtheta);

  }

  void jacobian(arguments_optim& x) {

    arma::mat I_p = arma::eye(p, p);
    arma::mat comm = dxt(p, q); // Commutation matrix

    arma::mat J_lambda = arma::kron(lambda_psi, I_p) +
      arma::kron(I_p, lambda_psi) * comm;
    arma::mat J_psi    = arma::kron(lambda, lambda);
    arma::mat J_theta  = arma::eye(p*p, p*p);

    arma::mat Dp = duplication(p, true, false);
    jacob = Dp.t() * arma::join_rows(J_lambda, J_psi, J_theta);

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(indices_lambda, indices_psi, indices_theta);
    x.vcov(indices_out, indices_out).zeros();
    x.vcov(indices_out, indices_out) += jacob * x.vcov(indices_in, indices_in) * jacob.t();

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    vectors.resize(1);
    vectors[0] = theta.diag();

    matrices.resize(1);
    matrices[0] = jacob;

  }

};

factor_cor* choose_factor_cor(const Rcpp::List& trans_setup) {

  factor_cor* mytrans = new factor_cor();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::uvec indices_lambda = indices_in[0];
  arma::uvec indices_psi = indices_in[1];
  arma::uvec indices_theta = indices_in[2];
  arma::mat grad_out(p, p, arma::fill::zeros);

  arma::mat lambda(p, q, arma::fill::zeros);
  arma::mat psi(q, q, arma::fill::zeros);
  arma::mat theta(p, p, arma::fill::zeros);

  arma::uvec lower_psi = arma::trimatl_ind(arma::size(psi));
  arma::uvec lower_theta = arma::trimatl_ind(arma::size(theta));
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));

  arma::uvec diag_psi   = arma::regspace<arma::uvec>(0, q - 1) * q
  + arma::regspace<arma::uvec>(0, q - 1);
  arma::uvec diag_theta = arma::regspace<arma::uvec>(0, p - 1) * p
  + arma::regspace<arma::uvec>(0, p - 1);

  mytrans->indices_lambda = indices_lambda;
  mytrans->indices_psi = indices_psi;
  mytrans->indices_theta = indices_theta;
  mytrans->indices_out = indices_out[0];
  mytrans->p = p;
  mytrans->q = q;
  mytrans->grad_out = grad_out;
  mytrans->lambda = lambda;
  mytrans->psi = psi;
  mytrans->theta = theta;
  mytrans->dpsi = psi;
  mytrans->dtheta = theta;
  mytrans->lower_psi = lower_psi;
  mytrans->lower_theta = lower_theta;
  mytrans->lower_diag = lower_diag;
  mytrans->diag_psi = diag_psi;
  mytrans->diag_theta = diag_theta;

  return mytrans;

}
