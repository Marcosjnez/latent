/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 10/11/2025
 */

// factor_cor transformation:

class factor_cor:public transformations {

public:

  int p, q;
  arma::mat X, grad_out, dX, dXtX, Dp;
  arma::uvec lower_diag;

  arma::mat R, Rhat, lambda, psi, theta, lambda_psi, glambda, gpsi, gtheta,
  dlambda, dglambda, dpsi, dRhat, dgpsi, dtheta, dgtheta;
  arma::uvec lambda_indices, psi_indices, theta_indices,
  diag_psi, diag_theta;
  arma::uvec lower_psi, lower_theta;

  void transform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(lambda_indices), p, q);
    psi.elem(lower_psi) = x.transparameters(psi_indices);
    theta.elem(lower_theta) = x.transparameters(theta_indices);
    psi = arma::symmatl(psi);
    theta = arma::symmatl(theta);

    Rhat = lambda * psi * lambda.t() + theta;
    x.transparameters(indices_out[0]) = arma::vectorise(Rhat.elem(lower_diag));

  }

  void update_grad(arguments_optim& x) {

    grad_out.elem(lower_diag) = x.grad(indices_out[0]);
    grad_out = arma::symmatl(grad_out);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    lambda_psi = lambda * psi;
    glambda = 2*grad_out * lambda_psi;

    gpsi = 2*lambda.t() * grad_out * lambda;
    gpsi.diag() *= 0.5;

    gtheta = 2*grad_out;
    gtheta.diag() *= 0.5;

    x.grad.elem(lambda_indices) += arma::vectorise(glambda);
    x.grad.elem(psi_indices) += arma::vectorise(gpsi(lower_psi));
    x.grad.elem(theta_indices) += arma::vectorise(gtheta(lower_theta));

    // arma::mat I_p = arma::eye(p, p);
    //
    // arma::mat J_lambda = arma::kron(lambda_psi.t(), I_p);
    // arma::mat J_psi    = arma::kron(lambda.t(), lambda.t());
    // arma::mat J_theta  = arma::eye(p*p, p*p);
    //
    // // Halve the diagonal entries:
    // J_psi.rows(diag_psi) *= 0.5;
    // J_theta.rows(diag_theta) *= 0.5;
    //
    // // Remove duplicated entries for psi and theta:
    // J_psi = J_psi.rows(lower_psi);
    // J_theta = J_theta.rows(lower_theta);
    //
    // // Multiply the jacobians by the duplication matrix:
    // jacob = (2 * arma::join_cols(J_lambda, J_psi, J_theta) * Dp).t();
    //
    // arma::uvec all_idx =
    //   arma::join_cols(
    //     arma::join_cols(lambda_indices, psi_indices),
    //     theta_indices
    //   );
    // x.grad.elem(all_idx) += jacob.t() * x.grad(indices_out[0]);

    // J_lambda *= Dp;
    // J_psi *= Dp;
    // J_theta *= Dp;
    // x.grad.elem(lambda_indices) += 2*J_lambda*x.grad(indices_out[0]);
    // x.grad.elem(psi_indices) += 2*J_psi*x.grad(indices_out[0]);
    // x.grad.elem(theta_indices) += 2*J_theta*x.grad(indices_out[0]);

  }

  void update_dparam(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(lambda_indices), p, q);
    dpsi.elem(lower_psi) = x.dtransparameters(psi_indices);
    dtheta.elem(lower_theta) = x.dtransparameters(theta_indices);
    dpsi = arma::symmatl(dpsi);
    dtheta = arma::symmatl(dtheta);

    dRhat = dlambda * psi * lambda.t() +
            lambda * psi * dlambda.t() +
            lambda * dpsi * lambda.t() +
            dtheta;

    x.dtransparameters(indices_out[0]) = arma::vectorise(dRhat.elem(lower_diag));

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out(p, p, arma::fill::zeros);
    dgrad_out.elem(lower_diag) = x.dgrad.elem(indices_out[0]);
    dgrad_out = arma::symmatl(dgrad_out);
    dgrad_out *= 0.50; // Do not double-count the symmetric part
    dgrad_out.diag() *= 2; // Restore the diagonal

    arma::mat dB = dlambda * psi + lambda * dpsi;
    dglambda = 2.0*(dgrad_out * lambda_psi + grad_out * dB);

    // dgpsi:
    dgpsi = 2.0*(dlambda.t() * grad_out * lambda +
                 lambda.t() * dgrad_out * lambda +
                 lambda.t() * grad_out * dlambda);
    dgpsi.diag() *= 0.5;

    // dgtheta:
    dgtheta = 2*dgrad_out;
    dgtheta.diag() *= 0.5;

    x.dgrad.elem(lambda_indices) += arma::vectorise(dglambda);
    x.dgrad.elem(psi_indices) += arma::vectorise(dgpsi(lower_psi));
    x.dgrad.elem(theta_indices) += arma::vectorise(dgtheta(lower_theta));

  }

  void jacobian(arguments_optim& x) {

    arma::mat I_p = arma::eye(p, p);

    arma::mat J_lambda = arma::kron(lambda_psi.t(), I_p);
    arma::mat J_psi    = arma::kron(lambda.t(), lambda.t());
    arma::mat J_theta  = arma::eye(p*p, p*p);

    // Halve the diagonal entries:
    J_psi.rows(diag_psi) *= 0.5;
    J_theta.rows(diag_theta) *= 0.5;

    // Remove duplicated entries for psi and theta:
    J_psi = J_psi.rows(lower_psi);
    J_theta = J_theta.rows(lower_theta);

    // Multiply the jacobians by a duplication matrix:
    jacob = (2 * arma::join_cols(J_lambda, J_psi, J_theta) * Dp).t();

  }

  void update_hess(arguments_optim& x) {

    Rf_error("sum_djacob not available");
    // sum_djacob = 2*arma::diagmat(arma::vectorise(grad_out));

  }

  void update_vcov(arguments_optim& x) {

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void outcomes(arguments_optim& x) {

    vectors.resize(1);
    vectors[0] = theta.diag();

    matrices.resize(5);
    matrices[0] = jacob;
    matrices[1] = lambda;
    matrices[2] = psi;
    matrices[3] = theta;
    matrices[4] = Rhat;

  }

};

factor_cor* choose_factor_cor(const Rcpp::List& trans_setup) {

  factor_cor* mytrans = new factor_cor();

  std::vector<arma::uvec> indices_in = trans_setup["indices_in"];
  std::vector<arma::uvec> indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::mat grad_out(p, p, arma::fill::zeros);

  arma::mat lambda(p, q, arma::fill::zeros);
  arma::mat psi(q, q, arma::fill::zeros);
  arma::mat theta(p, p, arma::fill::zeros);
  arma::mat Dp = duplication(p, true);

  arma::uvec lambda_indices = indices_in[1];
  arma::uvec psi_indices = indices_in[2];
  arma::uvec theta_indices = indices_in[3];

  arma::uvec lower_psi = arma::trimatl_ind(arma::size(psi));
  arma::uvec lower_theta = arma::trimatl_ind(arma::size(theta));
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));

  arma::uvec diag_psi   = arma::regspace<arma::uvec>(0, q - 1) * q
  + arma::regspace<arma::uvec>(0, q - 1);
  arma::uvec diag_theta = arma::regspace<arma::uvec>(0, p - 1) * p
  + arma::regspace<arma::uvec>(0, p - 1);

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->grad_out = grad_out;
  mytrans->lambda = lambda;
  mytrans->psi = psi;
  mytrans->theta = theta;
  mytrans->Dp = Dp;
  mytrans->dpsi = psi;
  mytrans->dtheta = theta;
  mytrans->lambda_indices = lambda_indices;
  mytrans->psi_indices = psi_indices;
  mytrans->theta_indices = theta_indices;
  mytrans->lower_psi = lower_psi;
  mytrans->lower_theta = lower_theta;
  mytrans->lower_diag = lower_diag;
  mytrans->diag_psi = diag_psi;
  mytrans->diag_theta = diag_theta;

  return mytrans;

}
