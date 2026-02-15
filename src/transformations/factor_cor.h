/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/02/2026
 */

// factor_cor transformation:

class factor_cor:public transformations {

public:

  arma::uvec indices_lambda, indices_psi, indices_theta, indices_in,
  indices_out, diag_psi, diag_theta, lower_psi, lower_theta, lower_diag;
  int p, q;
  arma::mat R, Rhat, lambda, psi, theta, lambda_psi, glambda, gpsi, gtheta,
  dlambda, dglambda, dpsi, dRhat, dgpsi, dtheta, dgtheta, grad_out, Dp, jacob;

  void transform(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(indices_lambda), p, q);
    psi.elem(lower_psi) = x.transparameters(indices_psi);
    theta.elem(lower_theta) = x.transparameters(indices_theta);
    psi = arma::symmatl(psi);
    theta = arma::symmatl(theta);

    Rhat = lambda * psi * lambda.t() + theta;
    x.transparameters(indices_out) = arma::vectorise(Rhat.elem(lower_diag));

  }

  void update_grad(arguments_optim& x) {

    grad_out.elem(lower_diag) = x.grad(indices_out);
    grad_out = arma::symmatl(grad_out);
    grad_out *= 0.50; // Do not double-count the symmetric part
    grad_out.diag() *= 2; // Restore the diagonal

    lambda_psi = lambda * psi;
    glambda = grad_out * (2*lambda_psi);

    gpsi = 2*lambda.t() * grad_out * lambda;
    gpsi.diag() *= 0.5;

    gtheta = 2*grad_out;
    gtheta.diag() *= 0.5;

    x.grad.elem(indices_lambda) += arma::vectorise(glambda);
    x.grad.elem(indices_psi) += arma::vectorise(gpsi(lower_psi));
    x.grad.elem(indices_theta) += arma::vectorise(gtheta(lower_theta));

    // arma::mat I_p = arma::eye(p, p);
    // arma::mat comm = dxt(p, q);
    //
    // arma::mat J_lambda = arma::kron(lambda_psi, I_p) +
    //   arma::kron(I_p, lambda_psi) * comm;
    // arma::mat J_psi    = arma::kron(lambda, lambda);
    // arma::mat J_theta  = arma::eye(p*p, p*p);
    //
    // // Halve the diagonal entries:
    // J_psi.cols(diag_psi) *= 0.5;
    // J_theta.cols(diag_theta) *= 0.5;
    //
    // // Remove duplicated entries for psi and theta:
    // J_psi = 2*J_psi.cols(lower_psi);
    // J_theta = 2*J_theta.cols(lower_theta);
    //
    // // Multiply the jacobians by the duplication matrix:
    // jacob = Dp.t() * arma::join_rows(J_lambda, J_psi, J_theta);
    // // The duplication matrix halves the duplicated elements in the
    // // transformed parameters´
    //
    // arma::uvec all_idx =
    //   arma::join_cols(
    //     arma::join_cols(indices_lambda, indices_psi),
    //     indices_theta
    //   );
    // x.grad.elem(all_idx) += jacob.t() * x.grad(indices_out[0]);

  }

  void dtransform(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(indices_lambda), p, q);
    dpsi.elem(lower_psi) = x.dtransparameters(indices_psi);
    dtheta.elem(lower_theta) = x.dtransparameters(indices_theta);
    dpsi = arma::symmatl(dpsi);
    dtheta = arma::symmatl(dtheta);

    dRhat = dlambda * psi * lambda.t() +
            lambda * psi * dlambda.t() +
            lambda * dpsi * lambda.t() +
            dtheta;

    x.dtransparameters(indices_out) = arma::vectorise(dRhat.elem(lower_diag));

  }

  void update_dgrad(arguments_optim& x) {

    arma::mat dgrad_out(p, p, arma::fill::zeros);
    dgrad_out.elem(lower_diag) = x.dgrad.elem(indices_out);
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

    x.dgrad.elem(indices_lambda) += arma::vectorise(dglambda);
    x.dgrad.elem(indices_psi) += arma::vectorise(dgpsi(lower_psi));
    x.dgrad.elem(indices_theta) += arma::vectorise(dgtheta(lower_theta));

  }

  void jacobian(arguments_optim& x) {

    arma::mat I_p = arma::eye(p, p);
    arma::mat comm = dxt(p, q); // Commutation matrix

    arma::mat J_lambda = arma::kron(lambda_psi, I_p) +
      arma::kron(I_p, lambda_psi) * comm;
    arma::mat J_psi    = arma::kron(lambda, lambda);
    arma::mat J_theta  = arma::eye(p*p, p*p);

    // Halve the diagonal entries:
    J_psi.cols(diag_psi) *= 0.5;
    J_theta.cols(diag_theta) *= 0.5;

    // Remove duplicated entries for psi and theta:
    J_psi = 2*J_psi.cols(lower_psi);
    J_theta = 2*J_theta.cols(lower_theta);

    // Multiply the jacobians by the duplication matrix:
    jacob = Dp.t() * arma::join_rows(J_lambda, J_psi, J_theta);
    // The duplication matrix halves the duplicated elements in the
    // transformed parameters´

  }

  void update_vcov(arguments_optim& x) {

    indices_in = arma::join_cols(indices_lambda, indices_psi, indices_theta);
    x.vcov(indices_out, indices_out) = jacob * x.vcov(indices_in, indices_in) * jacob.t();

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

  arma::uvec indices_lambda = trans_setup["indices_lambda"];
  arma::uvec indices_psi = trans_setup["indices_psi"];
  arma::uvec indices_theta = trans_setup["indices_theta"];
  arma::uvec indices_out = trans_setup["indices_out"];
  int p = trans_setup["p"];
  int q = trans_setup["q"];

  arma::mat grad_out(p, p, arma::fill::zeros);

  arma::mat lambda(p, q, arma::fill::zeros);
  arma::mat psi(q, q, arma::fill::zeros);
  arma::mat theta(p, p, arma::fill::zeros);
  arma::mat Dp = duplication(p, true);

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
  mytrans->lower_psi = lower_psi;
  mytrans->lower_theta = lower_theta;
  mytrans->lower_diag = lower_diag;
  mytrans->diag_psi = diag_psi;
  mytrans->diag_theta = diag_theta;

  return mytrans;

}
