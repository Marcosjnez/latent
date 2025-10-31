/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 30/10/2025
 */

// factor_cor transformation:

class factor_cor:public transformations {

public:

  int p, q;
  arma::mat X, grad_out, dX, dXtX, Dp;
  arma::uvec lower_diag;

  arma::mat R, Rhat, lambda, psi, theta, lambda_psi, glambda, gpsi, gtheta,
  dlambda, dglambda, dpsi, dRhat, dgpsi, dtheta, dgtheta;
  arma::uvec lambda_indices, psi_indices, theta_indices;
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

  void update_hess(arguments_optim& x) {

    // Rf_error("wrong sum_djacob");
    // arma::mat I(p, p, arma::fill::eye);
    // jacob = 2*Dp.t() * arma::kron(I, X.t());
    // int pp = p*p;
    // int q = 0.5*p*(p-1);
    //
    // sum_djacob.resize(pp, pp);
    // sum_djacob.zeros();
    // arma::mat dX(p, p, arma::fill::zeros);
    //
    // sum_djacob = 2*arma::diagmat(arma::vectorise(grad_out));

  }

  void update_vcov(arguments_optim& x) {

    // arma::mat I(p, p, arma::fill::eye);
    // jacob = 2*Dp.t() * arma::kron(I, X.t());

  }

  void dconstraints(arguments_optim& x) {

    constraints = false;

  }

  void M(arguments_optim& x) {

  }

  void outcomes(arguments_optim& x) {

    vectors.resize(1);
    vectors[0] = theta.diag();

    matrices.resize(4);
    matrices[0] = lambda;
    matrices[1] = psi;
    matrices[2] = theta;
    matrices[3] = Rhat;

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

  arma::uvec lambda_indices = indices_in[1];
  arma::uvec psi_indices = indices_in[2];
  arma::uvec theta_indices = indices_in[3];

  arma::uvec lower_psi = arma::trimatl_ind(arma::size(psi));
  arma::uvec lower_theta = arma::trimatl_ind(arma::size(theta));
  arma::uvec lower_diag = arma::trimatl_ind(arma::size(grad_out));

  mytrans->indices_in = indices_in;
  mytrans->indices_out = indices_out;
  mytrans->p = p;
  mytrans->q = q;
  mytrans->grad_out = grad_out;
  mytrans->lambda = lambda;
  mytrans->psi = psi;
  mytrans->theta = theta;
  mytrans->dpsi = psi;
  mytrans->dtheta = theta;
  mytrans->lambda_indices = lambda_indices;
  mytrans->psi_indices = psi_indices;
  mytrans->theta_indices = theta_indices;
  mytrans->lower_psi = lower_psi;
  mytrans->lower_theta = lower_theta;
  mytrans->lower_diag = lower_diag;

  return mytrans;

}
