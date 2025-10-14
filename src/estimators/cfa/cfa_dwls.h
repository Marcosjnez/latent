/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 14/10/2025
 */

/*
 * Confirmatory factor analysis (least-squares, weigthed and unweighted)
 */

class cfa_dwls: public estimators {

public:

  arma::mat R, Rhat, lambda, psi, theta, residuals, model, W,
  lambda_psi, W_residuals, glambda, W_residuals_lambda, gpsi, gtheta,
  dlambda, dglambda, dpsi, dgpsi, dtheta, dgtheta;
  arma::uvec lambda_indices, psi_indices, theta_indices;
  arma::uvec lower_psi, lower_theta;
  int p, q;

  void param(arguments_optim& x) {

    lambda = arma::reshape(x.transparameters(lambda_indices), p, q);
    psi.elem(lower_psi) = x.transparameters(psi_indices);
    theta.elem(lower_theta) = x.transparameters(theta_indices);
    psi = arma::symmatl(psi);
    theta = arma::symmatl(theta);

    Rhat = lambda * psi * lambda.t() + theta;
    residuals = R - Rhat;

  }

  void F(arguments_optim& x) {

    f = 0.5*arma::accu(residuals % residuals % W);
    x.f += f;

  }

  void G(arguments_optim& x) {

    lambda_psi = lambda * psi;
    W_residuals = W % residuals;
    glambda = -2* W_residuals * lambda_psi;

    W_residuals_lambda = W_residuals * lambda;
    gpsi = -2*lambda.t() * W_residuals_lambda;
    gpsi.diag() *= 0.5;

    gtheta = -2*W_residuals;
    gtheta.diag() *= 0.5;

    x.grad.elem(lambda_indices) += arma::vectorise(glambda);
    x.grad.elem(psi_indices) += arma::vectorise(gpsi(lower_psi));
    x.grad.elem(theta_indices) += arma::vectorise(gtheta(lower_theta));

  }

  void dG(arguments_optim& x) {

    dlambda = arma::reshape(x.dtransparameters(lambda_indices), p, q);
    dpsi.elem(lower_psi) = x.dtransparameters(psi_indices);
    dtheta.elem(lower_theta) = x.dtransparameters(theta_indices);
    dpsi = arma::symmatl(dpsi);
    dtheta = arma::symmatl(dtheta);

    // dglambda:
    // arma::mat dg1 = -2*W_residuals * dlambda * psi;
    // arma::mat W_dresiduals = -W % (dlambda * lambda_psi.t() +
    //   lambda_psi * dlambda.t());
    // arma::mat dg2 = -2*W_dresiduals * lambda_psi;
    // dglambda = dg1 + dg2;

    // dglambda:
    arma::mat dR = -( dlambda * psi * lambda.t() +
                      lambda * psi * dlambda.t() +
                      lambda * dpsi * lambda.t() +
                      dtheta );
    arma::mat dA = W % dR;
    arma::mat dB = dlambda * psi + lambda * dpsi;
    dglambda = -2.0 * ( dA * lambda_psi + W_residuals * dB );

    // dgpsi:
    // W_dresiduals = -W % (lambda * dpsi * lambda.t());
    // dgpsi = -2*lambda.t() * W_dresiduals * lambda;
    // dgpsi.diag() *= 0.5;

    // dgpsi:
    dgpsi = -2.0 * ( dlambda.t() * W_residuals * lambda +
                     lambda.t() * dA * lambda +
                     lambda.t() * W_residuals * dlambda );
    dgpsi.diag() *= 0.5;

    // dgtheta:
    // dgtheta = 2*W % dtheta;
    // dgtheta.diag() *= 0.5;

    // dgtheta:
    dgtheta = -2*(W % dR);
    dgtheta.diag() *= 0.5;

    x.dgrad.elem(lambda_indices) += arma::vectorise(dglambda);
    x.dgrad.elem(psi_indices) += arma::vectorise(dgpsi(lower_psi));
    x.dgrad.elem(theta_indices) += arma::vectorise(dgtheta(lower_theta));

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

    arma::mat Ip(p, p, arma::fill::eye);
    arma::mat Iq(q, q, arma::fill::eye);
    arma::vec w = arma::vectorise(W);
    arma::uvec indices_p = arma::trimatl_ind(arma::size(theta));
    arma::uvec indices_q = arma::trimatl_ind(arma::size(psi));

    arma::uvec indices_diag_q(q);
    for(int i=0; i < q; ++i) indices_diag_q[i] = i * q + i;
    arma::uvec indices_diag_p(p);
    for(int i=0; i < p; ++i) indices_diag_p[i] = i * p + i;

    // Lambda
    arma::mat dlambda_dRhat_W = gLRhat(lambda, psi);
    dlambda_dRhat_W.each_col() %= w;
    arma::mat lambda_psit_kron_Ip = arma::kron(lambda_psi.t(), Ip);
    arma::mat g1 = 2*lambda_psit_kron_Ip * dlambda_dRhat_W;
    arma::mat g2 = -2*arma::kron(psi, W_residuals);
    arma::mat hlambda = g1 + g2;
    x.hess(lambda_indices, lambda_indices) += hlambda;

    // Psi
    arma::mat LL = 2*arma::kron(lambda, lambda).t();
    arma::mat dpsi_dRhat_W = gPRhat(lambda, q);
    dpsi_dRhat_W.each_col() %= w;
    arma::mat hpsi = LL * dpsi_dRhat_W;
    hpsi.rows(indices_diag_q) *= 0.5;
    x.hess(psi_indices, psi_indices) += hpsi(lower_psi, lower_psi);

    // Theta
    arma::mat dtheta_dRhat_W = gURhat(p);
    dtheta_dRhat_W.each_col() %= w;
    arma::mat W2 = 2*W;
    W2.diag() *= 0.5;
    arma::vec w2 = arma::vectorise(W2);
    arma::mat htheta = arma::diagmat(w2);
    x.hess(theta_indices, theta_indices) += htheta(lower_theta, lower_theta);

    // Lambda & Psi
    g1 = 2*lambda_psit_kron_Ip * dpsi_dRhat_W;
    arma::mat g21 = -2*arma::kron(Iq, W_residuals_lambda);
    arma::mat dtg21 = g21 * dxt(q, q);
    g2 = g21 + dtg21;
    g2.cols(indices_diag_q) -= dtg21.cols(indices_diag_q);
    arma::mat dlambda_dpsi = g1 + g2;
    x.hess(lambda_indices, psi_indices) += dlambda_dpsi.cols(lower_psi);

    // Lambda & Theta
    arma::mat dlambda_dtheta = 2*lambda_psit_kron_Ip * dtheta_dRhat_W;
    x.hess(lambda_indices, theta_indices) += dlambda_dtheta.cols(lower_theta);

    // Psi & Theta
    arma::mat dtheta_dpsi = 2*arma::kron(lambda.t(), lambda.t());
    dlambda_dtheta.each_row() %= w.t();
    x.hess(psi_indices, theta_indices) += dtheta_dpsi(lower_psi, lower_theta);

    x.hess = arma::symmatu(x.hess);

    /*
     * Join all the derivatives such that
     * hLambda              dlambda_dpsi   dlambda_dtheta
     * dlambda_dpsi.t()     hpsi           dtheta_dpsi.t()
     * dlambda_dtheta.t()   dtheta_dpsi    htheta
     */

    // x.dLPU_dS.set_size(x.parameters.n_elem, x.p*x.p); x.dLPU_dS.zeros();
    //
    // // Rcpp::Rcout << "dlambda_dS" << std::endl;
    // arma::mat g1 = -2*x.lambda_phit_kron_Ip;
    // g1.each_row() %= x.w.t();
    // arma::mat g2 = g1 * dxt(x.p, x.p);
    // arma::mat g = g1 + g2;
    // g.cols(x.indexes_diag_p) *= 0.5;
    // x.dlambda_dS = g; //(x.lambda_indexes, x.S_indexes);
    // x.dLPU_dS.rows(x.lambda_indexes) += x.dlambda_dS.rows(x.target_indexes);
    //
    // // Rcpp::Rcout << "dphi_dS" << std::endl;
    // g1 = -2*arma::kron(x.lambda.t(), x.lambda.t());
    // g1.each_row() %= x.w.t();
    // g2 = g1 * dxt(x.p, x.p);
    // g = g1 + g2;
    // g.cols(x.indexes_diag_p) *= 0.5;
    // g.rows(x.indexes_diag_q) *= 0.5;
    // x.dphi_dS = g; //(x.phi_indexes, x.S_indexes);
    // x.dLPU_dS.rows(x.phi_indexes) += x.dphi_dS.rows(x.targetphi_indexes);
    //
    // // Rcpp::Rcout << "dpsi_dS" << std::endl;
    // g = -2*arma::diagmat(x.w);
    // g.cols(x.indexes_diag_p) *= 0.5;
    // x.dpsi_dS = g; //(x.psi_indexes, x.S_indexes);
    // x.dLPU_dS.rows(x.psi_indexes) += x.dpsi_dS.rows(x.targetpsi_indexes);

  }

  void outcomes(arguments_optim& x) {

    // x.uniquenesses = x.R.diag() - arma::diagvec(x.Rhat);
    // x.Rhat.diag() = x.R.diag();

    doubles.resize(1);

    vectors.resize(1);
    vectors[0] = theta.diag();

    matrices.resize(6);
    matrices[0] = lambda;
    matrices[1] = psi;
    matrices[2] = theta;
    matrices[3] = Rhat;
    matrices[4] = residuals;
    matrices[5] = W;

    cubes.resize(1);
    list_matrices.resize(1);

  };

};

cfa_dwls* choose_cfa_dwls(const Rcpp::List& estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat R = estimator_setup["R"];
  arma::mat W = estimator_setup["W"];
  int nfactors = estimator_setup["nfactors"];

  int p = R.n_rows;
  int q = nfactors;
  arma::mat lambda(p, q, arma::fill::zeros);
  arma::mat psi(q, q, arma::fill::zeros);
  arma::mat theta(p, p, arma::fill::zeros);

  arma::uvec lambda_indices = indices[1];
  arma::uvec psi_indices = indices[2];
  arma::uvec theta_indices = indices[3];

  arma::uvec lower_psi = arma::trimatl_ind(arma::size(psi));
  arma::uvec lower_theta = arma::trimatl_ind(arma::size(theta));

  myestimator->indices = indices;
  myestimator->R = R;
  myestimator->W = W;
  myestimator->lambda = lambda;
  myestimator->psi = psi;
  myestimator->theta = theta;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->lambda_indices = lambda_indices;
  myestimator->psi_indices = psi_indices;
  myestimator->theta_indices = theta_indices;
  myestimator->lower_psi = lower_psi;
  myestimator->lower_theta = lower_theta;
  myestimator->dlambda = lambda;
  myestimator->dpsi = psi;
  myestimator->dtheta = theta;

  return myestimator;

}
