/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/09/2025
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
    // psi = arma::reshape(x.transparameters(psi_indices), q, q);
    // theta = arma::reshape(x.transparameters(theta_indices), p, p);
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
    // The symmetric part will be counted twice in the += operation, so
    // the correct derivatives for the nondiagonal elements in psi and theta
    // are 2*gpsi and 2*gtheta

  }

  void dG(arguments_optim& x) {

    // dX = arma::reshape(dparameters, q, q);
    dg.set_size(transparameters.n_elem); dg.zeros();

    // dlambda.fill(dparameters(lambda_indices));
    // dphi.fill(dparameters(phi_indices));
    // dpsi.fill(dparameters(psi_indices));

    dlambda = arma::reshape(dparameters(lambda_indices), p, q);
    dpsi = arma::reshape(dparameters(psi_indices), q, q);
    dtheta = arma::reshape(dparameters(theta_indices), q, q);

    dpsi = arma::symmatl(dpsi);
    dtheta = arma::symmatl(dtheta);

    // dglambda:
    arma::mat dg1 = -2*W_residuals * dlambda * psi;
    arma::mat W_dresiduals = -W % (dlambda * lambda_psi.t() +
      lambda_psi * dlambda.t());
    arma::mat dg2 = -2*W_dresiduals * lambda_psi;
    dglambda = dg1 + dg2;

    // dgpsi:
    W_dresiduals = -W % (lambda * dpsi * lambda.t());
    dgpsi = -2*lambda.t() * W_dresiduals * lambda;
    dgpsi.diag() *= 0.5;

    // dgtheta:
    dgtheta = 2*W % dtheta;
    dgtheta.diag() *= 0.5;

    dg(lambda_indices) += arma::vectorise(dglambda);
    // if(positive) {
    //   dg(T_indexes) += dgphi.elem(targetT_indexes);
    // } else {
    dg(psi_indices) += arma::vectorise(dgpsi);
    // }
    dg(theta_indices) += arma::vectorise(dgtheta);

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void H(arguments_optim& x) {

    // x.hessian.set_size(x.parameters.n_elem, x.parameters.n_elem); x.hessian.zeros();
    //
    // x.w = arma::vectorise(x.W);
    // // Rcpp::Rcout << "hlambda" << std::endl;
    // // Lambda
    // x.dlambda_dRhat_W = gLRhat(x.lambda, x.phi);
    // x.dlambda_dRhat_W.each_col() %= x.w;
    // x.lambda_phit_kron_Ip = arma::kron(x.lambda_phi.t(), x.Ip);
    // arma::mat g1 = 2*x.lambda_phit_kron_Ip * x.dlambda_dRhat_W;
    // arma::mat g2 = -2*arma::kron(x.phi, x.W_residuals);
    // arma::mat hlambda = g1 + g2;
    // x.hlambda = hlambda; //(x.lambda_indexes, x.lambda_indexes);
    // x.hessian(x.lambda_indexes, x.lambda_indexes) += x.hlambda(x.target_indexes, x.target_indexes);
    //
    // // Rcpp::Rcout << "hphi" << std::endl;
    // // Phi
    // arma::mat LL = 2*arma::kron(x.lambda, x.lambda).t();
    // x.dphi_dRhat_W = gPRhat(x.lambda, x.phi, x.indexes_q);
    // x.dphi_dRhat_W.each_col() %= x.w;
    // arma::mat hphi_temp = LL * x.dphi_dRhat_W;
    // hphi_temp.rows(x.indexes_diag_q) *= 0.5;
    // x.hphi = hphi_temp; //(x.phi_indexes, x.phi_indexes);
    // x.hessian(x.phi_indexes, x.phi_indexes) += x.hphi(x.targetphi_indexes, x.targetphi_indexes);
    //
    // // Rcpp::Rcout << "hpsi" << std::endl;
    // // Psi
    // x.dpsi_dRhat_W = gURhat(x.psi);
    // x.dpsi_dRhat_W.each_col() %= x.w;
    // arma::mat W2 = 2*x.W;
    // W2.diag() *= 0.5;
    // arma::vec w2 = arma::vectorise(W2);
    // arma::mat hpsi = arma::diagmat(w2);
    // x.hpsi = hpsi; //(x.psi_indexes, x.psi_indexes);
    // x.hessian(x.psi_indexes, x.psi_indexes) += x.hpsi(x.targetpsi_indexes, x.targetpsi_indexes);
    //
    // // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // // Lambda & Phi
    // g1 = 2*x.lambda_phit_kron_Ip * x.dphi_dRhat_W;
    // arma::mat g21 = -2*arma::kron(x.Iq, x.W_residuals_lambda);
    // arma::mat dtg21 = g21 * dxt(x.q, x.q);
    // g2 = g21 + dtg21;
    // g2.cols(x.indexes_diag_q) -= dtg21.cols(x.indexes_diag_q);
    // arma::mat dlambda_dphi_temp = g1 + g2;
    // x.dlambda_dphi = dlambda_dphi_temp; //(x.lambda_indexes, x.phi_indexes);
    // x.hessian(x.lambda_indexes, x.phi_indexes) += x.dlambda_dphi(x.target_indexes, x.targetphi_indexes);
    //
    // // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // // Lambda & Psi
    // arma::mat dlambda_dpsi_temp = 2*x.lambda_phit_kron_Ip * x.dpsi_dRhat_W;
    // x.dlambda_dpsi = dlambda_dpsi_temp; //(x.lambda_indexes, x.psi_indexes);
    // x.hessian(x.lambda_indexes, x.psi_indexes) += x.dlambda_dpsi(x.target_indexes, x.targetpsi_indexes);
    //
    // // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // // Phi & Psi
    // arma::mat dpsi_dphi_temp = x.dphi_dRhat_W;
    // dpsi_dphi_temp.rows(x.indexes_diag_p) *= 0.5;
    // dpsi_dphi_temp *= 2;
    // x.dpsi_dphi = dpsi_dphi_temp; //(x.psi_indexes, x.phi_indexes);
    // x.hessian(x.psi_indexes, x.phi_indexes) += x.dpsi_dphi(x.targetpsi_indexes, x.targetphi_indexes);
    //
    // x.hessian = arma::symmatu(x.hessian);

    /*
     * Join all the derivatives such that
     * hLambda           dlambda_dphi   dlambda_dpsi
     * dlambda_dphi.t()  hphi           dpsi_dphi.t()
     * dlambda_dpsi.t()   dpsi_dphi      hpsi
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

  return myestimator;

}
