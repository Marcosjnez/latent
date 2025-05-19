/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Confirmatory factor analysis (maximum-likelihood)
 */

class cfa_ml: public estimators {

public:
  arma::mat R, Rhat, lambda, phi, psi, residuals, model, Rhat_inv,
  Ri_res_Ri, lambda_phi, glambda, gphi, gpsi,
  dlambda, dglambda, dphi, dgphi, dpsi, dgpsi;
  double logdetR;
  int nparam, p;
  arma::uvec lambda_indexes, phi_indexes, psi_indexes,
  target_indexes, targetphi_indexes, targetpsi_indexes;

  void param() {

    lambda.elem(target_indexes) = parameters(lambda_indexes);
    phi.elem(targetphi_indexes) = parameters(phi_indexes);
    psi.elem(targetpsi_indexes) = parameters(psi_indexes);
    phi = arma::symmatl(phi);
    psi = arma::symmatl(psi);

    Rhat = lambda * phi * lambda.t() + psi;
    if(!Rhat.is_sympd()) {
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, Rhat);
      arma::vec d = arma::clamp(eigval, 0.1, eigval.max());
      Rhat = eigvec * arma::diagmat(d) * eigvec.t();
    }
    Rhat_inv = arma::inv_sympd(Rhat);

  }

  void F() {

    f = arma::log_det_sympd(Rhat) - logdetR +
      arma::accu(R % Rhat_inv) - p;

  }

  void G() {

    g.set_size(nparam); g.zeros();

    residuals = R - Rhat;
    Ri_res_Ri = 2*Rhat_inv * -residuals * Rhat_inv;
    lambda_phi = lambda * phi;

    glambda = Ri_res_Ri * lambda_phi;
    gphi = lambda.t() * Ri_res_Ri * lambda;
    gphi.diag() *= 0.5;

    arma::mat dlogdetRhatdU = 2*Rhat_inv;
    dlogdetRhatdU.diag() *= 0.5;
    arma::mat dRhat_invdU = 2*(-Rhat_inv * R * Rhat_inv);
    dRhat_invdU.diag() *= 0.5;
    gpsi = dRhat_invdU + dlogdetRhatdU;

    g(lambda_indexes) += glambda.elem(target_indexes);
    // if(positive) {
    //   g(T_indexes) += gphi.elem(targetT_indexes);
    // } else {
    g(phi_indexes) += gphi.elem(targetphi_indexes);
    // }
    g(psi_indexes) += gpsi.elem(targetpsi_indexes);

  }

  void dG() {

    // dX = arma::reshape(dparameters, q, q);
    dg.set_size(nparam); dg.zeros();

    dlambda.elem(target_indexes) = dparameters(lambda_indexes);
    dphi.elem(targetphi_indexes) = dparameters(phi_indexes);
    dpsi.elem(targetpsi_indexes) = dparameters(psi_indexes);
    dphi = arma::symmatl(dphi);
    dpsi = arma::symmatl(dpsi);

    // dglambda:
    arma::mat dRhat = dlambda * lambda_phi.t() + lambda_phi * dlambda.t();
    arma::mat dresiduals = -dRhat;
    arma::mat dRhat_inv = -Rhat_inv * -dresiduals * Rhat_inv;
    arma::mat dRi_res_Ri = 2*(dRhat_inv * -residuals * Rhat_inv +
      Rhat_inv * -residuals * dRhat_inv + Rhat_inv * -dresiduals * Rhat_inv);
    dglambda = Ri_res_Ri * dlambda * phi + dRi_res_Ri * lambda_phi;

    // dgphi:
    dRhat = lambda * dphi * lambda.t();
    dresiduals = -dRhat;
    dRhat_inv = -Rhat_inv * -dresiduals * Rhat_inv;
    dRi_res_Ri = 2*(dRhat_inv * -residuals * Rhat_inv + Rhat_inv * -residuals * dRhat_inv +
      Rhat_inv * -dresiduals * Rhat_inv);
    dgphi = lambda.t() * dRi_res_Ri * lambda;
    dgphi.diag() *= 0.5;

    // dgpsi:
    dRhat = dpsi;
    dRhat_inv = -Rhat_inv * dRhat * Rhat_inv;
    arma::mat ddlogdetRhat = 2*dRhat_inv.t();
    ddlogdetRhat.diag() *= 0.5;
    arma::mat ddRhat_inv = 2*(-dRhat_inv * R * Rhat_inv + -Rhat_inv * R * dRhat_inv);
    ddRhat_inv.diag() *= 0.5;
    dgpsi = ddlogdetRhat + ddRhat_inv;

    dg(lambda_indexes) += dglambda.elem(target_indexes);
    // if(positive) {
    //   dg(T_indexes) += dgphi.elem(targetT_indexes);
    // } else {
    dg(phi_indexes) += dgphi.elem(targetphi_indexes);
    // }
    dg(psi_indexes) += dgpsi.elem(targetpsi_indexes);

  }

  void E() {}

  void M() {}

  void H() {

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

  void outcomes() {

    // x.uniquenesses = x.R.diag() - arma::diagvec(x.Rhat);
    // x.Rhat.diag() = x.R.diag();

    vectors.resize(1);
    vectors[0] = psi.diag();

    matrices.resize(5);
    matrices[0] = lambda;
    matrices[1] = phi;
    matrices[2] = psi;
    matrices[3] = Rhat;
    matrices[4] = residuals;

    // list_matrices.resize(1);
    // list_matrices[0] = conditionals;

  };

};

cfa_ml* choose_cfa_ml(Rcpp::List estimator_setup) {

  cfa_ml* myestimator = new cfa_ml();

  arma::mat R = estimator_setup["R"];
  arma::mat lambda = estimator_setup["lambda"];
  arma::mat phi = estimator_setup["phi"];
  arma::mat psi = estimator_setup["psi"];
  double logdetR = estimator_setup["logdetR"];
  int nparam = estimator_setup["nparam"];
  int p = estimator_setup["p"];
  arma::uvec lambda_indexes = estimator_setup["lambda_indexes"];
  arma::uvec phi_indexes = estimator_setup["phi_indexes"];
  arma::uvec psi_indexes = estimator_setup["psi_indexes"];
  arma::uvec target_indexes = estimator_setup["target_indexes"];
  arma::uvec targetphi_indexes = estimator_setup["targetphi_indexes"];
  arma::uvec targetpsi_indexes = estimator_setup["targetpsi_indexes"];
  arma::uvec indices = estimator_setup["indices"];

  myestimator->R = R;
  myestimator->logdetR = logdetR;
  myestimator->lambda = lambda;
  myestimator->phi = phi;
  myestimator->psi = psi;
  myestimator->nparam = nparam;
  myestimator->p = p;
  myestimator->lambda_indexes = lambda_indexes;
  myestimator->phi_indexes = phi_indexes;
  myestimator->psi_indexes = psi_indexes;
  myestimator->target_indexes = target_indexes;
  myestimator->targetphi_indexes = targetphi_indexes;
  myestimator->targetpsi_indexes = targetpsi_indexes;
  myestimator->indices = indices;

  return myestimator;

}
