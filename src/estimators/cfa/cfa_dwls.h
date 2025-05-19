/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Confirmatory factor analysis (least-squares, weigthed and unweighted)
 */

class cfa_dwls: public estimators {

public:

  arma::mat R, Rhat, lambda, phi, psi, residuals, model, W,
  lambda_phi, W_residuals, glambda, W_residuals_lambda, gphi, gpsi,
  dlambda, dglambda, dphi, dgphi, dpsi, dgpsi;
  int nparam;
  arma::uvec lambda_indexes, phi_indexes, psi_indexes,
  target_indexes, targetphi_indexes, targetpsi_indexes;

  void param() {

    // Rprintf("%u", parameters.n_elem);
    // Rprintf("\n");
    //
    // for (size_t i = 0; i < parameters.n_elem; ++i) {
    //   Rprintf("parameter ");
    //   Rprintf("%u", i);
    //   Rprintf(" = ");
    //   Rprintf("%f", parameters[i]);
    //   Rprintf("\n");
    // }

    lambda.elem(target_indexes) = parameters(lambda_indexes);
    phi.elem(targetphi_indexes) = parameters(phi_indexes);
    psi.elem(targetpsi_indexes) = parameters(psi_indexes);
    phi = arma::symmatl(phi);
    psi = arma::symmatl(psi);

    Rhat = lambda * phi * lambda.t() + psi;
    residuals = R - Rhat;

  }

  void F() {

    f = 0.5*arma::accu(residuals % residuals % W);

  }

  void G() {

    g.set_size(nparam); g.zeros();

    lambda_phi = lambda * phi;
    W_residuals = W % residuals;
    glambda = -2* W_residuals * lambda_phi;
    W_residuals_lambda = W_residuals * lambda;
    gphi = -2*lambda.t() * W_residuals_lambda;
    gphi.diag() *= 0.5;
    gpsi = -2*W_residuals;
    gpsi.diag() *= 0.5;

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
    arma::mat dg1 = -2*W_residuals * dlambda * phi;
    arma::mat W_dresiduals = -W % (dlambda * lambda_phi.t() +
      lambda_phi * dlambda.t());
    arma::mat dg2 = -2*W_dresiduals * lambda_phi;
    dglambda = dg1 + dg2;

    // dgphi:
    W_dresiduals = -W % (lambda * dphi * lambda.t());
    dgphi = -2*lambda.t() * W_dresiduals * lambda;
    dgphi.diag() *= 0.5;

    // dgpsi:
    dgpsi = 2*W % dpsi;
    dgpsi.diag() *= 0.5;

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

cfa_dwls* choose_cfa_dwls(Rcpp::List estimator_setup) {

  cfa_dwls* myestimator = new cfa_dwls();

  arma::mat R = estimator_setup["R"];
  arma::mat lambda = estimator_setup["lambda"];
  arma::mat phi = estimator_setup["phi"];
  arma::mat psi = estimator_setup["psi"];
  arma::mat W = estimator_setup["W"];
  int nparam = estimator_setup["nparam"];
  arma::uvec lambda_indexes = estimator_setup["lambda_indexes"];
  arma::uvec phi_indexes = estimator_setup["phi_indexes"];
  arma::uvec psi_indexes = estimator_setup["psi_indexes"];
  arma::uvec target_indexes = estimator_setup["target_indexes"];
  arma::uvec targetphi_indexes = estimator_setup["targetphi_indexes"];
  arma::uvec targetpsi_indexes = estimator_setup["targetpsi_indexes"];
  arma::uvec indices = estimator_setup["indices"];

  myestimator->R = R;
  myestimator->lambda = lambda;
  myestimator->phi = phi;
  myestimator->psi = psi;
  myestimator->W = W;
  myestimator->nparam = nparam;
  myestimator->lambda_indexes = lambda_indexes;
  myestimator->phi_indexes = phi_indexes;
  myestimator->psi_indexes = psi_indexes;
  myestimator->target_indexes = target_indexes;
  myestimator->targetphi_indexes = targetphi_indexes;
  myestimator->targetpsi_indexes = targetpsi_indexes;
  myestimator->indices = indices;

  return myestimator;

}
