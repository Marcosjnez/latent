/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Lower triangular DWLS
 */

class dwls_lt: public estimators {

public:

  arma::mat R, lambda, phi, psi, model, residuals, reduced_R, eigvec, Rhat,
  lambda_phi, W_residuals, W_residuals_lambda, W, dlambda_dRhat_W,
  lambda_phit_kron_Ip, Ip, Iq, hlambda, hphi, hpsi, dlambda_dS,
  dphi_dS, dpsi_dS, dphi_dRhat_W, dlambda_dphi, dlambda_dpsi, dpsi_dphi, dpsi_dRhat_W,
  Rhat_inv, Ri_res_Ri, glambda;
  arma::vec vector, sqrt_vector, eigval, w;
  int p, q;
  arma::uvec vector_indices, lambda_indices, targetlambda_indices, indices_q, indices_diag_q,
  indices_diag_p, phi_indices, targetphi_indices, psi_indices, targetpsi_indices;

  void param() {

    lambda.elem(targetlambda_indices) = parameters(lambda_indices);

  }

  void F() {

    // W is a matrix with the variance of the polychoric correlations
    // Only the variance, not the covariances, are considered
    Rhat = lambda * lambda.t();// + arma::diagmat(uniquenesses);
    // Rhat.diag().ones();
    residuals = R - Rhat;
    f = 0.5*arma::accu(residuals % residuals % W);

  }

  void G() {

    glambda = -2*(residuals % W) * lambda; // * Phi;
    // arma::mat DW_res = residuals % DW;
    // gU = -arma::diagvec(DW_res);
    arma::uvec alt_lower_indices = arma::trimatl_ind(arma::size(glambda), 0);
    g = glambda.elem(alt_lower_indices);

  }

  void dG() {

    // Rcpp::stop("dG not available");
    dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H() {

    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

    Rhat = lambda * phi * lambda.t() + psi;
    residuals = R - Rhat;
    lambda_phi = lambda * phi;
    W_residuals = W % residuals;
    W_residuals_lambda = W_residuals * lambda;

    w = arma::vectorise(W);
    // Rcpp::Rcout << "hlambda" << std::endl;
    // Lambda
    dlambda_dRhat_W = gLRhat(lambda, phi);
    dlambda_dRhat_W.each_col() %= w;
    lambda_phit_kron_Ip = arma::kron(lambda_phi.t(), Ip);
    arma::mat g1 = 2*lambda_phit_kron_Ip * dlambda_dRhat_W;
    arma::mat g2 = -2*arma::kron(phi, W_residuals);
    arma::mat hlambda = g1 + g2;
    hlambda = hlambda; //(lambda_indices, lambda_indices);
    hessian(lambda_indices, lambda_indices) += hlambda(targetlambda_indices,
            targetlambda_indices);

    // Rcpp::Rcout << "hphi" << std::endl;
    // Phi
    arma::mat LL = 2*arma::kron(lambda, lambda).t();
    dphi_dRhat_W = gPRhat(lambda, phi);
    dphi_dRhat_W.each_col() %= w;
    arma::mat hphi_temp = LL * dphi_dRhat_W;
    hphi_temp.rows(indices_diag_q) *= 0.5;
    hphi = hphi_temp; //(phi_indices, phi_indices);
    hessian(phi_indices, phi_indices) += hphi(targetphi_indices, targetphi_indices);

    // Rcpp::Rcout << "hpsi" << std::endl;
    // Psi
    dpsi_dRhat_W = gURhat(p);
    dpsi_dRhat_W.each_col() %= w;
    arma::mat W2 = 2*W;
    W2.diag() *= 0.5;
    arma::vec w2 = arma::vectorise(W2);
    arma::mat hpsi = arma::diagmat(w2);
    hpsi = hpsi; //(psi_indices, psi_indices);
    hessian(psi_indices, psi_indices) += hpsi(targetpsi_indices, targetpsi_indices);

    // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // Lambda & Phi
    g1 = 2*lambda_phit_kron_Ip * dphi_dRhat_W;
    arma::mat g21 = -2*arma::kron(Iq, W_residuals_lambda);
    arma::mat dtg21 = g21 * dxt(q, q);
    g2 = g21 + dtg21;
    g2.cols(indices_diag_q) -= dtg21.cols(indices_diag_q);
    dlambda_dphi = g1 + g2;
    hessian(lambda_indices, phi_indices) += dlambda_dphi(targetlambda_indices, targetphi_indices);

    // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // Lambda & Psi
    arma::mat dlambda_dpsi_temp = 2*lambda_phit_kron_Ip * dpsi_dRhat_W;
    dlambda_dpsi = dlambda_dpsi_temp; //(lambda_indices, psi_indices);
    hessian(lambda_indices, psi_indices) += dlambda_dpsi(targetlambda_indices, targetpsi_indices);

    // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // Phi & Psi
    arma::mat dpsi_dphi_temp = dphi_dRhat_W;
    dpsi_dphi_temp.rows(indices_diag_p) *= 0.5;
    dpsi_dphi_temp *= 2;
    dpsi_dphi = dpsi_dphi_temp; //(psi_indices, phi_indices);
    hessian(psi_indices, phi_indices) += dpsi_dphi(targetpsi_indices, targetphi_indices);

    hessian = arma::symmatu(hessian);

  }

  void E() {}

  void M() {}

  void outcomes() {

    model = lambda * lambda.t();
    uniquenesses = R.diag() - arma::diagvec(model);
    model.diag() = R.diag();
    residuals = R-model;
    phi.set_size(q, q); phi.eye();
    psi.set_size(p, p); psi.diag() = uniquenesses;

    /*
     * Compute the derivatives of the parameters wrt the correlation matrix
     */

    dparam_dS.set_size(nhessian, nS); dparam_dS.zeros();

    // Rcpp::Rcout << "dlambda_dS" << std::endl;
    arma::mat g1 = -2*lambda_phit_kron_Ip;
    g1.each_row() %= w.t();
    arma::mat g2 = g1 * dxt(p, p);
    arma::mat g = g1 + g2;
    g.cols(indices_diag_p) *= 0.5;
    dlambda_dS = g; //(lambda_indices, S_indices);
    dparam_dS.rows(lambda_indices) += dlambda_dS.rows(targetlambda_indices);

    // Rcpp::Rcout << "dphi_dS" << std::endl;
    g1 = -2*arma::kron(lambda.t(), lambda.t());
    g1.each_row() %= w.t();
    g2 = g1 * dxt(p, p);
    g = g1 + g2;
    g.cols(indices_diag_p) *= 0.5;
    g.rows(indices_diag_q) *= 0.5;
    dphi_dS = g; //(phi_indices, S_indices);
    dparam_dS.rows(phi_indices) += dphi_dS.rows(targetphi_indices);

    // Rcpp::Rcout << "dpsi_dS" << std::endl;
    g = -2*arma::diagmat(w);
    g.cols(indices_diag_p) *= 0.5;
    dpsi_dS = g; //(psi_indices, S_indices);
    dparam_dS.rows(psi_indices) += dpsi_dS.rows(targetpsi_indices);

  };

};

dwls_lt* choose_dwls_lt(const Rcpp::List& estimator_setup) {

 dwls_lt* myestimator = new dwls_lt(); // Creating dwls_lt object dynamically

  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  int nhessian = estimator_setup["nhessian"];
  int nS = estimator_setup["nS"];
  arma::uvec vp = consecutive(0, p-1);
  arma::uvec indices_diag_p = p*vp + vp;
  arma::uvec vq = consecutive(0, q-1);
  arma::uvec indices_diag_q = q*vq + vq;
  arma::mat Ip(p, p, arma::fill::eye); // Define identity pxp
  arma::mat Iq(q, q, arma::fill::eye); // Define identity qxq

  // Provide these:
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::uvec lambda_indices = estimator_setup["lambda_indices"];
  arma::uvec targetlambda_indices = estimator_setup["targetlambda_indices"];
  arma::uvec phi_indices = estimator_setup["phi_indices"];
  arma::uvec targetphi_indices = estimator_setup["targetphi_indices"];
  arma::uvec psi_indices = estimator_setup["psi_indices"];
  arma::uvec targetpsi_indices = estimator_setup["targetpsi_indices"];
  arma::mat R = estimator_setup["R"];
  arma::mat W = estimator_setup["W"];
  arma::mat lambda(p, q, arma::fill::zeros);

  myestimator->R = R;
  myestimator->W = W;
  myestimator->lambda = lambda;
  myestimator->p = p;
  myestimator->q = q;
  myestimator->nhessian = nhessian;
  myestimator->nS = nS;
  myestimator->indices_diag_p = indices_diag_p;
  myestimator->indices_diag_q = indices_diag_q;
  myestimator->Ip = Ip;
  myestimator->Iq = Iq;

  myestimator->indices = indices;
  myestimator->lambda_indices = lambda_indices;
  myestimator->phi_indices = phi_indices;
  myestimator->psi_indices = psi_indices;
  myestimator->targetlambda_indices = targetlambda_indices;
  myestimator->targetphi_indices = targetphi_indices;
  myestimator->targetpsi_indices = targetpsi_indices;

  return myestimator;

}
