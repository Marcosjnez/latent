/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 31/08/2025
 */

/*
 * EFA ML
 */

class ml_efa: public estimators {

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

  void param(arguments_optim& x) {

    vector = parameters;

  }

  void F(arguments_optim& x) {

    sqrt_vector = sqrt(vector);
    arma::mat sc = arma::diagmat(1/sqrt_vector);
    reduced_R = sc * R * sc;
    eig_sym(eigval, eigvec, reduced_R);
    arma::vec e = eigval(arma::span(0, p - q - 1));

    // double objective = -arma::accu(log(e) + 1/e - 1);
    f = -(arma::accu(log(e) - e) + p - q);
    // Rcpp::Rcout << f << std::endl;

  }

  void G(arguments_optim& x) {

    arma::mat A = eigvec(arma::span::all, arma::span(p-q, p-1));
    arma::vec eigenvalues = eigval(arma::span(p-q, p-1));
    g = ((A % A) * (eigenvalues - 1) + 1 - arma::diagvec(R)/vector)/vector;

  }

  void dG(arguments_optim& x) {

    // Rcpp::stop("dG not available");
    dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H(arguments_optim& x) {

    // Rcpp::stop("H not available");
    hess.set_size(parameters.n_elem, parameters.n_elem); hess.zeros();

  }

  void E(arguments_optim& x) {}

  void M(arguments_optim& x) {}

  void outcomes(arguments_optim& x) {

    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, q-1));
    arma::vec eigenvalues = eigval2(arma::span(0, q-1)) - 1;
    for(int i=0; i < q; ++i) {
      if(eigenvalues[i] < 0) eigenvalues[i] = 0;
    }
    arma::mat D = diagmat(sqrt(eigenvalues));
    arma::mat w = A * D;

    lambda = arma::diagmat(sqrt_vector) * w;
    model = lambda * lambda.t();
    uniquenesses = R.diag() - arma::diagvec(model);
    model.diag() = R.diag();
    residuals = R-model;
    phi.set_size(q, q); phi.eye();
    psi.set_size(p, p); psi.diag() = uniquenesses;

    /*
     * Compute the modified hessian (modhessian)
     */

    modhessian.set_size(nhessian, nhessian); modhessian.zeros();

    Rhat = lambda * phi * lambda.t() + psi;
    Rhat_inv = arma::inv_sympd(Rhat);
    residuals = R - Rhat;
    Ri_res_Ri = 2*Rhat_inv * -residuals * Rhat_inv;
    lambda_phi = lambda * phi;

    // Rcpp::Rcout << "hlambda" << std::endl;
    // Lambda
    arma::mat h1 = arma::kron(phi, Ri_res_Ri);
    arma::mat dRhat_dL = gLRhat(lambda, phi);
    arma::mat dRi_res_Ri_dRhat = 2*arma::kron(Rhat_inv, Rhat_inv) -
      arma::kron(Ri_res_Ri, Rhat_inv) - arma::kron(Rhat_inv, Ri_res_Ri);
    arma::mat dRi_res_Ri_dL = dRi_res_Ri_dRhat * dRhat_dL;
    arma::mat h2 = arma::kron(lambda_phi.t(), Ip) * dRi_res_Ri_dL;
    hlambda = h1 + h2;
    modhessian(lambda_indices, lambda_indices) += hlambda(targetlambda_indices,
               targetlambda_indices);

    // Rcpp::Rcout << "hphi" << std::endl;
    // Phi
    arma::mat dRhat_dP = gPRhat(lambda, phi);
    arma::mat dRi_res_Ri_dP = dRi_res_Ri_dRhat * dRhat_dP;
    hphi = arma::kron(lambda.t(), lambda.t()) * dRi_res_Ri_dP;
    hphi.rows(indices_diag_q) *= 0.5;
    modhessian(phi_indices, phi_indices) += hphi(targetphi_indices, targetphi_indices);

    // Rcpp::Rcout << "hpsi" << std::endl;
    // Psi
    arma::mat dRhat_dU = gURhat(p);
    arma::mat dRi_res_Ri_dU = dRi_res_Ri_dRhat * dRhat_dU;
    hpsi = dRi_res_Ri_dU;
    hpsi.rows(indices_diag_p) *= 0.5;
    modhessian(psi_indices, psi_indices) += hpsi(targetpsi_indices, targetpsi_indices);

    // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // Lambda & Phi
    h1 = arma::kron(lambda_phi.t(), Ip) * dRi_res_Ri_dP;
    arma::mat h21 = arma::kron(Iq, Ri_res_Ri * lambda);
    h2 = h21;
    h2 += h21 * dxt(q, q);
    h2.cols(indices_diag_q) = h21.cols(indices_diag_q);
    dlambda_dphi = h1 + h2; //(lambda_indices, phi_indices);
    modhessian(lambda_indices, phi_indices) += dlambda_dphi(targetlambda_indices, targetphi_indices);

    // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // Lambda & Psi
    dlambda_dpsi = arma::kron(lambda_phi.t(), Ip) * dRi_res_Ri_dU;
    modhessian(lambda_indices, psi_indices) += dlambda_dpsi(targetlambda_indices, targetpsi_indices);

    // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // Phi & Psi
    arma::mat DXT = dxt(q, q);
    arma::mat dRhat_invdP = -arma::kron((lambda.t() * Rhat_inv.t()).t(), Rhat_inv * lambda);
    h1 = dRhat_invdP + dRhat_invdP * DXT;
    h2 = (arma::kron((R * Rhat_inv).t(), Ip) + arma::kron(Ip, Rhat_inv * R)) * h1;
    dpsi_dphi = 2*(h1 - h2);
    dpsi_dphi.cols(indices_diag_q) *= 0.5;
    dpsi_dphi.rows(indices_diag_p) *= 0.5;
    modhessian(psi_indices, phi_indices) += dpsi_dphi(targetpsi_indices, targetphi_indices);

    modhessian = arma::symmatu(modhessian);

    /*
     * Compute the derivatives of the parameters wrt the correlation matrix
     */

    dparam_dS.set_size(nhessian, nS); dparam_dS.zeros();

    // Rcpp::Rcout << "dlambda_dS" << std::endl;
    arma::mat DXTS = dxt(p, p);
    arma::mat dRi_res_Ri_dS = -2*arma::kron(Rhat_inv, Rhat_inv);
    dRi_res_Ri_dS += dRi_res_Ri_dS * DXTS;
    arma::mat h = arma::kron(lambda_phi.t(), Ip) * dRi_res_Ri_dS;
    h.cols(indices_diag_p) *= 0.5;
    dlambda_dS = h; //(lambda_indices, S_indices);
    dparam_dS.rows(lambda_indices) += dlambda_dS.rows(targetlambda_indices);

    // Rcpp::Rcout << "dphi_dS" << std::endl;
    h = arma::kron(lambda.t(), lambda.t()) * dRi_res_Ri_dS;
    h.cols(indices_diag_p) *= 0.5;
    h.rows(indices_diag_q) *= 0.5;
    dphi_dS = h; //(phi_indices, S_indices);
    dparam_dS.rows(phi_indices) += dphi_dS.rows(targetphi_indices);

    // Rcpp::Rcout << "dpsi_dS" << std::endl;
    h = dRi_res_Ri_dS;
    h.cols(indices_diag_p) *= 0.5;
    h.rows(indices_diag_p) *= 0.5;
    dpsi_dS = h; //(psi_indices, S_indices);
    dparam_dS.rows(psi_indices) += dpsi_dS.rows(targetpsi_indices);

  };

};

ml_efa* choose_ml_efa(const Rcpp::List& estimator_setup) {

  ml_efa* myestimator = new ml_efa();

  int p = estimator_setup["p"];
  int q = estimator_setup["q"];
  int nhessian = estimator_setup["nhessian"];
  int nS = estimator_setup["nS"];
  arma::uvec vp = consecutive(0, p-1);
  arma::uvec indices_diag_p = p*vp + vp;
  arma::uvec vq = consecutive(0, q-1);
  arma::uvec indices_diag_q = q*vq + vq;
  arma::mat Ip(p, p, arma::fill::eye); // Define identity
  arma::mat Iq(q, q, arma::fill::eye); // Define identity

  // Provide these:
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::uvec lambda_indices = estimator_setup["lambda_indices"];
  arma::uvec targetlambda_indices = estimator_setup["targetlambda_indices"];
  arma::uvec phi_indices = estimator_setup["phi_indices"];
  arma::uvec targetphi_indices = estimator_setup["targetphi_indices"];
  arma::uvec psi_indices = estimator_setup["psi_indices"];
  arma::uvec targetpsi_indices = estimator_setup["targetpsi_indices"];
  arma::mat R = estimator_setup["R"];

  myestimator->R = R;
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
