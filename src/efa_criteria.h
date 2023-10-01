/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 28/05/2022
 *
 */

// Criteria for factor extraction

class efa_criterion {

public:

  virtual void F(arguments_efa& x) = 0;

  virtual void G(arguments_efa& x) = 0;

  virtual void dG(arguments_efa& x) = 0;

  virtual void H(arguments_efa& x) = 0;

  virtual void H2(arguments_cfa& x) = 0;

  virtual void outcomes(arguments_efa& x) = 0;

};

/*
 * ULS
 */

class uls: public efa_criterion {

public:

  void F(arguments_efa& x) {

    // x.psi2 = 0.5*(x.lower + x.upper) + 0.5*abs(x.upper - x.lower) % sin(x.parameters);
    x.reduced_R = x.R - arma::diagmat(x.psi2);
    // x.reduced_R = x.R - arma::diagmat(x.parameters);
    eig_sym(x.eigval, x.eigvec, x.reduced_R);
    arma::vec e = x.eigval(arma::span(0, x.p - x.q - 1));

    x.f = 0.5*arma::accu(e % e);

  }

  void G(arguments_efa& x) {

    arma::vec e_values = x.eigval(arma::span(0, x.p - x.q - 1));
    arma::mat e_vectors = x.eigvec(arma::span::all, arma::span(0, x.p - x.q - 1));
    x.g_psi2 = -arma::diagvec(e_vectors * arma::diagmat(e_values) * e_vectors.t());

  }

  void dG(arguments_efa& x) {}

  void H(arguments_efa& x) {

    x.hessian.set_size(x.parameters.n_elem, x.parameters.n_elem); x.hessian.zeros();
    x.Rhat = x.lambda * x.phi * x.lambda.t() + x.psi;
    x.residuals = x.R - x.Rhat;
    x.lambda_phi = x.lambda * x.phi;
    x.W_residuals = x.W % x.residuals;
    x.W_residuals_lambda = x.W_residuals * x.lambda;

    x.w = arma::vectorise(x.W);
    // Rcpp::Rcout << "hlambda" << std::endl;
    // Lambda
    x.dlambda_dRhat_W = gLRhat(x.lambda, x.phi);
    x.dlambda_dRhat_W.each_col() %= x.w;
    x.lambda_phit_kron_Ip = arma::kron(x.lambda_phi.t(), x.Ip);
    arma::mat g1 = 2*x.lambda_phit_kron_Ip * x.dlambda_dRhat_W;
    arma::mat g2 = -2*arma::kron(x.phi, x.W_residuals);
    arma::mat hlambda = g1 + g2;
    x.hlambda = hlambda; //(x.lambda_indexes, x.lambda_indexes);
    x.hessian(x.lambda_indexes, x.lambda_indexes) += x.hlambda(x.target_indexes, x.target_indexes);

    // Rcpp::Rcout << "hphi" << std::endl;
    // Phi
    arma::mat LL = 2*arma::kron(x.lambda, x.lambda).t();
    x.dphi_dRhat_W = gPRhat(x.lambda, x.phi, x.indexes_q);
    x.dphi_dRhat_W.each_col() %= x.w;
    arma::mat hphi_temp = LL * x.dphi_dRhat_W;
    hphi_temp.rows(x.indexes_diag_q) *= 0.5;
    x.hphi = hphi_temp; //(x.phi_indexes, x.phi_indexes);
    x.hessian(x.phi_indexes, x.phi_indexes) += x.hphi(x.targetphi_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "hpsi" << std::endl;
    // Psi
    x.dpsi_dRhat_W = gURhat(x.psi);
    x.dpsi_dRhat_W.each_col() %= x.w;
    arma::mat W2 = 2*x.W;
    W2.diag() *= 0.5;
    arma::vec w2 = arma::vectorise(W2);
    arma::mat hpsi = arma::diagmat(w2);
    x.hpsi = hpsi; //(x.psi_indexes, x.psi_indexes);
    x.hessian(x.psi_indexes, x.psi_indexes) += x.hpsi(x.targetpsi_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // Lambda & Phi
    g1 = 2*x.lambda_phit_kron_Ip * x.dphi_dRhat_W;
    arma::mat g21 = -2*arma::kron(x.Iq, x.W_residuals_lambda);
    arma::mat dtg21 = g21 * dxt(x.q, x.q);
    g2 = g21 + dtg21;
    g2.cols(x.indexes_diag_q) -= dtg21.cols(x.indexes_diag_q);
    arma::mat dlambda_dphi_temp = g1 + g2;
    x.dlambda_dphi = dlambda_dphi_temp; //(x.lambda_indexes, x.phi_indexes);
    x.hessian(x.lambda_indexes, x.phi_indexes) += x.dlambda_dphi(x.target_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // Lambda & Psi
    arma::mat dlambda_dpsi_temp = 2*x.lambda_phit_kron_Ip * x.dpsi_dRhat_W;
    x.dlambda_dpsi = dlambda_dpsi_temp; //(x.lambda_indexes, x.psi_indexes);
    x.hessian(x.lambda_indexes, x.psi_indexes) += x.dlambda_dpsi(x.target_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // Phi & Psi
    arma::mat dpsi_dphi_temp = x.dphi_dRhat_W;
    dpsi_dphi_temp.rows(x.indexes_diag_p) *= 0.5;
    dpsi_dphi_temp *= 2;
    x.dpsi_dphi = dpsi_dphi_temp; //(x.psi_indexes, x.phi_indexes);
    x.hessian(x.psi_indexes, x.phi_indexes) += x.dpsi_dphi(x.targetpsi_indexes, x.targetphi_indexes);

    x.hessian = arma::symmatu(x.hessian);

  }

  void H2(arguments_cfa& x) {

    x.dLPU_dS.set_size(x.parameters.n_elem, x.p*x.p); x.dLPU_dS.zeros();

    // Rcpp::Rcout << "dlambda_dS" << std::endl;
    arma::mat g1 = -2*x.lambda_phit_kron_Ip;
    g1.each_row() %= x.w.t();
    arma::mat g2 = g1 * dxt(x.p, x.p);
    arma::mat g = g1 + g2;
    g.cols(x.indexes_diag_p) *= 0.5;
    x.dlambda_dS = g; //(x.lambda_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.lambda_indexes) += x.dlambda_dS.rows(x.target_indexes);

    // Rcpp::Rcout << "dphi_dS" << std::endl;
    g1 = -2*arma::kron(x.lambda.t(), x.lambda.t());
    g1.each_row() %= x.w.t();
    g2 = g1 * dxt(x.p, x.p);
    g = g1 + g2;
    g.cols(x.indexes_diag_p) *= 0.5;
    g.rows(x.indexes_diag_q) *= 0.5;
    x.dphi_dS = g; //(x.phi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.phi_indexes) += x.dphi_dS.rows(x.targetphi_indexes);

    // Rcpp::Rcout << "dpsi_dS" << std::endl;
    g = -2*arma::diagmat(x.w);
    g.cols(x.indexes_diag_p) *= 0.5;
    x.dpsi_dS = g; //(x.psi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.psi_indexes) += x.dpsi_dS.rows(x.targetpsi_indexes);

    // std::vector<arma::uvec> list_targetindices = {x.target_indexes,
    //                                               x.p*x.q + x.targetphi_indexes,
    //                                               x.p*x.q + 0.5*x.q*(x.q+1L) + x.targetpsi_indexes};
    // arma::uvec target_indices = list_to_vector(list_targetindices);
    // std::vector<arma::uvec> list_indices = {x.lambda_indexes, x.phi_indexes, x.psi_indexes};
    // arma::uvec indices = list_to_vector(list_indices);
    //
    // arma::mat dLPU_dS = arma::join_cols(x.dlambda_dS, x.dphi_dS, x.dpsi_dS);
    // x.dLPU_dS.rows(indices) += dLPU_dS.rows(target_indices);

  }

  void outcomes(arguments_efa& x) {

    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, x.reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, x.q-1));
    arma::vec eigenvalues = eigval2(arma::span(0, x.q-1));
    for(int i=0; i < x.q; ++i) {
      if(eigenvalues(i) < 0) eigenvalues(i) = 0;
    }
    arma::mat D = arma::diagmat(sqrt(eigenvalues));

    x.lambda = A * D;
    x.Rhat = x.lambda * x.lambda.t();
    x.uniquenesses = 1 - arma::diagvec(x.Rhat);// x.R.diag() - arma::diagvec(x.Rhat) FIX
    x.Rhat.diag() = x.R.diag();

  };

};

/*
 * ml
 */

class ml: public efa_criterion {

public:

  void F(arguments_efa& x) {

    // x.psi2 = 0.5*(x.lower + x.upper) + 0.5*abs(x.upper - x.lower) % sin(x.parameters);
    x.sqrt_psi = sqrt(x.psi2);
    // x.sqrt_psi = sqrt(x.parameters);
    arma::mat sc = arma::diagmat(1/x.sqrt_psi);
    x.reduced_R = sc * x.R * sc;
    eig_sym(x.eigval, x.eigvec, x.reduced_R);
    arma::vec e = x.eigval(arma::span(0, x.p - x.q - 1));

    // double objective = -arma::accu(log(e) + 1/e - 1);
    x.f = -(arma::accu(log(e) - e) + x.p - x.q);

  }

  void G(arguments_efa& x) {

    arma::mat A = x.eigvec(arma::span::all, arma::span(x.p-x.q, x.p-1));
    arma::vec eigenvalues = x.eigval(arma::span(x.p-x.q, x.p-1));
    // x.g = ((A % A) * (eigenvalues - 1) + 1 - arma::diagvec(x.R)/x.parameters)/x.parameters;
    x.g_psi2 = ((A % A) * (eigenvalues - 1) + 1 - arma::diagvec(x.R)/x.psi2)/x.psi2;

  }

  void dG(arguments_efa& x) {}

  void H(arguments_efa& x) {

    x.hessian.set_size(x.parameters.n_elem, x.parameters.n_elem); x.hessian.zeros();
    x.Rhat = x.lambda * x.phi * x.lambda.t() + x.psi;
    x.Rhat_inv = arma::inv_sympd(x.Rhat);
    x.residuals = x.R - x.Rhat;
    x.Ri_res_Ri = 2*x.Rhat_inv * -x.residuals * x.Rhat_inv;
    x.lambda_phi = x.lambda * x.phi;

    // Rcpp::Rcout << "hlambda" << std::endl;
    // Lambda
    arma::mat h1 = arma::kron(x.phi, x.Ri_res_Ri);
    arma::mat dRhat_dL = gLRhat(x.lambda, x.phi);
    arma::mat dRi_res_Ri_dRhat = 2*arma::kron(x.Rhat_inv, x.Rhat_inv) -
      arma::kron(x.Ri_res_Ri, x.Rhat_inv) - arma::kron(x.Rhat_inv, x.Ri_res_Ri);
    arma::mat dRi_res_Ri_dL = dRi_res_Ri_dRhat * dRhat_dL;
    arma::mat h2 = arma::kron(x.lambda_phi.t(), x.Ip) * dRi_res_Ri_dL;
    x.hlambda = h1 + h2;
    x.hessian(x.lambda_indexes, x.lambda_indexes) += x.hlambda(x.target_indexes, x.target_indexes);

    // Rcpp::Rcout << "hphi" << std::endl;
    // Phi
    arma::mat dRhat_dP = gPRhat(x.lambda, x.phi, x.indexes_diag_q);
    arma::mat dRi_res_Ri_dP = dRi_res_Ri_dRhat * dRhat_dP;
    x.hphi = arma::kron(x.lambda.t(), x.lambda.t()) * dRi_res_Ri_dP;
    x.hphi.rows(x.indexes_diag_q) *= 0.5;
    x.hessian(x.phi_indexes, x.phi_indexes) += x.hphi(x.targetphi_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "hpsi" << std::endl;
    // Psi
    arma::mat dRhat_dU = gURhat(x.psi);
    arma::mat dRi_res_Ri_dU = dRi_res_Ri_dRhat * dRhat_dU;
    x.hpsi = dRi_res_Ri_dU;
    x.hpsi.rows(x.indexes_diag_p) *= 0.5;
    x.hessian(x.psi_indexes, x.psi_indexes) += x.hpsi(x.targetpsi_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // Lambda & Phi
    h1 = arma::kron(x.lambda_phi.t(), x.Ip) * dRi_res_Ri_dP;
    arma::mat h21 = arma::kron(x.Iq, x.Ri_res_Ri * x.lambda);
    h2 = h21;
    h2 += h21 * dxt(x.q, x.q);
    h2.cols(x.indexes_diag_q) = h21.cols(x.indexes_diag_q);
    x.dlambda_dphi = h1 + h2; //(x.lambda_indexes, x.phi_indexes);
    x.hessian(x.lambda_indexes, x.phi_indexes) += x.dlambda_dphi(x.target_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // Lambda & Psi
    x.dlambda_dpsi = arma::kron(x.lambda_phi.t(), x.Ip) * dRi_res_Ri_dU;
    x.hessian(x.lambda_indexes, x.psi_indexes) += x.dlambda_dpsi(x.target_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // Phi & Psi
    arma::mat DXT = dxt(x.q, x.q);
    arma::mat dRhat_invdP = -arma::kron((x.lambda.t() * x.Rhat_inv.t()).t(), x.Rhat_inv * x.lambda);
    h1 = dRhat_invdP + dRhat_invdP * DXT;
    h2 = (arma::kron((x.R * x.Rhat_inv).t(), x.Ip) + arma::kron(x.Ip, x.Rhat_inv * x.R)) * h1;
    x.dpsi_dphi = 2*(h1 - h2);
    x.dpsi_dphi.cols(x.indexes_diag_q) *= 0.5;
    x.dpsi_dphi.rows(x.indexes_diag_p) *= 0.5;
    x.hessian(x.psi_indexes, x.phi_indexes) += x.dpsi_dphi(x.targetpsi_indexes, x.targetphi_indexes);

    x.hessian = arma::symmatu(x.hessian);

  }

  void H2(arguments_cfa& x) {

    x.dLPU_dS.set_size(x.parameters.n_elem, x.p*x.p); x.dLPU_dS.zeros();

    // Rcpp::Rcout << "dlambda_dS" << std::endl;
    arma::mat DXTS = dxt(x.p, x.p);
    arma::mat dRi_res_Ri_dS = -2*arma::kron(x.Rhat_inv, x.Rhat_inv);
    dRi_res_Ri_dS += dRi_res_Ri_dS * DXTS;
    arma::mat h = arma::kron(x.lambda_phi.t(), x.Ip) * dRi_res_Ri_dS;
    h.cols(x.indexes_diag_p) *= 0.5;
    x.dlambda_dS = h; //(x.lambda_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.lambda_indexes) += x.dlambda_dS.rows(x.target_indexes);

    // Rcpp::Rcout << "dphi_dS" << std::endl;
    h = arma::kron(x.lambda.t(), x.lambda.t()) * dRi_res_Ri_dS;
    h.cols(x.indexes_diag_p) *= 0.5;
    h.rows(x.indexes_diag_q) *= 0.5;
    x.dphi_dS = h; //(x.phi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.phi_indexes) += x.dphi_dS.rows(x.targetphi_indexes);

    // Rcpp::Rcout << "dpsi_dS" << std::endl;
    h = dRi_res_Ri_dS;
    h.cols(x.indexes_diag_p) *= 0.5;
    h.rows(x.indexes_diag_p) *= 0.5;
    x.dpsi_dS = h; //(x.psi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.psi_indexes) += x.dpsi_dS.rows(x.targetpsi_indexes);

  }

  void outcomes(arguments_efa& x) {

    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, x.reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, x.q-1));
    arma::vec eigenvalues = eigval2(arma::span(0, x.q-1)) - 1;
    for(int i=0; i < x.q; ++i) {
      if(eigenvalues[i] < 0) eigenvalues[i] = 0;
    }
    arma::mat D = diagmat(sqrt(eigenvalues));
    arma::mat w = A * D;

    x.lambda = diagmat(x.sqrt_psi) * w;
    x.Rhat = x.lambda * x.lambda.t();
    x.uniquenesses = 1 - diagvec(x.Rhat);// x.R.diag() - arma::diagvec(x.Rhat) FIX
    x.Rhat.diag() = x.R.diag();

  };

};

/*
 * DWLS
 */

class dwls: public efa_criterion {

public:

  void F(arguments_efa& x) {

    // W is a matrix with the variance of the polychoric correlations
    // Only the variance, not the covariances, are considered
    x.Rhat = x.lambda * x.lambda.t();// + arma::diagmat(x.uniquenesses);
    // x.Rhat.diag().ones();
    x.residuals = x.R - x.Rhat;
    x.f = 0.5*arma::accu(x.residuals % x.residuals % x.W);

  }

  void G(arguments_efa& x) {

    x.gL = -2*(x.residuals % x.W) * x.lambda; // * x.Phi;
    // arma::mat DW_res = x.residuals % x.DW;
    // x.gU = -arma::diagvec(DW_res);

  }

  void dG(arguments_efa& x) {}

  void H(arguments_efa& x) {

    x.hessian.set_size(x.parameters.n_elem, x.parameters.n_elem); x.hessian.zeros();
    x.Rhat = x.lambda * x.phi * x.lambda.t() + x.psi;
    x.residuals = x.R - x.Rhat;
    x.lambda_phi = x.lambda * x.phi;
    x.W_residuals = x.W % x.residuals;
    x.W_residuals_lambda = x.W_residuals * x.lambda;

    x.w = arma::vectorise(x.W);
    // Rcpp::Rcout << "hlambda" << std::endl;
    // Lambda
    x.dlambda_dRhat_W = gLRhat(x.lambda, x.phi);
    x.dlambda_dRhat_W.each_col() %= x.w;
    x.lambda_phit_kron_Ip = arma::kron(x.lambda_phi.t(), x.Ip);
    arma::mat g1 = 2*x.lambda_phit_kron_Ip * x.dlambda_dRhat_W;
    arma::mat g2 = -2*arma::kron(x.phi, x.W_residuals);
    arma::mat hlambda = g1 + g2;
    x.hlambda = hlambda; //(x.lambda_indexes, x.lambda_indexes);
    x.hessian(x.lambda_indexes, x.lambda_indexes) += x.hlambda(x.target_indexes, x.target_indexes);

    // Rcpp::Rcout << "hphi" << std::endl;
    // Phi
    arma::mat LL = 2*arma::kron(x.lambda, x.lambda).t();
    x.dphi_dRhat_W = gPRhat(x.lambda, x.phi, x.indexes_q);
    x.dphi_dRhat_W.each_col() %= x.w;
    arma::mat hphi_temp = LL * x.dphi_dRhat_W;
    hphi_temp.rows(x.indexes_diag_q) *= 0.5;
    x.hphi = hphi_temp; //(x.phi_indexes, x.phi_indexes);
    x.hessian(x.phi_indexes, x.phi_indexes) += x.hphi(x.targetphi_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "hpsi" << std::endl;
    // Psi
    x.dpsi_dRhat_W = gURhat(x.psi);
    x.dpsi_dRhat_W.each_col() %= x.w;
    arma::mat W2 = 2*x.W;
    W2.diag() *= 0.5;
    arma::vec w2 = arma::vectorise(W2);
    arma::mat hpsi = arma::diagmat(w2);
    x.hpsi = hpsi; //(x.psi_indexes, x.psi_indexes);
    x.hessian(x.psi_indexes, x.psi_indexes) += x.hpsi(x.targetpsi_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // Lambda & Phi
    g1 = 2*x.lambda_phit_kron_Ip * x.dphi_dRhat_W;
    arma::mat g21 = -2*arma::kron(x.Iq, x.W_residuals_lambda);
    arma::mat dtg21 = g21 * dxt(x.q, x.q);
    g2 = g21 + dtg21;
    g2.cols(x.indexes_diag_q) -= dtg21.cols(x.indexes_diag_q);
    arma::mat dlambda_dphi_temp = g1 + g2;
    x.dlambda_dphi = dlambda_dphi_temp; //(x.lambda_indexes, x.phi_indexes);
    x.hessian(x.lambda_indexes, x.phi_indexes) += x.dlambda_dphi(x.target_indexes, x.targetphi_indexes);

    // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // Lambda & Psi
    arma::mat dlambda_dpsi_temp = 2*x.lambda_phit_kron_Ip * x.dpsi_dRhat_W;
    x.dlambda_dpsi = dlambda_dpsi_temp; //(x.lambda_indexes, x.psi_indexes);
    x.hessian(x.lambda_indexes, x.psi_indexes) += x.dlambda_dpsi(x.target_indexes, x.targetpsi_indexes);

    // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // Phi & Psi
    arma::mat dpsi_dphi_temp = x.dphi_dRhat_W;
    dpsi_dphi_temp.rows(x.indexes_diag_p) *= 0.5;
    dpsi_dphi_temp *= 2;
    x.dpsi_dphi = dpsi_dphi_temp; //(x.psi_indexes, x.phi_indexes);
    x.hessian(x.psi_indexes, x.phi_indexes) += x.dpsi_dphi(x.targetpsi_indexes, x.targetphi_indexes);

    x.hessian = arma::symmatu(x.hessian);

  }

  void H2(arguments_cfa& x) {

    x.dLPU_dS.set_size(x.parameters.n_elem, x.p*x.p); x.dLPU_dS.zeros();

    // Rcpp::Rcout << "dlambda_dS" << std::endl;
    arma::mat g1 = -2*x.lambda_phit_kron_Ip;
    g1.each_row() %= x.w.t();
    arma::mat g2 = g1 * dxt(x.p, x.p);
    arma::mat g = g1 + g2;
    g.cols(x.indexes_diag_p) *= 0.5;
    x.dlambda_dS = g; //(x.lambda_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.lambda_indexes) += x.dlambda_dS.rows(x.target_indexes);

    // Rcpp::Rcout << "dphi_dS" << std::endl;
    g1 = -2*arma::kron(x.lambda.t(), x.lambda.t());
    g1.each_row() %= x.w.t();
    g2 = g1 * dxt(x.p, x.p);
    g = g1 + g2;
    g.cols(x.indexes_diag_p) *= 0.5;
    g.rows(x.indexes_diag_q) *= 0.5;
    x.dphi_dS = g; //(x.phi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.phi_indexes) += x.dphi_dS.rows(x.targetphi_indexes);

    // Rcpp::Rcout << "dpsi_dS" << std::endl;
    g = -2*arma::diagmat(x.w);
    g.cols(x.indexes_diag_p) *= 0.5;
    x.dpsi_dS = g; //(x.psi_indexes, x.S_indexes);
    x.dLPU_dS.rows(x.psi_indexes) += x.dpsi_dS.rows(x.targetpsi_indexes);

    // std::vector<arma::uvec> list_targetindices = {x.target_indexes,
    //                                               x.p*x.q + x.targetphi_indexes,
    //                                               x.p*x.q + 0.5*x.q*(x.q+1L) + x.targetpsi_indexes};
    // arma::uvec target_indices = list_to_vector(list_targetindices);
    // std::vector<arma::uvec> list_indices = {x.lambda_indexes, x.phi_indexes, x.psi_indexes};
    // arma::uvec indices = list_to_vector(list_indices);
    //
    // arma::mat dLPU_dS = arma::join_cols(x.dlambda_dS, x.dphi_dS, x.dpsi_dS);
    // x.dLPU_dS.rows(indices) += dLPU_dS.rows(target_indices);

  }

  void outcomes(arguments_efa& x) {

    x.Rhat = x.lambda * x.lambda.t();
    x.uniquenesses = 1 - arma::diagvec(x.Rhat); // x.R.diag() - arma::diagvec(x.Rhat) FIX
    x.Rhat.diag() = x.R.diag();

  };

};

// Choose the estimator:

efa_criterion* choose_efa_criterion(std::string estimator) {

  efa_criterion *criterion;

  if (estimator == "uls") {

    criterion = new uls();

  } else if(estimator == "ml") {

    criterion = new ml();

  } else if(estimator == "dwls") {

    criterion = new dwls();

  } else if(estimator == "gls") {

    Rcpp::stop("The 'gls' estimator is not available yet");

  } else if(estimator == "pa") {

  } else if(estimator == "minrank") {

    Rcpp::stop("The 'minrank' estimator is not available yet");

  } else {

    Rcpp::stop("Available estimators: \n uls, ml, dwls, gls, minrank, pa");

  }

  return criterion;

}
