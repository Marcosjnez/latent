/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 31/08/2025
 */

/*
 * EFA ULS
 */

class uls_efa: public estimators {

public:

  arma::mat R, lambda, phi, psi, model, residuals, reduced_R, eigvec, Rhat,
  lambda_phi, W_residuals, W_residuals_lambda, W, dlambda_dRhat_W,
  lambda_phit_kron_Ip, Ip, Iq, hlambda, hphi, hpsi, dlambda_dS,
  dphi_dS, dpsi_dS, dphi_dRhat_W, dlambda_dphi, dlambda_dpsi, dpsi_dphi, dpsi_dRhat_W,
  Rhat_inv, Ri_res_Ri, glambda;
  arma::vec vector, sqrt_vector, eigval, w;
  // double f;
  int p, q;
  arma::uvec vector_indices, lambda_indices, targetlambda_indices, indices_q, indices_diag_q,
  indices_diag_p, phi_indices, targetphi_indices, psi_indices, targetpsi_indices;

  void param(arguments_optim& x) {

    vector = transparameters;

  }

  void F(arguments_optim& x) {

    reduced_R = R - arma::diagmat(vector);
    eig_sym(eigval, eigvec, reduced_R);
    arma::vec e = eigval(arma::span(0, p - q - 1));

    f = 0.5*arma::accu(e % e);

  }

  void G(arguments_optim& x) {

    grad.set_size(transparameters.n_elem);
    grad.zeros();

    arma::vec e_values = eigval(arma::span(0, p - q - 1));
    arma::mat e_vectors = eigvec(arma::span::all, arma::span(0, p - q - 1));
    grad = -arma::diagvec(e_vectors * arma::diagmat(e_values) * e_vectors.t());

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

    // arma::vec eigval;
    // arma::mat eigvec;
    // eig_sym(eigval, eigvec, reduced_R);
    //
    // arma::vec eigval2 = reverse(eigval);
    // arma::mat eigvec2 = reverse(eigvec, 1);
    //
    // arma::mat A = eigvec2(arma::span::all, arma::span(0, q-1));
    // arma::vec eigenvalues = eigval2(arma::span(0, q-1));
    // for(int i=0; i < q; ++i) {
    //   if(eigenvalues(i) < 0) eigenvalues(i) = 0;
    // }
    // arma::mat D = arma::diagmat(sqrt(eigenvalues));
    //
    // lambda = A * D;
    // model = lambda * lambda.t();
    // uniquenesses = R.diag() - arma::diagvec(model);
    // model.diag() = R.diag();
    // residuals = R-model;
    // phi.set_size(q, q); phi.eye();
    // psi.set_size(p, p); psi.diag() = uniquenesses;

    /*
     * Compute the modified hessian (modhessian)
     */

    // modhessian.set_size(nhessian, nhessian); modhessian.zeros();
    //
    // Rhat = lambda * phi * lambda.t() + psi;
    // residuals = R - Rhat;
    // lambda_phi = lambda * phi;
    // W.set_size(p, p); W.ones();
    // W_residuals = W % residuals;
    // W_residuals_lambda = W_residuals * lambda;
    //
    // w = arma::vectorise(W);
    // // Rcpp::Rcout << "hlambda" << std::endl;
    // // Lambda
    // dlambda_dRhat_W = gLRhat(lambda, phi);
    // dlambda_dRhat_W.each_col() %= w;
    // lambda_phit_kron_Ip = arma::kron(lambda_phi.t(), Ip);
    // arma::mat g1 = 2*lambda_phit_kron_Ip * dlambda_dRhat_W;
    // arma::mat g2 = -2*arma::kron(phi, W_residuals);
    // arma::mat hlambda = g1 + g2;
    // hlambda = hlambda; //(lambda_indices, lambda_indices);
    // modhessian(lambda_indices, lambda_indices) += hlambda(targetlambda_indices, targetlambda_indices);
    //
    // // Rcpp::Rcout << "hphi" << std::endl;
    // // Phi
    // arma::mat LL = 2*arma::kron(lambda, lambda).t();
    // dphi_dRhat_W = gPRhat(lambda, phi);
    // dphi_dRhat_W.each_col() %= w;
    // arma::mat hphi_temp = LL * dphi_dRhat_W;
    // hphi_temp.rows(indices_diag_q) *= 0.5;
    // hphi = hphi_temp; //(phi_indices, phi_indices);
    // modhessian(phi_indices, phi_indices) += hphi(targetphi_indices, targetphi_indices);
    //
    // // Rcpp::Rcout << "hpsi" << std::endl;
    // // Psi
    // dpsi_dRhat_W = gURhat(p);
    // dpsi_dRhat_W.each_col() %= w;
    // arma::mat W2 = 2*W;
    // W2.diag() *= 0.5;
    // arma::vec w2 = arma::vectorise(W2);
    // arma::mat hpsi = arma::diagmat(w2);
    // hpsi = hpsi; //(psi_indices, psi_indices);
    // modhessian(psi_indices, psi_indices) += hpsi(targetpsi_indices, targetpsi_indices);
    //
    // // Rcpp::Rcout << "dlambda_dphi" << std::endl;
    // // Lambda & Phi
    // g1 = 2*lambda_phit_kron_Ip * dphi_dRhat_W;
    // arma::mat g21 = -2*arma::kron(Iq, W_residuals_lambda);
    // arma::mat dtg21 = g21 * dxt(q, q);
    // g2 = g21 + dtg21;
    // g2.cols(indices_diag_q) -= dtg21.cols(indices_diag_q);
    // dlambda_dphi = g1 + g2;
    // modhessian(lambda_indices, phi_indices) += dlambda_dphi(targetlambda_indices, targetphi_indices);
    //
    // // Rcpp::Rcout << "dlambda_dpsi" << std::endl;
    // // Lambda & Psi
    // arma::mat dlambda_dpsi_temp = 2*lambda_phit_kron_Ip * dpsi_dRhat_W;
    // dlambda_dpsi = dlambda_dpsi_temp; //(lambda_indices, psi_indices);
    // modhessian(lambda_indices, psi_indices) += dlambda_dpsi(targetlambda_indices, targetpsi_indices);
    //
    // // Rcpp::Rcout << "dpsi_dphi" << std::endl;
    // // Phi & Psi
    // arma::mat dpsi_dphi_temp = dphi_dRhat_W;
    // dpsi_dphi_temp.rows(indices_diag_p) *= 0.5;
    // dpsi_dphi_temp *= 2;
    // dpsi_dphi = dpsi_dphi_temp; //(psi_indices, phi_indices);
    // modhessian(psi_indices, phi_indices) += dpsi_dphi(targetpsi_indices, targetphi_indices);
    //
    // modhessian = arma::symmatu(modhessian);
    //
    // /*
    //  * Compute the derivatives of the parameters wrt the correlation matrix
    //  */
    //
    // dparam_dS.set_size(nhessian, nS); dparam_dS.zeros();
    //
    // // Rcpp::Rcout << "dlambda_dS" << std::endl;
    // g1 = -2*lambda_phit_kron_Ip;
    // g1.each_row() %= w.t();
    // g2 = g1 * dxt(p, p);
    // arma::mat g = g1 + g2;
    // g.cols(indices_diag_p) *= 0.5;
    // dlambda_dS = g; //(lambda_indices, S_indices);
    // dparam_dS.rows(lambda_indices) += dlambda_dS.rows(targetlambda_indices);
    //
    // // Rcpp::Rcout << "dphi_dS" << std::endl;
    // g1 = -2*arma::kron(lambda.t(), lambda.t());
    // g1.each_row() %= w.t();
    // g2 = g1 * dxt(p, p);
    // g = g1 + g2;
    // g.cols(indices_diag_p) *= 0.5;
    // g.rows(indices_diag_q) *= 0.5;
    // dphi_dS = g; //(phi_indices, S_indices);
    // dparam_dS.rows(phi_indices) += dphi_dS.rows(targetphi_indices);
    //
    // // Rcpp::Rcout << "dpsi_dS" << std::endl;
    // g = -2*arma::diagmat(w);
    // g.cols(indices_diag_p) *= 0.5;
    // dpsi_dS = g; //(psi_indices, S_indices);
    // dparam_dS.rows(psi_indices) += dpsi_dS.rows(targetpsi_indices);

  };

};

uls_efa* choose_uls_efa(const Rcpp::List& estimator_setup) {

  uls_efa* myestimator = new uls_efa();

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
