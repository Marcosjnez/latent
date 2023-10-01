/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 17/09/2023
 *
 */

typedef struct arguments_rotate{

  /*
   * Structure of arguments for rotation
   */

  int p, q, iteration = 0;
  std::vector<int> pi, qi;
  double w = 1, f, q2;
  arma::vec k = {0}, gamma = {0}, epsilon = {0.01}, clf_epsilon = {0.01};
  double ss = 1, inprod = 1, ng = 1;
  bool convergence = false;
  int maxit = 1e04;
  double eps = 1e-05;
  std::string optim = "newtonTR";
  std::string normalization = "none";

  Rcpp::Nullable<arma::mat> nullable_Target, nullable_Weight,
  nullable_PhiTarget, nullable_PhiWeight;
  Rcpp::Nullable<arma::uvec> nullable_oblq_factors, nullable_block_weights;
  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks; // NEW
  Rcpp::Nullable<Rcpp::List> nullable_rot_control;

  arma::mat lambda, T, L, Phi, Inv_T, dL, dP, Inv_T_dt, dT, g,
  gL, gP, dg, dgL, dgP, hL, hP, d_constr, d_constr_temp, rg, A,
  dH, f1, f2, L2, L2N, LoL2, ML2, IgCL2N, I_gamma_C, N, M, H, HL2, h, L2h,
  Target, Weight, Phi_Target, Phi_Weight, Weight2, Phi_Weight2, S,
  Ls, Lg, L2s, L2g, exp_aL2g, g_exp_aL2g, gL2g, gL2s, C, logC, logCN,
  gC1, gC, glogC, glogCN, gexplogCN, exp_lCN, gL1, gL2, I, I1, I2, Ng,
  dxtL, dxt_L2s, dmudL, dc2dL, dmudP, dc2dP, LtLxI, dxtP, expmu, expmmu, dir;

  arma::vec Lvec;
  arma::uvec lower, larger;
  std::vector<arma::uvec> loweri, largeri;

  arma::rowvec var, varq;
  std::vector<arma::rowvec> vari, varqi;
  arma::colvec mu, c2, dc2dmu;
  std::vector<arma::colvec> mui, c2i, expmui, expmmui, dc2dmui;
  arma::mat dvarqdL;
  std::vector<arma::mat> dvarqdLi, dmudLi, dc2dLi, dmudPi, dc2dPi, LtLxIi,
  dxtLi, dxtPi;
  double prodvar, prodvarq, prodc2;
  std::vector<double> prodvari, prodvarqi, prodc2i;

  arma::vec term;
  arma::uvec oblq_indexes, loblq_indexes, orth_indexes;
  std::vector<arma::uvec> list_oblq_indexes;
  std::vector<std::vector<arma::uvec>> blocks;
  std::vector<arma::uvec> rows_list, cols_list;
  arma::vec block_weights;
  int n_blocks = 1, n_rotations = 1, i, n_loads;

  std::vector<std::string> rotations = {"oblimin"};
  std::string projection = "oblq";
  std::vector<arma::mat> Li, Phii, Li2, Mi, Ni, HLi2, LoLi2, IgCL2Ni, V,
  f1i, f2i, Hi, Ii, HL2i, I_gamma_Ci;
  std::vector<arma::vec> termi;

  Rcpp::Nullable<arma::mat> nullable_indexes1, nullable_indexes2;
  arma::mat indexes1, indexes2;

} args;

typedef struct arguments_efast{

  int p, q;
  std::string cor = "pearson", estimator = "uls", projection = "oblq";
  std::vector<std::string> rotation = {"oblimin"};
  Rcpp::Nullable<int> nullable_nobs = R_NilValue;
  Rcpp::Nullable<arma::mat> nullable_Target = R_NilValue, nullable_Weight = R_NilValue,
    nullable_PhiTarget = R_NilValue, nullable_PhiWeight = R_NilValue;
  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks = R_NilValue;
  Rcpp::Nullable<arma::vec> nullable_block_weights = R_NilValue;
  Rcpp::Nullable<arma::uvec> nullable_oblq_factors = R_NilValue;
  std::string normalization = "none";
  double w = 1;
  arma::vec k = {0}, gamma = {0}, epsilon = {0.01}, clf_epsilon = {0.1};
  int random_starts = 10L, cores = 1L;
  Rcpp::Nullable<arma::vec> nullable_init = R_NilValue;
  Rcpp::Nullable<Rcpp::List> nullable_efa_control = R_NilValue,
    nullable_rot_control = R_NilValue;

} args_efast;

typedef struct arguments_efa{

  arma::mat init;
  arma::mat R, loadings, Rhat, residuals, Phi, W, X, gL, gU, g, smoothed;
  int q, p, nobs = 0;
  double f = 0;
  arma::mat lambda, phi, reduced_R, eigvec;
  arma::vec u, uniquenesses, eigval, parameters, sqrt_psi, psi2, g_psi2;
  std::string estimator = "uls";
  double efa_factr = 1e07;
  std::string optim = "gradient", cor = "pearson", missing = "pairwise.complete.cases",
    std_error = "normal";
  arma::vec lower = {0.005}, upper = {0.995};
  arma::mat rg, dir, dparameters, dH;
  int iteration, iterations = 0L, maxit = 1000L;
  bool convergence = false, heywood = false;
  std::string manifold = "box";
  int random_starts = 1L, cores = 1L;
  double ss = 1, inprod = 1, ng = 1, eps = 1e-05,
    c1 = 10e-04, c2 = 0.5, rho = 0.5;
  int M = 5L, armijo_maxit = 10L;
  std::string search = "back";
  std::string normalization = "none";
  int lambda_parameters;
  arma::uvec lower_tri_ind;
  Rcpp::Nullable<Rcpp::List> nullable_efa_control, nullable_first_efa,
  nullable_second_efa, nullable_init;
  Rcpp::List correlation_result;

  arma::mat hessian, psi, lambda_phi, W_residuals, W_residuals_lambda,
  dlambda_dRhat_W, lambda_phit_kron_Ip, hlambda, dphi_dRhat, dpsi_dRhat, hphi,
  hpsi, dphi_dRhat_W, Iq, dlambda_dphi, dlambda_dpsi, dpsi_dphi, Ip,
  dpsi_dRhat_W, Rhat_inv, Ri_res_Ri;
  arma::vec w;
  arma::uvec lambda_indexes, target_indexes, phi_indexes, targetphi_indexes,
  psi_indexes, targetpsi_indexes, indexes_q, indexes_diag_q, indexes_diag_p;


} args_efa;

typedef struct arguments_cor{

  int nobs, p, q, iteration = 0L, maxit = 1e04, cores = 1L;
  double f = 0.00, eps = 1e-05, ng = 1, ss = 1, inprod = 1, n_pairs;
  bool convergence = false;
  std::vector<size_t> s;
  std::vector<std::vector<double>> taus, mvphi;
  std::vector<std::vector<std::vector<int>>> n;
  arma::mat T, correlation, dT, dir;
  arma::mat g, rg;
  arma::mat dg, dH;
  arma::mat dcor, gcor, dgcor;

  std::string estimator = "uls", cor = "pearson", missing = "pairwise.complete.cases",
    std_error = "normal";
  arma::mat X, R, W;
  Rcpp::List correlation_result;

} args_cor;

typedef struct arguments_cfa{

  bool confirmatory = false;
  // CFA;
  arma::mat borrar;
  double f = 0.00, f_null = 0.00, logdetR;
  std::string estimator = "gls", projection = "id", missing = "pairwise.complete.cases",
    cor = "pearson";
  int nobs, p, q, n_lambda, n_phi, n_psi;
  arma::mat Ip, Iq; // Precompute

  arma::mat W_residuals, W_residuals_lambda, lambda_phit_kron_Ip,
  dlambda_dRhat_W, dphi_dRhat_W, dpsi_dRhat_W; // Repeated

  arma::mat X, W, R, Rhat, residuals, lambda_phi,
  lambda, phi, psi, dlambda, dphi, dpsi,
  glambda, gphi, gpsi, dglambda, dgphi, dgpsi,
  hlambda, hphi, hpsi, dlambda_dphi, dlambda_dpsi, dpsi_dphi, hessian,
  dlambda_dS, dphi_dS, dpsi_dS, dLPU_dS;

  arma::mat Rhat_inv, Ri_res_Ri;

  arma::vec parameters, transformed, dparameters, gradient, dgradient, uniquenesses, w;
  arma::uvec lambda_indexes, phi_indexes, psi_indexes, S_indexes,
  target_indexes, targetphi_indexes, targetpsi_indexes;
  arma::uvec indexes_diag_q, indexes_diag_p, indexes_diag_q2, indexes_p, indexes_q;
  arma::uvec target_positive;

  // Partially oblique in CFA
  arma::mat T, dT, gT, dgT, Phi_Target, A;
  arma::uvec oblq_indexes, free_indices_phi, T_indexes, targetT_indexes;
  bool positive = false;

  // Manifold stuff:
  arma::vec g, dg, rg, dH;

  // Optim stuff:
  arma::vec dir;
  double c1 = 10e-04, c2 = 0.5, rho = 0.5, eps = 1e-05, ng = 1, ss = 1, inprod = 1;
  int M = 5L, armijo_maxit = 10L, iteration = 0L, maxit = 1000L,
    random_starts= 1L, cores = 1L;
  std::string search = "back";
  bool convergence = false;

  // Checks:
  Rcpp::Nullable<Rcpp::List> nullable_control = R_NilValue;
  std::string optim = "L-BFGS", std_error = "normal";
  int df = 0L, df_null = 0L;

  // EFA
  arma::mat init;
  arma::mat loadings, Phi, gL, gU, smoothed;
  arma::mat reduced_R, eigvec;
  arma::vec u, eigval, sqrt_psi, psi2, g_psi2;
  double efa_factr = 1e07;
  arma::vec lower = {0.005}, upper = {0.995};
  int iterations = 0L;
  bool heywood = false;
  std::string normalization = "none";
  int lambda_parameters;
  arma::uvec lower_tri_ind;
  Rcpp::Nullable<Rcpp::List> nullable_efa_control, nullable_first_efa,
  nullable_second_efa, nullable_init;
  Rcpp::List correlation_result;

} args_cfa;

typedef struct arguments_optim{

  int nblocks;
  double f = 0.00, f_null = 0.00;
  // Optim stuff:
  double c1 = 10e-04, c2 = 0.5, rho = 0.5, eps = 1e-05, ng = 1, ss = 1, inprod = 1;
  int M = 5L, armijo_maxit = 10L, iteration = 0L, maxit = 10000L,
    random_starts = 1L, cores = 1L;
  std::string search = "back";
  bool convergence = false;
  arma::vec parameters, dparameters, gradient, dgradient, g, dg, rg, dH, dir;
  arma::mat hessian, B;
  std::vector<arma::mat> dLPU_dS;
  arma::vec se;

  // Checks:
  Rcpp::Nullable<Rcpp::List> nullable_control = R_NilValue;
  std::string optim = "L-BFGS", std_error = "normal";
  arma::uvec lower, upper;

  // Output:
  Rcpp::List lambda, phi, psi, Rhat, residuals, R;
  std::vector<double> fs;
  int df = 0L, df_null = 0L, total_nobs = 0L;
  std::vector<int> nobs, p, q;
  Rcpp::CharacterVector cor, estimator, projection;
  std::vector<arma::mat> Phi_Target;
  std::vector<arma::uvec> oblq_indexes;

} args_opt;

typedef std::tuple<arma::mat, arma::vec, arma::mat, double, int, bool> efa_NTR;
typedef std::tuple<arma::mat, arma::mat, arma::mat, double, int, bool> NTR;
typedef std::tuple<arma::mat, arma::mat, double, int, bool> cor_NTR;
typedef std::tuple<arma::vec, double, int, bool> cfa_NTR;
