/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 17/09/2023
 *
 */

Rcpp::List cfa(arma::vec parameters,
               std::vector<arma::mat> X,
               arma::ivec nfactors,
               arma::ivec nobs,
               std::vector<arma::mat> lambda,
               std::vector<arma::mat> phi,
               std::vector<arma::mat> psi,
               std::vector<arma::uvec> lambda_indexes,
               std::vector<arma::uvec> phi_indexes,
               std::vector<arma::uvec> psi_indexes,
               std::vector<arma::uvec> target_indexes,
               std::vector<arma::uvec> targetphi_indexes,
               std::vector<arma::uvec> targetpsi_indexes,
               std::vector<arma::uvec> free_indices_phi,
               Rcpp::CharacterVector char_cor,
               Rcpp::CharacterVector char_estimator,
               Rcpp::CharacterVector char_projection,
               Rcpp::CharacterVector char_missing,
               int random_starts, int cores,
               Rcpp::Nullable<Rcpp::List> nullable_control) {


  Rcpp::Timer timer;

  std::vector<std::string> cor = Rcpp::as<std::vector<std::string>>(char_cor);
  std::vector<std::string> estimator = Rcpp::as<std::vector<std::string>>(char_estimator);
  std::vector<std::string> projection = Rcpp::as<std::vector<std::string>>(char_projection);
  std::vector<std::string> missing = Rcpp::as<std::vector<std::string>>(char_missing);

  int nblocks = X.size();
  std::vector<arguments_cor> xcor(nblocks);
  std::vector<arguments_cfa> xcfa(nblocks);

  for(int i=0; i < nblocks; ++i) {

    // cor structure:
    xcor[i].X = X[i];
    xcor[i].cor = cor[i];
    xcor[i].estimator = estimator[i];
    xcor[i].p = X[i].n_cols;
    xcor[i].q = nfactors[i];
    xcor[i].missing = missing[i];
    xcor[i].cores = cores;
    xcor[i].nobs = nobs[i];

    check_cor(xcor[i]);

    // cfa structure:
    xcfa[i].parameters = parameters;
    xcfa[i].estimator = xcor[i].estimator;
    xcfa[i].projection = projection[i];
    // xcfa.optim = optim; // FIX optim in check_cfa
    // if(i == 1L) Rcpp::stop("4");
    xcfa[i].X = xcor[i].X;
    xcfa[i].p = xcor[i].p;
    xcfa[i].q = xcor[i].q;
    xcfa[i].R = xcor[i].R;
    xcfa[i].cor = xcor[i].cor;
    xcfa[i].missing = xcor[i].missing;
    xcfa[i].cores = xcor[i].cores;
    xcfa[i].nobs = xcor[i].nobs;
    xcfa[i].random_starts = random_starts;
    xcfa[i].nullable_control = nullable_control;
    xcfa[i].W = xcor[i].W;
    xcfa[i].std_error = xcor[i].std_error;
    xcfa[i].lambda = lambda[i];
    xcfa[i].phi = phi[i];
    // Rcpp::Rcout << xcfa[i].phi << std::endl;
    xcfa[i].psi = psi[i];
    xcfa[i].lambda_indexes = lambda_indexes[i]-1;
    xcfa[i].phi_indexes = phi_indexes[i]-1;
    xcfa[i].psi_indexes = psi_indexes[i]-1;
    xcfa[i].target_indexes = target_indexes[i]-1;
    xcfa[i].targetphi_indexes = targetphi_indexes[i]-1;
    xcfa[i].targetpsi_indexes = targetpsi_indexes[i]-1;
    // Rcpp::stop("1");
    if(!free_indices_phi[i].is_empty()){
      xcfa[i].free_indices_phi = free_indices_phi[i]-1;
    }
    // Rcpp::Rcout << xcfa[i].free_indices_phi << std::endl;

    check_cfa(xcfa[i]);

    // Rcpp::Rcout << xcfa[i].phi << std::endl;

  }

  // Select the optimizer:
  arguments_optim opt;
  opt.nblocks = nblocks;
  opt.parameters = parameters;
  opt.nullable_control = nullable_control;
  check_opt(opt);
  cfa_optim* algorithm = choose_cfa_optim(opt.optim);

  algorithm->optim(opt, xcfa);
  cfa_criterion2* cfa_criterion = new ultimate_criterion();
  cfa_criterion->outcomes(opt, xcfa);

  // Standard errors:
  cfa_criterion->H(opt, xcfa);
  cfa_criterion->H2(opt, xcfa);
  arma::mat inv_hessian;
  if(!opt.hessian.is_sympd()) {
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, opt.hessian);
    arma::vec d = arma::clamp(eigval, 0.01, eigval.max());
    inv_hessian = eigvec * arma::diagmat(1/d) * eigvec.t();
  } else {
    inv_hessian = arma::inv_sympd(opt.hessian);
  }
  // arma::mat inv_hessian = arma::inv(opt.hessian);
  arma::uvec indexes = trimatl_ind(arma::size(xcfa[0].R), 0);
  arma::mat full_Sigma = asymptotic_normal(xcfa[0].R);
  arma::mat Sigma = full_Sigma(indexes, indexes);
  arma::mat dLPU_dS = opt.dLPU_dS[0].cols(indexes);
  for(int i=1; i < opt.nblocks; ++i) {
    arma::uvec indexes = trimatl_ind(arma::size(xcfa[i].R), 0);
    arma::mat full_Sigma = asymptotic_normal(xcfa[i].R);
    Sigma = diagConc(Sigma, full_Sigma(indexes, indexes));
    dLPU_dS = arma::join_rows(dLPU_dS, opt.dLPU_dS[i].cols(indexes));
  }
  arma::mat B = dLPU_dS * Sigma * dLPU_dS.t();
  arma::mat COV = inv_hessian * B * inv_hessian;
  double denominator = (opt.total_nobs-1L)/opt.nblocks;
  opt.se = sqrt(arma::diagvec(COV)/denominator);

  Rcpp::List correlation;
  correlation["correlation"] = opt.R;

  Rcpp::List cfa;
  cfa["f"] = opt.f;
  cfa["iterations"] = opt.iteration;
  cfa["convergence"] = opt.convergence;
  cfa["lambda"] = opt.lambda;
  cfa["phi"] = opt.phi;
  cfa["psi"] = opt.psi;
  cfa["Rhat"] = opt.Rhat;
  cfa["fs"] = opt.fs;
  cfa["parameters"] = opt.parameters;
  cfa["se"] = opt.se;
  // cfa["hessian"] = opt.hessian;
  // cfa["H0"] = opt.B;
  // cfa["dLPU_dS"] = opt.dLPU_dS;
  // cfa["Phi_Target"] = opt.Phi_Target;
  // cfa["oblq_indexes"] = opt.oblq_indexes;

  Rcpp::List modelInfo;
  modelInfo["cor"] = opt.cor;
  modelInfo["estimator"] = opt.estimator;
  modelInfo["projection"] = opt.projection;
  modelInfo["nvars"] = opt.p;
  modelInfo["nfactors"] = opt.q;
  modelInfo["nobs"] = opt.nobs;
  modelInfo["df"] = opt.df;
  modelInfo["df_null"] = opt.df_null;
  modelInfo["f_null"] = opt.f_null;
  modelInfo["lower"] = opt.lower;
  modelInfo["upper"] = opt.upper;

  Rcpp::List result;
  result["correlation"] = correlation;
  result["cfa"] = cfa;
  result["modelInfo"] = modelInfo;

  timer.step("elapsed");
  result["elapsed"] = timer;

  result.attr("class") = "cfa";
  return result;

}

// [[Rcpp::export]]
Rcpp::List cfa_test(arma::mat R,
                    arma::mat lambda,
                    arma::mat phi,
                    arma::mat psi,
                    arma::mat dlambda,
                    arma::mat dphi,
                    arma::mat dpsi,
                    arma::mat W,
                    std::string estimator,
                    std::string projection) {

  arguments_cfa xcfa;

  xcfa.parameters = {1.00};
  xcfa.estimator = estimator;
  xcfa.projection = projection;
  xcfa.R = R;
  xcfa.p = lambda.n_rows;
  xcfa.q = lambda.n_cols;
  xcfa.lambda = lambda;
  xcfa.phi = phi;
  xcfa.psi = psi;
  xcfa.lambda_indexes = {0};
  xcfa.phi_indexes = {0};
  xcfa.psi_indexes = {0};
  xcfa.target_indexes = {0};
  xcfa.targetphi_indexes = {0};
  xcfa.targetpsi_indexes = {0};
  check_cfa(xcfa);
  xcfa.W = W;
  xcfa.dlambda = dlambda;
  xcfa.dphi = dphi;
  xcfa.dpsi = dpsi;

  cfa_criterion* cfa_criterion = choose_cfa_criterion(estimator);
  cfa_manifold* cfa_manifold = choose_cfa_manifold(projection);
  cfa_criterion->F(xcfa);
  cfa_criterion->G(xcfa);
  cfa_criterion->dG(xcfa);
  cfa_manifold->grad(xcfa);
  cfa_manifold->dgrad(xcfa);
  cfa_criterion->H(xcfa);
  cfa_criterion->H2(xcfa);

  Rcpp::List cfa;
  cfa["f"] = xcfa.f;
  cfa["lambda"] = xcfa.lambda;
  cfa["phi"] = xcfa.phi;
  cfa["psi"] = xcfa.psi;
  cfa["glambda"] = xcfa.glambda;
  cfa["gphi"] = xcfa.gphi;
  cfa["gpsi"] = xcfa.gpsi;
  cfa["dglambda"] = xcfa.dglambda;
  cfa["dgphi"] = xcfa.dgphi;
  cfa["dgpsi"] = xcfa.dgpsi;
  cfa["hlambda"] = xcfa.hlambda;
  cfa["dlambda_dphi"] = xcfa.dlambda_dphi;
  cfa["dlambda_dpsi"] = xcfa.dlambda_dpsi;
  cfa["hphi"] = xcfa.hphi;
  cfa["dpsi_dphi"] = xcfa.dpsi_dphi;
  cfa["hpsi"] = xcfa.hpsi;
  cfa["dlambda_dS"] = xcfa.dlambda_dS;
  cfa["dphi_dS"] = xcfa.dphi_dS;
  cfa["dpsi_dS"] = xcfa.dpsi_dS;

  return cfa;

}
