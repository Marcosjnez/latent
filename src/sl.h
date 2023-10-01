Rcpp::List sl(arma::mat X, int n_generals, int n_groups,
              std::string cor, std::string estimator,
              std::string missing,
              Rcpp::Nullable<int> nullable_nobs,
              Rcpp::Nullable<Rcpp::List> first_efa,
              Rcpp::Nullable<Rcpp::List> second_efa, int cores) {

  Rcpp::Timer timer;
  Rcpp::List result;

  arguments_cor xcor;
  xcor.X = X;
  xcor.cor = cor;
  xcor.estimator = estimator;
  xcor.p = X.n_cols;
  xcor.q = n_generals + n_groups;
  xcor.missing = missing;

  check_cor(xcor);
  Rcpp::List correlation_result = xcor.correlation_result;

  result["correlation"] = correlation_result;

  if(nullable_nobs.isNotNull()) {
    xcor.nobs = Rcpp::as<int>(nullable_nobs);
  }

  Rcpp::List first, second;

  if(first_efa.isNotNull()) {
    first = first_efa;
  }

  if(second_efa.isNotNull()) {
    second = second_efa;
  }

  int nfactors = xcor.q;
  int n_items = xcor.p;

  // Arguments to pass to first efa in SL:

  arguments_efast x1;

  // Check inputs:

  pass_to_efast(first, x1);

  // Arguments to pass to second efa in SL:

  arguments_efast x2;

  // Check inputs:

  pass_to_efast(second, x2);

  // First efa:

  Rcpp::List first_order_efa = efast(xcor.R, n_groups, x1.cor, x1.estimator, x1.rotation,
                                     x1.projection, "none", nullable_nobs,
                                     x1.nullable_Target, x1.nullable_Weight,
                                     x1.nullable_PhiTarget, x1.nullable_PhiWeight,
                                     x1.nullable_blocks,
                                     x1.nullable_block_weights,
                                     x1.nullable_oblq_factors,
                                     x1.gamma, x1.epsilon, x1.k, x1.w,
                                     x1.random_starts, x1.cores,
                                     x1.nullable_init,
                                     x1.nullable_efa_control, x1.nullable_rot_control);

  Rcpp::List efa_model_rotation = first_order_efa["rotation"];

  arma::mat loadings_1 = efa_model_rotation["lambda"];
  arma::mat Phi_1 = efa_model_rotation["phi"];
  std::vector<std::string> rotation_none(1); rotation_none[0] = "none";

  if ( n_generals == 1 ) {

    Rcpp::List efa_result = efast(Phi_1, n_generals, x2.cor, x2.estimator, rotation_none,
                                  "none", "none", nullable_nobs,
                                  x2.nullable_Target, x2.nullable_Weight,
                                  x2.nullable_PhiTarget, x2.nullable_PhiWeight,
                                  x2.nullable_blocks,
                                  x2.nullable_block_weights,
                                  x2.nullable_oblq_factors,
                                  x2.gamma, x2.epsilon, x2.k, x2.w,
                                  x2.random_starts, x2.cores,
                                  x2.nullable_init,
                                  x2.nullable_efa_control, x2.nullable_rot_control);

    Rcpp::List efa_result_efa = efa_result["efa"];
    arma::mat loadings_2 = efa_result_efa["lambda"];
    arma::vec uniquenesses_2 = efa_result_efa["uniquenesses"];

    arma::mat L = join_rows(loadings_2, diagmat(sqrt(uniquenesses_2)));
    arma::mat SL_loadings = loadings_1 * L;

    for (int j=0; j < SL_loadings.n_cols; ++j) {
      if (sum(SL_loadings.col(j)) < 0) {
        SL_loadings.col(j) *= -1;
      }
    }

    arma::mat Hierarchical_Phi(1, 1, arma::fill::eye);
    efa_result_efa["phi"] = Hierarchical_Phi;

    arma::mat Rhat = SL_loadings * SL_loadings.t();
    arma::vec uniquenesses = 1 - diagvec(Rhat);
    Rhat.diag().ones();

    result["lambda"] = SL_loadings;
    result["first_order_solution"] = first_order_efa;
    result["second_order_solution"] = efa_result_efa;
    result["uniquenesses"] = uniquenesses;
    result["Rhat"] = Rhat;

  } else {

    Rcpp::List efa_result = efast(Phi_1, n_generals, x2.cor, x2.estimator, x2.rotation,
                                  x2.projection, "none", nullable_nobs,
                                  x2.nullable_Target, x2.nullable_Weight,
                                  x2.nullable_PhiTarget, x2.nullable_PhiWeight,
                                  x2.nullable_blocks,
                                  x2.nullable_block_weights,
                                  x2.nullable_oblq_factors,
                                  x2.gamma, x2.epsilon, x2.k, x2.w,
                                  x2.random_starts, x2.cores,
                                  x2.nullable_init,
                                  x2.nullable_efa_control, x2.nullable_rot_control);

    Rcpp::List efa_result_rotation = efa_result["rotation"];
    arma::mat loadings_2 = efa_result_rotation["lambda"];

    arma::vec uniquenesses_2 = efa_result_rotation["uniquenesses"];

    arma::mat Phi2 = efa_result_rotation["phi"];
    arma::mat SL_Phi(nfactors, nfactors, arma::fill::eye);
    SL_Phi(arma::span(0, n_generals-1), arma::span(0, n_generals-1)) = Phi2;
    // arma::mat sqrt_Phi2 = arma::sqrtmat_sympd(Phi2);

    arma::mat loadings_12 = loadings_1 * loadings_2;
    arma::mat sqrt_uniquenesses_2 = diagmat(sqrt(uniquenesses_2));
    arma::mat lu = loadings_1 * sqrt_uniquenesses_2;

    // arma::mat A = join_rows(loadings_12 * sqrt_Phi2, lu);
    arma::mat SL_loadings = join_rows(loadings_12, lu);
    arma::mat Rhat = SL_loadings * SL_Phi * SL_loadings.t();
    arma::vec uniquenesses = 1 - diagvec(Rhat);
    Rhat.diag().ones();

    result["lambda"] = SL_loadings;
    result["phi"] = SL_Phi;
    result["first_order_solution"] = first_order_efa;
    result["second_order_solution"] = efa_result;
    result["uniquenesses"] = uniquenesses;
    result["Rhat"] = Rhat;

  }

  Rcpp::List modelInfo;
  modelInfo["R"] = xcor.R;
  modelInfo["n_generals"] = n_generals;
  modelInfo["n_groups"] = n_groups;
  modelInfo["nullable_nobs"] = nullable_nobs;
  modelInfo["first_efa"] = first_efa;
  modelInfo["second_efa"] = second_efa;
  result["modelInfo"] = modelInfo;

  timer.step("elapsed");
  result["elapsed"] = timer;

  result.attr("class") = "SL";
  return result;

}

