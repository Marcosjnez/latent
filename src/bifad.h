Rcpp::List bifad(arma::mat R, int n_generals, int n_groups,
                 std::string projection,
                 Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                 double cutoff,
                 std::string normalization,
                 Rcpp::Nullable<int> nullable_nobs,
                 Rcpp::Nullable<Rcpp::List> first_efa,
                 Rcpp::Nullable<Rcpp::List> second_efa,
                 Rcpp::Nullable<Rcpp::List> nullable_rot_control,
                 int random_starts, int cores) {

  Rcpp::List first, second;

  if(first_efa.isNotNull()) {
    first = first_efa;
  }
  if(second_efa.isNotNull()) {
    second = second_efa;
  }

  int nfactors = n_generals + n_groups;

  // Arguments to pass to first efa:
  arguments_efast x1;
  // Check inputs:
  pass_to_efast(first, x1);

  // First EFA (group factors):

  Rcpp::List first_order_efa = efast(R, n_groups, x1.cor, x1.estimator, x1.rotation,
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

  Rcpp::List efa_model = first_order_efa["efa"];
  Rcpp::List efa_model_rotation = first_order_efa["rotation"];

  arma::mat unrotated = efa_model["lambda"];
  arma::mat loadings_1 = efa_model_rotation["lambda"];
  arma::mat Phi_1 = efa_model_rotation["phi"];

  arma::mat add;
  Rcpp::List second_order_efa;
  if(n_generals > 1) {

    // Arguments to pass to second efa:
    arguments_efast x2;
    // Check inputs:
    pass_to_efast(second, x2);

    // Second EFA (general factors):

    second_order_efa = efast(R, n_generals, x2.cor, x2.estimator, x2.rotation,
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

    Rcpp::List rot = second_order_efa["rotation"];

    arma::mat loadings_2 = rot["lambda"];
    arma::mat Phi_2 = rot["phi"];
    add = get_target(loadings_2, Phi_2, cutoff);

  } else {

    add.set_size(R.n_rows, 1);
    add.ones();

  }

  // Rotation (partially specified target with augmented loading matrix):

  unrotated.insert_cols(0, n_generals);

  // Create the augmented target matrix:
  arma::mat Target = get_target(loadings_1, Phi_1, cutoff);
  Target.insert_cols(0, add);
  SEXP Target_ = Rcpp::wrap(Target);
  Rcpp::Nullable<arma::mat> nullable_Target = Target_;
  // arma::mat Weight(Target.n_rows, Target.n_cols, arma::fill::eye);
  // SEXP Weight_ = Rcpp::wrap(Weight);
  // Rcpp::Nullable<arma::mat> nullable_Weight = Weight_;

  Rcpp::List rot = rotate(unrotated, "target", projection,
                          {0}, {0}, {0}, 0,
                          nullable_Target, R_NilValue, //nullable_Weight,
                          R_NilValue, R_NilValue,
                          R_NilValue, R_NilValue,
                          nullable_oblq_factors,
                          normalization,
                          nullable_rot_control,
                          random_starts, cores);

  Rcpp::List result;
  result["first_order"] = first_order_efa;
  result["second_order"] = second_order_efa;
  result["bifactor"] = rot;
  result["Target"] = Target;

  Rcpp::List modelInfo;
  modelInfo["R"] = R;
  modelInfo["n_generals"] = n_generals;
  modelInfo["n_groups"] = n_groups;
  modelInfo["projection"] = projection;
  modelInfo["cutoff"] = cutoff;
  modelInfo["normalization"] = normalization;
  modelInfo["nullable_nobs"] = nullable_nobs;
  modelInfo["first_efa"] = first_efa;
  modelInfo["second_efa"] = second_efa;
  modelInfo["nullable_rot_control"] = nullable_rot_control;

  result.attr("class") = "BIFAD";
  return result;

}

