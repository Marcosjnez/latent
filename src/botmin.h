Rcpp::List botmin(arma::mat R, int n_generals, int n_groups,
                  std::string estimator, std::string projection,
                  Rcpp::Nullable<int> nullable_nobs,
                  Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                  double cutoff, int random_starts, int cores,
                  Rcpp::Nullable<Rcpp::List> nullable_efa_control,
                  Rcpp::Nullable<Rcpp::List> nullable_rot_control) {

  if(cutoff < 0) {
    Rcpp::stop("The cutoff needs to be nonnegative.");
  }

  // Arguments to pass to the first efa:

  arguments_efast x1;
  x1.estimator = estimator;
  x1.projection = "oblq";
  x1.nullable_oblq_factors = R_NilValue;
  x1.nullable_efa_control = nullable_efa_control;
  x1.nullable_rot_control = nullable_rot_control;
  x1.rotation = {"oblimin"};
  // x1.normalization; x1.gamma

  // First-order efa:

  Rcpp::List first_order_efa = efast(R, n_groups, x1.cor, x1.estimator, x1.rotation,
                                     x1.projection, "none", x1.nullable_nobs,
                                     x1.nullable_Target, x1.nullable_Weight,
                                     x1.nullable_PhiTarget, x1.nullable_PhiWeight,
                                     x1.nullable_blocks,
                                     x1.nullable_block_weights,
                                     x1.nullable_oblq_factors,
                                     x1.gamma, x1.epsilon, x1.k, x1.w,
                                     x1.random_starts, x1.cores,
                                     x1.nullable_init,
                                     x1.nullable_efa_control, x1.nullable_rot_control);

  Rcpp::List first_rot = first_order_efa["rotation"];
  arma::mat L = first_rot["lambda"];
  arma::mat Phi = first_rot["phi"];

  // Final efa:

  // Arguments to pass to the final efa:

  arguments_efast x2;
  x2.estimator = estimator;
  x2.projection = projection;
  x2.nullable_oblq_factors = nullable_oblq_factors;
  x2.nullable_efa_control = nullable_efa_control;
  x2.nullable_rot_control = nullable_rot_control;
  x2.rotation = {"oblimin", "target"};

  // Build the target using the cutoff:
  arma::mat Target = get_target(L, Phi, cutoff);
  SEXP Target_ = Rcpp::wrap(Target);
  x2.nullable_Target = Target_;
  // Set up the block list:
  unsigned g = n_generals; unsigned s = n_groups; unsigned p = R.n_rows;
  arma::uvec blocks_vector1 = {p, p};
  std::vector<arma::uvec> blocks_list1 = vector_to_list3(blocks_vector1);
  arma::uvec blocks_vector2 = {g, s};
  std::vector<arma::uvec> blocks_list2 = vector_to_list4(blocks_vector2);
  std::vector<std::vector<arma::uvec>> blocks_list(2);
  blocks_list[0] = blocks_list1;
  blocks_list[1] = blocks_list2;

  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks_list_ = Rcpp::wrap(blocks_list);
  // SEXP blocks_list_ = Rcpp::wrap(blocks_list);
  x2.nullable_blocks = blocks_list_;

  // conversion from R to C++
  // Rcpp::as<T>();
  // conversion from C++ to R
  // Rcpp::wrap();

  int nfactors = n_generals + n_groups;
  Rcpp::List final_efa = efast(R, nfactors, x2.cor, x2.estimator, x2.rotation,
                               x2.projection, "none", x2.nullable_nobs,
                               x2.nullable_Target, x2.nullable_Weight,
                               x2.nullable_PhiTarget, x2.nullable_PhiWeight,
                               x2.nullable_blocks,
                               x2.nullable_block_weights,
                               x2.nullable_oblq_factors,
                               x2.gamma, x2.epsilon, x2.k, x2.w,
                               x2.random_starts, x2.cores,
                               x2.nullable_init,
                               x2.nullable_efa_control, x2.nullable_rot_control);

  Rcpp::List efa_result = final_efa["efa"];
  Rcpp::List rotation_result = final_efa["rotation"];
  Rcpp::List modelInfo = final_efa["modelInfo"];
  Rcpp::List results;
  results["efa"] = efa_result;
  results["bifactor"] = rotation_result;
  results["modelInfo"] = modelInfo;
  results["first_order_solution"] = first_order_efa;

  return results;

}

