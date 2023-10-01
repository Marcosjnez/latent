Rcpp::List GSLiD(arma::mat R, int n_generals, int n_groups,
                 std::string estimator, std::string projection,
                 Rcpp::Nullable<int> nullable_nobs,
                 Rcpp::Nullable<arma::mat> nullable_Target,
                 Rcpp::Nullable<arma::mat> nullable_PhiTarget,
                 Rcpp::Nullable<arma::mat> nullable_PhiWeight,
                 double w, int maxit, int random_starts, int cores,
                 Rcpp::Nullable<arma::vec> nullable_init,
                 Rcpp::Nullable<Rcpp::List> nullable_efa_control,
                 Rcpp::Nullable<Rcpp::List> nullable_rot_control,
                 Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks,
                 Rcpp::Nullable<arma::vec> nullable_block_weights,
                 Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                 double cutoff, bool verbose) {

  std::vector<std::string> rotation;

  arguments_efa xefa;
  xefa.estimator = estimator;
  xefa.R = R;
  xefa.p = R.n_cols;
  xefa.q = n_generals + n_groups;
  xefa.upper = arma::diagvec(xefa.R);
  xefa.nullable_efa_control = nullable_efa_control;
  xefa.nullable_init = nullable_init;
  // Check inputs for efa:
  check_efa(xefa);

  // Structure of rotation arguments:

  arguments_rotate x;
  x.p = R.n_rows, x.q = n_generals + n_groups;
  x.lambda.set_size(x.p, x.q);
  x.Phi.set_size(x.q, x.q); x.Phi.eye();
  x.w = w;
  x.rotations = rotation;
  x.projection = projection;
  x.nullable_Target = nullable_Target;
  x.nullable_Weight = R_NilValue;
  x.nullable_PhiTarget = nullable_PhiTarget;
  x.nullable_PhiWeight = nullable_PhiWeight;
  x.nullable_blocks = nullable_blocks;
  x.nullable_oblq_factors = nullable_oblq_factors;
  x.nullable_block_weights = nullable_block_weights;
  x.nullable_rot_control = nullable_rot_control;

  // Choose manifold and rotation criteria:

  if(x.projection == "poblq" || x.projection == "orth") {
    x.rotations = {"target"};
  } else if(x.projection == "oblq") {
    x.rotations = {"xtarget"};
  } else {
    Rcpp::stop("Unkown projection method");
  }

  // Model Info:

  int p = R.n_cols;
  int q = x.q;
  double df_null = p*(p-1)/2;
  double df = p*(p+1)/2 - (p*q + p - q*(q-1)/2);

  double f_null;
  if(estimator == "uls" || estimator == "pa") {
    f_null = 0.5*(arma::accu(R % R) - R.n_cols);
  } else if(estimator == "dwls") {
    f_null = 0.5*arma::accu(R % R % xefa.W);
  } else if(estimator == "ml") {
    f_null = -arma::log_det_sympd(R);
  } else if(estimator == "minrank") {
    f_null = 0;
  }

  Rcpp::List modelInfo;
  modelInfo["R"] = R;
  modelInfo["estimator"] = estimator;
  modelInfo["projection"] = projection;
  modelInfo["rotation"] = x.rotations;
  modelInfo["nvars"] = R.n_cols;
  modelInfo["nfactors"] = q;
  modelInfo["nobs"] = nullable_nobs;
  modelInfo["df"] = df;
  modelInfo["df_null"] = df_null;
  modelInfo["f_null"] = f_null;
  modelInfo["w"] = x.w;
  modelInfo["k"] = x.k;
  modelInfo["gamma"] = x.gamma;
  modelInfo["epsilon"] = x.epsilon;
  modelInfo["w"] = x.w;
  modelInfo["normalization"] = x.normalization;
  modelInfo["Target"] = x.nullable_Target;
  modelInfo["Weight"] = x.nullable_Weight;
  modelInfo["PhiTarget"] = x.nullable_PhiTarget;
  modelInfo["PhiWeight"] = x.nullable_PhiWeight;
  modelInfo["blocks"] = x.nullable_blocks;
  modelInfo["block_weights"] = x.nullable_block_weights;
  modelInfo["oblq_factors"] = x.nullable_oblq_factors;

  // Check inputs and compute constants for rotation criteria:

  Rcpp::Nullable<arma::mat> nullable_Weight = R_NilValue;

  check_rotate(x, random_starts, cores);

  rotation_manifold* manifold = choose_manifold(projection);
  rotation_criterion *criterion = choose_criterion(x.rotations, x.projection, x.cols_list);

  arma::mat new_Target = x.Target;

  // Unrotated efa:

  // Select one manifold:
  efa_manifold* efa_manifold = choose_efa_manifold(xefa.manifold);
  // Select the estimator:
  efa_criterion* efa_criterion = choose_efa_criterion(xefa.estimator);

  Rcpp::List efa_result = efa(xefa, efa_manifold, efa_criterion,
                              xefa.random_starts, xefa.cores);

  xefa.heywood = efa_result["heywood"];
  arma::mat unrotated_loadings = efa_result["lambda"];
  arma::mat loadings = unrotated_loadings;

  Rcpp::List result, rotation_result;
  rotation_result["unrotated_loadings"] = unrotated_loadings;

  // Initialize stuff:

  x.lambda = unrotated_loadings; // x.lambda will always correspond to the unrotated loadings

  arma::vec congruence;
  arma::cube Targets(x.p, x.q, maxit, arma::fill::zeros);
  Targets.slice(0) = new_Target;
  arma::vec max_abs_diffs(maxit), min_congruences(maxit);
  int Target_discrepancies;
  bool Target_convergence = true;
  arma::mat old_Target;

  if (verbose) Rcpp::Rcout << "Rotating..." << std::endl;

  int i = 0;

  do{

    Rcpp::checkUserInterrupt();

    old_Target = new_Target;

    rotation_result = rotate_efa(x, manifold, criterion, random_starts,
                                 cores);

    arma::mat new_loadings = rotation_result["lambda"];

    // congruence = tucker_congruence(loadings, new_loadings);

    // min_congruences[i] = congruence.min();
    // max_abs_diffs[i] = arma::abs(loadings - new_loadings).max();

    loadings = new_loadings;
    arma::mat new_Phi = rotation_result["phi"];

    // update target
    update_target(n_generals, x.p, x.q, loadings, new_Phi, cutoff, new_Target);

    x.Target = new_Target;
    x.Weight = 1 - new_Target;
    x.Weight2 = x.Weight % x.Weight;

    Target_discrepancies = arma::accu(arma::abs(old_Target - new_Target));

    bool check = is_duplicate(Targets, new_Target, i);
    Targets.slice(i) = new_Target;

    ++i;

    if (verbose) Rcpp::Rcout << "\r" << "  Iteration " << i
                             << "  Target discrepancies = " << Target_discrepancies << "   \r";

    if(check) break;

  } while (i < maxit);

  if(i == maxit && Target_discrepancies != 0) {

    Rcpp::Rcout << "\n" << std::endl;
    Rcpp::warning("Maximum iteration reached without convergence");

    Target_convergence = false;

  } else if(Target_discrepancies != 0) {

    Rcpp::Rcout << "\n" << std::endl;
    Rcpp::warning("Recursive Target iterates. The last result of the iteration is returned");

    Target_convergence = false;

  }

  arma::mat L = loadings;
  arma::mat Phi = rotation_result["phi"];

  for (int j=0; j < x.q; ++j) {
    if (sum(L.col(j)) < 0) {
      L.col(j)   *= -1;
      Phi.col(j) *= -1;
      Phi.row(j) *= -1;
    }
  }

  arma::vec propVar = arma::diagvec(Phi * L.t() * L)/x.p;

  rotation_result["lambda"] = L;
  rotation_result["phi"] = Phi;
  arma::mat Rhat = L * Phi * L.t();
  rotation_result["uniquenesses"] = 1 - diagvec(Rhat);
  Rhat.diag().ones();
  rotation_result["Rhat"] = Rhat;
  rotation_result["Target"] = x.Target;
  rotation_result["Weights"] = x.Weight;
  rotation_result["Target_iterations"] = i;
  rotation_result["Target_convergence"] = Target_convergence;
  // rotation_result["min_congruences"] = min_congruences.head(i);
  // rotation_result["max_abs_diffs"] = max_abs_diffs.head(i);
  // Targets  = Targets(arma::span::all, arma::span::all, arma::span(0, i-1));
  // rotation_result["Targets"] = Targets;
  result["efa"] = efa_result;
  result["bifactor"] = rotation_result;
  result["modelInfo"] = modelInfo;

  return result;

}

