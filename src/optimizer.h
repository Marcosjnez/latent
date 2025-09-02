/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 02/09/2025
 */

// By default, this optimizer performs minimization

Rcpp::List optimizer(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer) {

  /*
   * parameters: list of vectors of starting values
   */

  /* control_manifold: list of manifolds
   * Each list projects a set of parameters onto a manifold
   * NOTE: each list must contain nonoverlapping parameters
   * manifold_indices: list of indices that relate sets of parameters to manifolds
   */

  /* control_estimator: list of estimators
   * estimator_indices: list of indices that relate sets of parameters to estimators
   * estimator_target: list of indices that relate estimators to parameters
   */

  /* control_optimizer: list of parameters for the optimizer
   */

  /*
   * rstarts: number of independent repetitions of the optimization
   */

  /*
   * cores: number of cores
   */

  // Rf_error("STOP");

  Rcpp::List result;
  arguments_optim x;

  // Starting values for parameters:
  std::vector<arma::vec> parameters;
  if(control_optimizer.containsElementNamed("parameters")) {
    std::vector<arma::vec> params = control_optimizer["parameters"];
    parameters = params;
  }
  // Starting values for transparameters:
  std::vector<arma::vec> transparameters;
  if(control_optimizer.containsElementNamed("transparameters")) {
    std::vector<arma::vec> params = control_optimizer["transparameters"];
    transparameters = params;
  }
  // Starting values for posterior:
  std::vector<arma::mat> posterior;
  if(control_optimizer.containsElementNamed("posterior")) {
    std::vector<arma::mat> post = control_optimizer["posterior"];
    posterior = post;
  }

  int rstarts = control_optimizer["rstarts"];
  int cores = control_optimizer["cores"];

  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  optim* algorithm = choose_optim(x, control_optimizer);

  std::vector<std::vector<manifolds*>> xmanifolds(rstarts,
                                                  std::vector<manifolds*>(x.nmanifolds));
  std::vector<std::vector<transformations*>> xtransforms(rstarts,
                                                         std::vector<transformations*>(x.ntransforms));
  std::vector<std::vector<estimators*>> xestimators(rstarts,
                                                    std::vector<estimators*>(x.nestimators));

  for(int j=0; j < rstarts; ++j) {

    for(int i=0; i < x.nmanifolds; ++i) {
      xmanifolds[j][i] = choose_manifold(control_manifold[i]);
    }

    for(int i=0; i < x.ntransforms; ++i) {
      xtransforms[j][i] = choose_transform(control_transform[i]);
    }

    for(int i=0; i < x.nestimators; ++i) {
      xestimators[j][i] = choose_estimator(control_estimator[i]);

    }

  }

  /*
   * OPTIMIZATION
   */

  // Perform multiple random starts:

  arma::vec xf(rstarts); // vector containing the cost function of each random start
  std::vector<optim_result> x2(rstarts);
  std::vector<arguments_optim> args(rstarts);

  auto start = std::chrono::high_resolution_clock::now();

  #ifdef _OPENMP
    omp_set_num_threads(cores);
    #pragma omp parallel for
  #endif
  for (int i=0; i < rstarts; ++i) {

    args[i] = x;

    if(!parameters.empty()) {
      args[i].parameters = parameters[i];
    }
    if(!transparameters.empty()) {
      args[i].transparameters = transparameters[i];
    }
    if(!posterior.empty()) {
      args[i].posterior = posterior[i];
    }

    x2[i] = algorithm->optimize(args[i], xtransforms[i], xmanifolds[i], xestimators[i]);
    xf[i] = std::get<2>(x2[i]);

  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::nanoseconds elapsed_ns = end - start;
  double elapsed_sec = elapsed_ns.count() / 1e9;  // Convert to seconds

  // Pick up the best x.pick rstarts:
  if (x.pick > 0L) {

    Rcpp::List result(x.pick);
    arma::uvec indices = arma::sort_index(xf);

    for (int i = 0; i < x.pick; ++i) {
      const auto &t = x2[ indices(i) ];  // assume tuple<arma::vec, arma::vec, double, int, bool>

      Rcpp::List it = Rcpp::List::create(
        Rcpp::Named("parameters")       = std::get<0>(t),
        Rcpp::Named("transparameters")  = std::get<1>(t),
        Rcpp::Named("f")                = std::get<2>(t),
        Rcpp::Named("iterations")       = std::get<3>(t),
        Rcpp::Named("convergence")      = std::get<4>(t)
      );

      result[i] = it;  // store this iteration

    }

    return result;

  }

  // Choose the optimization process with the smallest objective value:
  arma::uword index_minimum = index_min(xf);
  optim_result x1 = x2[index_minimum];

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;
  // Rf_error("141");
  final_manifold->outcomes(args[index_minimum], xmanifolds[index_minimum]);
  final_transform->outcomes(args[index_minimum], xtransforms[index_minimum]);
  final_estimator->outcomes(args[index_minimum], xestimators[index_minimum]);
  // Rf_error("143");

  arma::vec opt_param = std::get<0>(x1);
  arma::vec opt_transparam = std::get<1>(x1);
  double f = std::get<2>(x1);
  int iterations = std::get<3>(x1);
  bool convergence = std::get<4>(x1);
  double ng = std::get<5>(x1);
  arma::mat rg = std::get<6>(x1);
  result["parameters"] = opt_param;
  result["transparameters"] = opt_transparam;
  result["f"] = f;
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  result["ng"] = ng;
  result["rg"] = rg;
  // result["doubles"] = args[index_minimum].doubles;
  // result["vectors"] = args[index_minimum].vectors;
  // result["matrices"] = args[index_minimum].matrices;
  // result["list_matrices"] = args[index_minimum].list_matrices;
  result["posterior"] = args[index_minimum].posterior;
  result["elapsed"] = elapsed_sec;

  Rcpp::List outputsMani;
  outputsMani["doubles"] = std::get<0>(args[index_minimum].outputs_manifold);
  outputsMani["vectors"] = std::get<1>(args[index_minimum].outputs_manifold);
  outputsMani["matrices"] = std::get<2>(args[index_minimum].outputs_manifold);
  outputsMani["cubes"] = std::get<3>(args[index_minimum].outputs_manifold);
  outputsMani["list_vectors"] = std::get<4>(args[index_minimum].outputs_manifold);
  outputsMani["list_matrices"] = std::get<5>(args[index_minimum].outputs_manifold);

  Rcpp::List outputsTrans;
  outputsTrans["doubles"] = std::get<0>(args[index_minimum].outputs_transform);
  outputsTrans["vectors"] = std::get<1>(args[index_minimum].outputs_transform);
  outputsTrans["matrices"] = std::get<2>(args[index_minimum].outputs_transform);
  outputsTrans["cubes"] = std::get<3>(args[index_minimum].outputs_transform);
  outputsTrans["list_vectors"] = std::get<4>(args[index_minimum].outputs_transform);
  outputsTrans["list_matrices"] = std::get<5>(args[index_minimum].outputs_transform);

  Rcpp::List outputsEst;
  outputsEst["doubles"] = std::get<0>(args[index_minimum].outputs_estimator);
  outputsEst["vectors"] = std::get<1>(args[index_minimum].outputs_estimator);
  outputsEst["matrices"] = std::get<2>(args[index_minimum].outputs_estimator);
  outputsEst["cubes"] = std::get<3>(args[index_minimum].outputs_estimator);
  outputsEst["list_vectors"] = std::get<4>(args[index_minimum].outputs_estimator);
  outputsEst["list_matrices"] = std::get<5>(args[index_minimum].outputs_estimator);

  Rcpp::List outputs;
  outputs["manifolds"] = outputsMani;
  outputs["transformations"] = outputsTrans;
  outputs["estimators"] = outputsEst;
  result["outputs"] = outputs;

  return result;

};
