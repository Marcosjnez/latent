/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

// By default, this optimizer performs minimization

Rcpp::List optimizer(Rcpp::List control_manifold,
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
  // Starting values for posterior:
  std::vector<arma::mat> posterior;
  if(control_optimizer.containsElementNamed("posterior")) {
    std::vector<arma::mat> post = control_optimizer["posterior"];
    posterior = post;
  }

  int rstarts = control_optimizer["rstarts"];
  int cores = control_optimizer["cores"];

  x.nmanifolds = control_manifold.size();
  x.nestimators = control_estimator.size();
  std::vector<std::vector<manifolds*>> xmanifolds(rstarts,
                                                  std::vector<manifolds*>(x.nmanifolds));
  std::vector<std::vector<estimators*>> xestimators(rstarts,
                                                    std::vector<estimators*>(x.nestimators));

  optim* algorithm = choose_optim(x, control_optimizer);

  for(int j=0; j < rstarts; ++j) {

    for(int i=0; i < x.nmanifolds; ++i) {
      xmanifolds[j][i] = choose_manifold(control_manifold[i], xmanifolds[j][i]);
    }

    for(int i=0; i < x.nestimators; ++i) {
      xestimators[j][i] = choose_estimator(control_estimator[i], xestimators[j][i]);

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
    if(!posterior.empty()) {
      args[i].posterior = posterior[i];
    }

    x2[i] = algorithm->optimize(args[i], xmanifolds[i], xestimators[i]);
    xf[i] = std::get<1>(x2[i]);

  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::nanoseconds elapsed_ns = end - start;
  double elapsed_sec = elapsed_ns.count() / 1e9;  // Convert to seconds

  // Choose the optimization process with the smallest objective value:
  arma::uword index_minimum = index_min(xf);
  optim_result x1 = x2[index_minimum];

  product_estimator* final_estimator;
  final_estimator->outcomes(args[index_minimum], xestimators[index_minimum]);

  arma::vec opt_param = std::get<0>(x1);
  double f = std::get<1>(x1);
  int iterations = std::get<2>(x1);
  bool convergence = std::get<3>(x1);
  double ng = std::get<4>(x1);
  arma::mat rg = std::get<5>(x1);
  result["parameters"] = opt_param;
  result["f"] = f;
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  result["ng"] = ng;
  result["rg"] = rg;
  result["doubles"] = args[index_minimum].doubles;
  result["vectors"] = args[index_minimum].vectors;
  result["matrices"] = args[index_minimum].matrices;
  result["list_matrices"] = args[index_minimum].list_matrices;
  result["posterior"] = args[index_minimum].posterior;
  result["elapsed"] = elapsed_sec;

  return result;

};
