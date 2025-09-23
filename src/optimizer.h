/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/09/2025
 */

// By default, this optimizer performs minimization of a loss function

Rcpp::List optimizer(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer) {

  Rcpp::List result;
  arguments_optim x;

  // Starting values for parameters:
  std::vector<arma::vec> parameters = control_optimizer["parameters"];

  // Starting values for transparameters:
  std::vector<arma::vec> transparameters = control_optimizer["transparameters"];

  // Number of different starting values and cores:
  int rstarts = control_optimizer["rstarts"];
  int cores = control_optimizer["cores"];

  // Size of model structures:
  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  // Check parameters for the optimizer:
  optim* algorithm = choose_optim(x, control_optimizer);

  // Initialize as many model structures as random starts:
  std::vector<std::vector<manifolds*>> xmanifolds(rstarts,
                                                  std::vector<manifolds*>(x.nmanifolds));
  std::vector<std::vector<transformations*>> xtransforms(rstarts,
                                                         std::vector<transformations*>(x.ntransforms));
  std::vector<std::vector<estimators*>> xestimators(rstarts,
                                                    std::vector<estimators*>(x.nestimators));

  // Check and setup each model structure:
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

  // Perform the optimization from multiple random starting values:

  // Vector containing the loss function of each random start:
  arma::vec xf(rstarts);
  // We will select the optimization procedure with minimum loss value
  std::vector<optim_result> x2(rstarts);
  std::vector<arguments_optim> args(rstarts);

  // Start the clock:
  auto start = std::chrono::high_resolution_clock::now();

  // Set up the parallelization:
  #ifdef _OPENMP
    omp_set_num_threads(cores);
    #pragma omp parallel for
  #endif
  for(int i=0; i < rstarts; ++i) {

    args[i] = x;

    args[i].parameters = parameters[i];
    args[i].transparameters = transparameters[i];

    x2[i] = algorithm->optimize(args[i], xtransforms[i], xmanifolds[i], xestimators[i]);
    xf[i] = std::get<2>(x2[i]);

  }

  // Stop the clock:
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::nanoseconds elapsed_ns = end - start;
  double elapsed_sec = elapsed_ns.count() / 1e9;  // Convert to seconds

  // Return the best x.pick rstarts:
  if (x.pick > 0L) {

    Rcpp::List result(x.pick);
    arma::uvec indices = arma::sort_index(xf);

    for (int i = 0; i < x.pick; ++i) {

      const auto &t = x2[ indices(i) ];

      Rcpp::List it = Rcpp::List::create(
        Rcpp::Named("parameters")       = std::get<0>(t),
        Rcpp::Named("transparameters")  = std::get<1>(t),
        Rcpp::Named("f")                = std::get<2>(t),
        Rcpp::Named("iterations")       = std::get<3>(t),
        Rcpp::Named("convergence")      = std::get<4>(t)
      );

      result[i] = it;  // store this iteration
      // tuple<arma::vec, arma::vec, double, int, bool>
    }

    result["elapsed"] = elapsed_sec;
    return result;

  }

  // Choose the optimization procedure with the smallest objective value:
  arma::uword index_minimum = index_min(xf);
  optim_result x1 = x2[index_minimum];

  // Compute some outcomes:
  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // Store the outputs generated from each function type and class:
  final_manifold->outcomes(args[index_minimum], xmanifolds[index_minimum]);
  final_transform->outcomes(args[index_minimum], xtransforms[index_minimum]);
  final_estimator->outcomes(args[index_minimum], xestimators[index_minimum]);

  // Extract the information from the fitted object:
  arma::vec opt_param = std::get<0>(x1); // Parameter estimates
  arma::vec opt_transparam = std::get<1>(x1); // Transformed parameters
  double f = std::get<2>(x1); // Minimum value of the loss function
  int iterations = std::get<3>(x1); // Number of iterations
  bool convergence = std::get<4>(x1); // Convergence
  double ng = std::get<5>(x1); // Norm of the gradient
  arma::mat rg = std::get<6>(x1); // Gradient vector
  result["parameters"] = opt_param;
  result["transparameters"] = opt_transparam;
  result["f"] = f;
  result["iterations"] = iterations;
  result["convergence"] = convergence;
  result["ng"] = ng;
  result["rg"] = rg;
  result["posterior"] = args[index_minimum].posterior;
  result["elapsed"] = elapsed_sec;

  // Store all the outputs of each manifold class:
  Rcpp::List outputsMani;
  outputsMani["doubles"] = std::get<0>(args[index_minimum].outputs_manifold);
  outputsMani["vectors"] = std::get<1>(args[index_minimum].outputs_manifold);
  outputsMani["matrices"] = std::get<2>(args[index_minimum].outputs_manifold);
  outputsMani["cubes"] = std::get<3>(args[index_minimum].outputs_manifold);
  outputsMani["list_vectors"] = std::get<4>(args[index_minimum].outputs_manifold);
  outputsMani["list_matrices"] = std::get<5>(args[index_minimum].outputs_manifold);

  // Store all the outputs of each transformation class:
  Rcpp::List outputsTrans;
  outputsTrans["doubles"] = std::get<0>(args[index_minimum].outputs_transform);
  outputsTrans["vectors"] = std::get<1>(args[index_minimum].outputs_transform);
  outputsTrans["matrices"] = std::get<2>(args[index_minimum].outputs_transform);
  outputsTrans["cubes"] = std::get<3>(args[index_minimum].outputs_transform);
  outputsTrans["list_vectors"] = std::get<4>(args[index_minimum].outputs_transform);
  outputsTrans["list_matrices"] = std::get<5>(args[index_minimum].outputs_transform);

  // Store all the outputs of each estimator class:
  Rcpp::List outputsEst;
  outputsEst["doubles"] = std::get<0>(args[index_minimum].outputs_estimator);
  outputsEst["vectors"] = std::get<1>(args[index_minimum].outputs_estimator);
  outputsEst["matrices"] = std::get<2>(args[index_minimum].outputs_estimator);
  outputsEst["cubes"] = std::get<3>(args[index_minimum].outputs_estimator);
  outputsEst["list_vectors"] = std::get<4>(args[index_minimum].outputs_estimator);
  outputsEst["list_matrices"] = std::get<5>(args[index_minimum].outputs_estimator);

  // Store all the outputs in a list:
  Rcpp::List outputs;
  outputs["manifolds"] = outputsMani;
  outputs["transformations"] = outputsTrans;
  outputs["estimators"] = outputsEst;
  result["outputs"] = outputs;

  return result;

};
