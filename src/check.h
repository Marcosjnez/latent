/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/02/2025
 */

Rcpp::List check(arma::vec parameters,
                 arma::vec dparameters,
                 Rcpp::List control_manifold,
                 Rcpp::List control_estimator,
                 Rcpp::List control_optimizer,
                 double eps) {

  /*
   * parameters: vector of parameters
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

  arguments_optim x;

  x.parameters = parameters;

  x.dparameters = dparameters;
  x.nmanifolds = control_manifold.size();
  x.nestimators = control_estimator.size();

  product_manifold* final_manifold;
  product_estimator* final_estimator;

  std::vector<manifolds*> xmanifolds(x.nmanifolds);
  std::vector<estimators*> xestimators(x.nestimators);

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i], xmanifolds[i]);
  }

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i], xestimators[i]);
  }

  optim* algorithm = choose_optim(x, control_optimizer);

  /*
   * Computations
   */

  Rcpp::List computations;

  Rcpp::Rcout << "0" << std::endl;
  final_estimator->param(x, xestimators);
  Rcpp::Rcout << "1" << std::endl;
  final_estimator->F(x, xestimators);
  Rcpp::Rcout << "2" << std::endl;
  final_estimator->G(x, xestimators);
  Rcpp::Rcout << "3" << std::endl;
  final_manifold->param(x, xmanifolds);
  Rcpp::Rcout << "4" << std::endl;
  final_manifold->retr(x, xmanifolds);
  Rcpp::Rcout << "5" << std::endl;
  final_manifold->proj(x, xmanifolds);
  Rcpp::Rcout << "6" << std::endl;
  // final_estimator->dG(x, xestimators);
  // Rcpp::Rcout << "7" << std::endl;
  // final_manifold->hess(x, xmanifolds);
  // Rcpp::Rcout << "8" << std::endl;
  // final_estimator->outcomes(x, xestimators);
  // Rcpp::Rcout << "9" << std::endl;
  // final_estimator->H(x, xestimators);
  // Rcpp::Rcout << "10" << std::endl;

  final_estimator->param(x, xestimators);
  final_estimator->outcomes(x, xestimators);

  computations["parameters"] = x.parameters;
  computations["f"] = x.f;
  computations["g"] = x.g;
  computations["rg"] = x.rg;
  computations["dg"] = x.dg;
  computations["dH"] = x.dH;
  computations["hessian"] = x.hessian;
  computations["modhessian"] = x.modhessian;
  computations["dparam_dS"] = x.dparam_dS;
  // computations["doubles"] = x.doubles;
  // computations["vectors"] = x.vectors;
  // computations["matrices"] = x.matrices;
  // computations["list_matrices"] = x.list_matrices;

  /*
   * DERIVATIVES
   */

  Rcpp::List derivatives;

  product_estimator* F1;
  // double eps = 1e-04;

  // Check the gradient
  arma::vec numgrad(parameters.n_elem);
  for(int i = 0; i < parameters.n_elem; ++i) {
    x.parameters = parameters; x.parameters[i] += eps;
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f1 = x.f;
    x.parameters = parameters; x.parameters[i] -= eps;
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f0 = x.f;
    numgrad[i] = (f1-f0) / (2*eps);
  }

  x.parameters = parameters;
  F1->param(x, xestimators);
  F1->F(x, xestimators);
  double f = x.f;
  F1->G(x, xestimators);
  arma::vec gradient = x.g;

  derivatives["f"] = f;
  derivatives["gradient"] = gradient;
  derivatives["numgrad"] = numgrad;

  // Check the differential of the gradient
  // x.parameters = parameters; x.parameters += eps*dparameters;
  // F1->param(x, xestimators);
  // F1->G(x, xestimators);
  // arma::vec g1 = x.g;
  // x.parameters = parameters; x.parameters -= eps*dparameters;
  // F1->param(x, xestimators);
  // F1->G(x, xestimators);
  // arma::vec g0 = x.g;
  // arma::vec dnumgrad = (g1-g0) / (2*eps);
  //
  // x.parameters = parameters;
  // F1->param(x, xestimators);
  // F1->G(x, xestimators);
  // F1->dG(x, xestimators);
  // arma::vec dgradient = x.dg;
  //

  // derivatives["dgradient"] = dgradient;
  // derivatives["dnumgrad"] = dnumgrad;


  // Get the Hessian:
  arma::mat numhess(parameters.n_elem, parameters.n_elem);
  for(int i = 0; i < parameters.n_elem; ++i) {
    x.parameters = parameters; x.parameters[i] += eps;
    F1->param(x, xestimators);
    F1->G(x, xestimators);
    arma::vec g1 = x.g;

    x.parameters = parameters; x.parameters[i] -= eps;
    F1->param(x, xestimators);
    F1->G(x, xestimators);
    arma::vec g0 = x.g;

    numhess.col(i) = (g1-g0) / (2*eps);
  }

  derivatives["numhess"] = numhess;

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("computations") = computations,
    Rcpp::Named("derivatives") = derivatives
  );

  return result;

}
