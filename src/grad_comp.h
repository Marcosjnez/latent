/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/05/2025
 */

Rcpp::List grad_comp(arma::vec parameters,
                     Rcpp::List control_manifold,
                     Rcpp::List control_transform,
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

  x.nmanifolds = control_manifold.size();
  x.ntransforms = control_transform.size();
  x.nestimators = control_estimator.size();

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  std::vector<manifolds*> xmanifolds(x.nmanifolds);
  std::vector<transformations*> xtransforms(x.ntransforms);
  std::vector<estimators*> xestimators(x.nestimators);

  optim* algorithm = choose_optim(x, control_optimizer);
  x.parameters = parameters;

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i], xmanifolds[i]);
  }

  for(int i=0; i < x.ntransforms; ++i) {
    xtransforms[i] = choose_transform(control_transform[i], xtransforms[i]);
  }

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i], xestimators[i]);
  }

  /*
   * Computations
   */

  Rcpp::List computations;

  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  // Rf_error("67");
  final_transform->transform(x, xtransforms);
  // Rf_error("68");
  final_estimator->param(x, xestimators);
  // Rf_error("70");
  final_estimator->F(x, xestimators);
  // Rf_error("72");
  final_estimator->G(x, xestimators);
  // Rf_error("74");
  final_transform->jacobian(x, xtransforms);
  final_manifold->proj(x, xmanifolds);
  // Rf_error("78");

  Rcpp::List result;
  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;
  result["f"] = x.f;
  result["grad"] = x.grad;
  result["g"] = x.g;
  result["rg"] = x.rg;

  /*
   * DERIVATIVES
   */

  product_transform* T1;
  product_estimator* F1;
  // double eps = 1e-04;

  // Check the gradient of the transformed parameters:
  arma::vec transparameters = x.transparameters;
  arma::vec numg(transparameters.n_elem);
  for(int i = 0; i < transparameters.n_elem; ++i) {
    x.transparameters = transparameters; x.transparameters[i] += eps;
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f1 = x.f;
    x.transparameters = transparameters; x.transparameters[i] -= eps;
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f0 = x.f;
    numg[i] = (f1-f0) / (2*eps);
  }

  result["numg"] = numg;

  // Check the gradient:
  arma::vec numgrad(parameters.n_elem);
  for(int i = 0; i < parameters.n_elem; ++i) {
    x.parameters = parameters; x.parameters[i] += eps;
    T1->transform(x, xtransforms);
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f1 = x.f;
    x.parameters = parameters; x.parameters[i] -= eps;
    T1->transform(x, xtransforms);
    F1->param(x, xestimators);
    F1->F(x, xestimators);
    double f0 = x.f;
    numgrad[i] = (f1-f0) / (2*eps);
  }

  result["numgrad"] = numgrad;

  // final_estimator->E(x, xestimators);
  // final_estimator->M(x, xestimators);

  return result;

}
