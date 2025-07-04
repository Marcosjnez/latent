/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/06/2025
 */

Rcpp::List grad_comp(arma::vec parameters,
                     Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer,
                     std::string compute,
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

  Rcpp::List result;
  result["parameters"] = x.parameters;

  if(compute == "pre_setup") return result;

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i]);
  }

  if(compute == "choose_manifold") return result;

  for(int i=0; i < x.ntransforms; ++i) {
    xtransforms[i] = choose_transform(control_transform[i]);
  }

  if(compute == "choose_transform") return result;

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i]);
  }

  if(compute == "choose_estimator") return result;

  /*
   * Computations
   */

  Rcpp::List computations;

  if(compute == "setup") return result;
  final_manifold->param(x, xmanifolds);
  if(compute == "mani_param") return result;
  final_manifold->retr(x, xmanifolds);
  if(compute == "mani_retr") return result;
  final_manifold->param(x, xmanifolds);
  if(compute == "mani") return result;
  final_transform->transform(x, xtransforms);
  result["transparameters"] = x.transparameters;
  if(compute == "trans") return result;
  final_estimator->param(x, xestimators);
  if(compute == "param") return result;
  final_estimator->F(x, xestimators);
  result["f"] = x.f;
  if(compute == "f") return result;
  final_estimator->G(x, xestimators);
  result["grad"] = x.grad;
  if(compute == "grad") return result;
  final_transform->jacobian(x, xtransforms);
  // final_transform->outcomes(x, xtransforms);
  result["g"] = x.g;
  // result["matrices"] = xtransforms[x.ntransforms-1L]->matrices;
  if(compute == "g") return result;
  final_manifold->proj(x, xmanifolds);

  result["dg"] = x.dg;
  result["rg"] = x.rg;
  result["J"] = xtransforms[0]->jacob;
  result["h"] = x.h;

  /*
   * NUMERICAL DERIVATIVES
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
