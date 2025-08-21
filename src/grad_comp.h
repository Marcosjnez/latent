/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 21/08/2025
 */

Rcpp::List grad_comp(Rcpp::List control_manifold,
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

  Rcpp::List result;
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

  if(compute == "pre_opt") return result;
  optim* algorithm = choose_optim(x, control_optimizer);
  if(compute == "choose_optim") return result;

  std::vector<arma::vec> parameters_list = control_optimizer["parameters"];
  std::vector<arma::vec> transparameters_list = control_optimizer["transparameters"];
  arma::vec parameters = parameters_list[0];
  arma::vec transparameters = transparameters_list[0];
  x.parameters = parameters;
  x.transparameters = transparameters;

  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;

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
  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;
  if(compute == "mani_retr") return result;

  final_manifold->param(x, xmanifolds);
  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;
  if(compute == "mani") return result;

  final_transform->transform(x, xtransforms);
  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;
  if(compute == "trans") return result;

  final_estimator->param(x, xestimators);
  result["parameters"] = x.parameters;
  result["transparameters"] = x.transparameters;
  if(compute == "param") return result;

  final_estimator->F(x, xestimators);
  result["f"] = x.f;
  if(compute == "f") return result;

  final_estimator->G(x, xestimators);
  result["grad"] = x.grad;
  if(compute == "grad") return result;

  final_transform->update_grad(x, xtransforms);
  result["grad"] = x.grad;
  result["g"] = x.g;
  if(compute == "g") return result;

  final_manifold->proj(x, xmanifolds);
  result["rg"] = x.rg;
  if(compute == "rg") return result;

  final_estimator->H(x, xestimators);
  result["hess"] = x.hess;
  if(compute == "hess") return result;

  final_transform->update_hess(x, xtransforms);
  result["hess"] = x.hess;
  result["h"] = x.h;
  if(compute == "h") return result;

  final_transform->dconstraints(x, xtransforms);
  result["dconstr"] = x.mat_dconstraints;
  if(compute == "dconstr") return result;

  final_estimator->outcomes(x, xestimators);
  result["doubles"] = std::get<0>(x.outputs_estimator);
  result["vectors"] = std::get<1>(x.outputs_estimator);
  result["matrices"] = std::get<2>(x.outputs_estimator);
  result["cubes"] = std::get<3>(x.outputs_estimator);
  result["list_vectors"] = std::get<4>(x.outputs_estimator);
  result["list_matrices"] = std::get<5>(x.outputs_estimator);
  if(compute == "outcomes") return result;

  result["dg"] = x.dg;
  // result["J"] = xtransforms[0]->jacob;

  /*
   * NUMERICAL DERIVATIVES
   */

  product_transform* T1;
  product_estimator* F1;
  // double eps = 1e-04;

  // Check the gradient of the transformed parameters:
  // transparameters = x.transparameters;
  // arma::vec numgrad(transparameters.n_elem);
  // for(int i = 0; i < transparameters.n_elem; ++i) {
  //   x.transparameters = transparameters; x.transparameters[i] += eps;
  //   // T1->transform(x, xtransforms);
  //   F1->param(x, xestimators);
  //   F1->F(x, xestimators);
  //   double f1 = x.f;
  //   x.transparameters = transparameters; x.transparameters[i] -= eps;
  //   // T1->transform(x, xtransforms);
  //   F1->param(x, xestimators);
  //   F1->F(x, xestimators);
  //   double f0 = x.f;
  //   numgrad[i] = (f1-f0) / (2*eps);
  // }
  //
  // result["numgrad"] = numgrad;

  // Check the gradient:
  arma::vec numg(parameters.n_elem);
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
    numg[i] = (f1-f0) / (2*eps);
  }

  result["numg"] = numg;

  // final_estimator->E(x, xestimators);
  // final_estimator->M(x, xestimators);

  return result;

}
