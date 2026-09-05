/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

Rcpp::List get_dgrad(Rcpp::S4 fit,
                     arma::vec dparameters) {

  Rcpp::List modelInfo = fit.slot("modelInfo");
  Rcpp::List control_manifold = modelInfo["control_manifold"];
  Rcpp::List control_transform = modelInfo["control_transform"];
  Rcpp::List control_estimator = modelInfo["control_estimator"];
  Rcpp::List control_optimizer = modelInfo["control_optimizer"];

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

  optim* algorithm = choose_optim(x, control_optimizer);

  for(int i=0; i < x.nmanifolds; ++i) {
    xmanifolds[i] = choose_manifold(control_manifold[i]);
  }

  for(int i=0; i < x.ntransforms; ++i) {
    xtransforms[i] = choose_transform(control_transform[i]);
  }

  for(int i=0; i < x.nestimators; ++i) {
    xestimators[i] = choose_estimator(control_estimator[i]);
  }

  /*
   * Computations
   */

  Rcpp::List Optim = fit.slot("Optim");
  arma::vec parameters = Optim["parameters"];
  arma::vec transparameters = Optim["transparameters"];
  x.parameters = parameters;
  x.transparameters = transparameters;
  x.dparameters = dparameters;

  Rcpp::List computations;

  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);

  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  final_estimator->G(x, xestimators);
  final_transform->update_grad(x, xtransforms);
  final_transform->dtransform(x, xtransforms);
  final_estimator->dG(x, xestimators);
  final_transform->update_dgrad(x, xtransforms);

  result["f"] = x.f;
  result["g"] = x.g;
  result["dg"] = x.dg;
  result["grad"] = x.grad;
  result["dgrad"] = x.dgrad;

  return result;

}
