/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 28/08/2025
 */

Rcpp::List vcov_all(Rcpp::List control_manifold,
                    Rcpp::List control_transform,
                    Rcpp::List control_estimator,
                    Rcpp::List control_optimizer,
                    arma::mat H) {

  /*
   * parameters: vector of parameters
   */

  /* control_manifold: list of manifolds
   * Each list projects a set of parameters onto a manifold
   * NOTE: each list must contain nonoverlapping parameters
   */

  /* control_estimator: list of estimators
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

  optim* algorithm = choose_optim(x, control_optimizer);

  std::vector<arma::vec> parameters_list = control_optimizer["parameters"];
  std::vector<arma::vec> transparameters_list = control_optimizer["transparameters"];
  arma::vec parameters = parameters_list[0];
  arma::vec transparameters = transparameters_list[0];
  x.parameters = parameters;
  x.transparameters = transparameters;

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

  Rcpp::List computations;

  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  final_estimator->G(x, xestimators);
  final_transform->update_grad(x, xtransforms);
  final_manifold->proj(x, xmanifolds);
  final_estimator->H(x, xestimators);
  final_transform->update_hess(x, xtransforms);
  if(!H.is_empty()) {
    x.h = H;
  }
  final_transform->update_vcov(x, xtransforms);
  result["vcov"] = x.vcov;
  result["se"] = x.se;
  result["h"] = x.h;

  return result;

}
