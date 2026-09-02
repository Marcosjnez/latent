/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/09/2026
 */

Rcpp::List get_f(Rcpp::List control_manifold,
                 Rcpp::List control_transform,
                 Rcpp::List control_estimator,
                 Rcpp::List control_optimizer) {

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

  Rcpp::List computations;

  // final_manifold->param(x, xmanifolds);
  // final_manifold->retr(x, xmanifolds);
  // final_manifold->param(x, xmanifolds);

  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);

  result["f"] = x.f;

  return result;

}
