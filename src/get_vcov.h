/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

Rcpp::List get_vcov(Rcpp::S4 fit,
                    arma::mat vcov,
                    int cores) {

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

  product_transform final_transform;
  product_estimator final_estimator;

  std::vector<manifolds*> xmanifolds(x.nmanifolds);
  std::vector<transformations*> xtransforms(x.ntransforms);
  std::vector<estimators*> xestimators(x.nestimators);

  pointer_vector_guard<manifolds> manifold_guard(xmanifolds);
  pointer_vector_guard<transformations> transform_guard(xtransforms);
  pointer_vector_guard<estimators> estimator_guard(xestimators);

  std::unique_ptr<optim> algorithm(choose_optim(x, fit));

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

  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);

  x.v = vcov;
  final_transform.update_vcov(x, xtransforms);

  result["vcov"] = x.vcov;
  result["se"] = x.se;
  result["jacob"] = x.jacob;

  return result;

}
