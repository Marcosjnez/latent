/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

Rcpp::List get_df(Rcpp::S4 fit,
                  arma::vec dparameters,
                  double eps) {

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

  product_manifold final_manifold;
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

  arma::vec parameters = x.parameters;

  if(dparameters.n_elem != parameters.n_elem) {
    Rcpp::stop("dparameters must have one element per free parameter.");
  }
  if(!std::isfinite(eps) || eps <= 0.00) {
    Rcpp::stop("eps must be positive and finite.");
  }

  Rcpp::List computations;

  x.parameters = parameters + eps*dparameters;

  final_manifold.param(x, xmanifolds);
  final_manifold.retr(x, xmanifolds);
  final_manifold.param(x, xmanifolds);

  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  double f1 = x.f;

  x.parameters = parameters - eps*dparameters;

  final_manifold.param(x, xmanifolds);
  final_manifold.retr(x, xmanifolds);
  final_manifold.param(x, xmanifolds);

  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  double f2 = x.f;

  double df = (f1-f2)/(2*eps);

  result["df"] = df;

  return result;

}
