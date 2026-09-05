/*
 * Author: Marcos Jiménez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

Rcpp::List get_dconstr(Rcpp::S4 fit) {

  Rcpp::List modelInfo = fit.slot("modelInfo");
  Rcpp::List control_manifold = modelInfo["control_manifold"];
  Rcpp::List control_transform = modelInfo["control_transform"];
  Rcpp::List control_estimator = modelInfo["control_estimator"];
  Rcpp::List control_optimizer = modelInfo["control_optimizer"];

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

  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);
  final_manifold.param(x, xmanifolds);

  x.dconstr.set_size(x.transparameters.n_elem, 0L);
  x.d2constraints.set_size(0L, 0L);
  final_manifold.dconstraints(x, xmanifolds);
  final_transform.dconstraints(x, xtransforms);

  const arma::uword expected_dimension =
    x.transparameters.n_elem*x.dconstr.n_cols;

  if(x.d2constraints.n_rows != expected_dimension ||
     x.d2constraints.n_cols != expected_dimension) {
    Rcpp::stop("The first- and second-constraint derivative objects are not aligned.");
  }

  Rcpp::List result;
  result["dconstr"] = x.dconstr;
  result["d2constr"] = x.d2constraints;

  (void)algorithm;

  return result;

}
