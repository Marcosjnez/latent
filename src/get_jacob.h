/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

Rcpp::List get_jacob(Rcpp::S4 fit) {

  Rcpp::List modelInfo = fit.slot("modelInfo");
  Rcpp::List control_transform = modelInfo["control_transform"];

  Rcpp::List result;
  arguments_optim x;
  x.ntransforms = control_transform.size();

  product_transform final_transform;
  std::vector<transformations*> xtransforms(x.ntransforms);
  pointer_vector_guard<transformations> transform_guard(xtransforms);

  // The Jacobian does not require a random directional derivative.
  std::unique_ptr<optim> algorithm(choose_optim(x, fit, false));

  for(int i=0; i < x.ntransforms; ++i) {
    xtransforms[i] = choose_transform(control_transform[i]);
  }

  /*
   * Computations
   */

  final_transform.transform(x, xtransforms);
  final_transform.jacobian(x, xtransforms);

  result["jacob"] = x.jacob;

  return result;

}
