/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

// Gradient descent algorithm:

optim_result gd(arguments_optim x,
                std::vector<transformations*>& xtransforms,
                std::vector<manifolds*>& xmanifolds,
                std::vector<estimators*>& xestimators) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  std::unique_ptr<stepsize> step(choose_stepsize(x.step));

  // Ensure initial parameters are ok:
  final_manifold.param(x, xmanifolds);
  final_manifold.retr(x, xmanifolds);
  final_manifold.param(x, xmanifolds);
  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);

  // Parameterization

  final_estimator.F(x, xestimators);
  // update gradient
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);
  // Riemannian gradient
  final_manifold.proj(x, xmanifolds);

  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);

  do {

    step->update(x, xtransforms, xmanifolds, xestimators);

    // update gradient
    final_estimator.G(x, xestimators);
    final_transform.update_grad(x, xtransforms);
    // Riemannian gradient
    final_manifold.param(x, xmanifolds);
    final_manifold.proj(x, xmanifolds);

    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iterations;

    if(x.print) {
      Rcpp::Rcout << "Iteration = " << x.iterations << std::endl;
      Rcpp::Rcout << "f = " << x.f << std::endl;
      Rcpp::Rcout << "step size = " << x.ss << std::endl;
      Rcpp::Rcout << "step iters = " << x.step_iteration << std::endl;
      Rcpp::Rcout << "ng = " << x.ng << std::endl;
      Rcpp::Rcout << "dif = " << std::sqrt(x.df*x.df) << std::endl;
      Rcpp::Rcout << "" << std::endl;
    }
    if ((x.ng < x.eps) | (std::sqrt(x.df*x.df) < x.df_eps)) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.g, x.dir);

  return result;

}
