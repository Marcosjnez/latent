/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

// Line-search algorithm satisfying the Wolfe conditions:

void wolfe(arguments_optim& x,
           std::vector<transformations*>& xtransforms,
           std::vector<manifolds*>& xmanifolds,
           std::vector<estimators*>& xestimators) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  x.ss = 1;
  double max_dir = arma::abs(x.dir).max();
  if(max_dir > 1.0) {
    x.ss /= max_dir;
  }

  double f0 = x.f;
  arma::vec parameters = x.parameters;
  x.inprod = arma::dot(x.dir, x.rg);

  x.step_iteration = 0L;
  double lower = 0.0, upper = arma::datum::inf;
  double best_step = 0.0, best_f = f0;
  const double roundoff = 8*std::numeric_limits<double>::epsilon()*
    std::max(1.0, std::abs(f0));

  do {

    ++x.step_iteration;
    x.parameters = parameters + x.ss*x.dir;
    final_manifold.param(x, xmanifolds);
    final_manifold.retr(x, xmanifolds);
    final_transform.transform(x, xtransforms);
    final_estimator.param(x, xestimators);
    final_estimator.F(x, xestimators);
    x.df = x.f - f0;

    const double required = x.c1*x.ss*x.inprod;
    const bool sufficient = x.df <= required;
    // When a decrease is smaller than roundoff, use the gradient to resolve a tie.
    const bool rounded = std::abs(required) <= roundoff && x.df <= roundoff;
    if(!std::isfinite(x.f) || (!sufficient && !rounded)) {
      upper = x.ss;
    } else {
      final_estimator.G(x, xestimators);
      final_transform.update_grad(x, xtransforms);
      final_manifold.param(x, xmanifolds);
      final_manifold.proj(x, xmanifolds);
      const double slope = arma::dot(x.dir, x.rg);
      if(!std::isfinite(slope)) {
        upper = x.ss;
      } else if(slope >= x.c2*x.inprod) {
        if(sufficient || arma::norm(x.rg, 2) < x.ng) return;
        upper = x.ss;
      } else {
        if(sufficient && x.f < best_f) {
          best_step = x.ss;
          best_f = x.f;
        }
        lower = x.ss;
      }
    }

    x.ss = std::isfinite(upper) ? 0.5*(lower + upper) : x.ss*x.ss_fac;
  } while(x.step_iteration < x.step_maxit && x.ss >= x.step_eps);

  // Keep the best sufficient decrease when the curvature search is exhausted.
  x.ss = best_step;
  x.parameters = parameters + best_step*x.dir;
  final_manifold.param(x, xmanifolds);
  if(best_step > 0.0) final_manifold.retr(x, xmanifolds);
  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  x.df = x.f - f0;

}
