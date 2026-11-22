/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

// L-BFGS algorithm:

optim_result lbfgs(arguments_optim& x,
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
  x.iterations = 0L;

  // Parameterization
  // final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  // update the gradient
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);
  // Riemannian gradient
  // final_manifold.param(x, xmanifolds);
  final_manifold.proj(x, xmanifolds);
  x.dir = -x.rg;
  x.ng = arma::norm(x.rg, 2);
  x.convergence = false;

  // Keep only the recent parameter and gradient differences used by L-BFGS.
  const int memory = std::min(x.M, x.maxit);
  std::vector<arma::vec> s(memory), y(memory);
  std::vector<double> rho(memory), alpha(memory);
  int next = 0, count = 0;
  const double curvature_eps = std::sqrt(std::numeric_limits<double>::epsilon());

  while(std::isfinite(x.f) && std::isfinite(x.ng)) {

    if(x.rg.is_empty() || x.ng <= x.eps*std::sqrt(x.rg.n_elem)) {
      x.convergence = true;
      break;
    }
    if(x.iterations >= x.maxit) break;

    arma::vec old_parameters = x.parameters;
    arma::vec old_rg = x.rg;
    step->update(x, xtransforms, xmanifolds, xestimators);
    ++x.iterations;

    final_estimator.G(x, xestimators);
    final_transform.update_grad(x, xtransforms);
    final_manifold.param(x, xmanifolds);
    final_manifold.proj(x, xmanifolds);
    x.ng = arma::norm(x.rg, 2);

    if(x.print && x.rstarts == 1L && x.print_interval > 0 &&
       x.iterations % x.print_interval == 0) {
      Rprintf("iter = %d  f = %.8f  ng = %.8f\r", x.iterations, x.f, x.ng);
      R_FlushConsole();
      R_ProcessEvents();
    }
    if(x.ss == 0.0) break;

    x.dir = -x.rg;
    if(memory == 0) continue;

    arma::vec delta = x.parameters - old_parameters;
    arma::vec grad_delta = x.rg - old_rg;
    const double ys = arma::dot(grad_delta, delta);
    const double yy = arma::dot(grad_delta, grad_delta);
    if(!std::isfinite(ys) || !std::isfinite(yy) || yy <= 0.0 ||
       ys <= curvature_eps*arma::norm(delta, 2)*std::sqrt(yy)) continue;

    const double gamma = ys/yy;
    if(!std::isfinite(1.0/ys) || !std::isfinite(gamma) || gamma <= 0.0) continue;
    s[next] = std::move(delta);
    y[next] = std::move(grad_delta);
    rho[next] = 1.0/ys;
    next = (next + 1) % memory;
    count = std::min(count + 1, memory);

    // Two-loop recursion: newest to oldest, then oldest to newest.
    arma::vec z = x.rg;
    for(int j=1; j <= count; ++j) {
      const int i = (next + memory - j) % memory;
      alpha[i] = rho[i]*arma::dot(s[i], z);
      z -= alpha[i]*y[i];
    }
    z *= gamma;
    for(int j=count; j > 0; --j) {
      const int i = (next + memory - j) % memory;
      const double beta = rho[i]*arma::dot(y[i], z);
      z += s[i]*(alpha[i] - beta);
    }
    if(z.is_finite() && arma::dot(z, x.rg) > 0.0) x.dir = -z;
    else count = 0;

  }

  optim_result result = std::make_tuple(x.parameters,
                                        x.transparameters,
                                        x.f,
                                        x.iterations,
                                        x.convergence,
                                        x.ng, x.rg, x.g, x.dir);

  return result;

}
