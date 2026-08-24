/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

// L-BFGS algorithm:

optim_result lbfgs(arguments_optim x,
                   std::vector<transformations*>& xtransforms,
                   std::vector<manifolds*>& xmanifolds,
                   std::vector<estimators*>& xestimators) {

  product_manifold* final_manifold;
  product_transform* final_transform;
  product_estimator* final_estimator;

  // Ensure initial parameters are ok:
  final_manifold->param(x, xmanifolds);
  final_manifold->retr(x, xmanifolds);
  final_manifold->param(x, xmanifolds);
  final_transform->transform(x, xtransforms);
  final_estimator->param(x, xestimators);
  x.iterations = 0L;

  // Parameterization
  // final_estimator->param(x, xestimators);
  final_estimator->F(x, xestimators);
  // update the gradient
  final_estimator->G(x, xestimators);
  final_transform->update_grad(x, xtransforms);
  // Riemannian gradient
  // final_manifold->param(x, xmanifolds);
  final_manifold->proj(x, xmanifolds);
  x.dir = -x.rg;
  x.inprod = arma::dot(-x.dir, x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p1 = x.parameters.size();
  arma::mat B(p1, p1, arma::fill::eye);

  std::vector<arma::vec> s(x.maxit), y(x.maxit);
  std::vector<double> p(x.maxit), alpha(x.maxit), beta(x.maxit);

  x.convergence = false;

  do {

    // x.ss *= 2;

    int k = x.iterations;
    ++x.iterations;

    arma::uvec seq(2);
    seq[0] = x.M; seq[1] = k;
    int min = seq.min();
    arma::vec max(2);
    max[0] = min; max[1] = 0;
    int m = max.max();

    arma::vec old_parameters = x.parameters;
    arma::vec old_rg = x.rg;

    if(!x.dir.is_finite() || !x.rg.is_finite()) {
      x.convergence = false;
      break;
    }

    double dir_dot_grad = arma::dot(x.dir, x.rg);
    if(!std::isfinite(dir_dot_grad) || dir_dot_grad >= 0.0) {
      x.dir = -x.rg;
    }

    // Start every L-BFGS line search from a unit step:
    x.ss = 1;
    // armijo(x, final_manifold, final_estimator, xmanifolds, xestimators);
    wolfe(x, xtransforms, xmanifolds, xestimators);

    final_estimator->param(x, xestimators); // Necessary¿?
    final_estimator->G(x, xestimators);
    final_transform->update_grad(x, xtransforms);
    // After the retraction in armijo you need to param the manifolds:
    final_manifold->param(x, xmanifolds);
    final_manifold->proj(x, xmanifolds);

    // Find a new direction:
    arma::vec q = arma::vectorise(x.rg);
    s[k] = arma::vectorise(x.parameters - old_parameters);
    y[k] = arma::vectorise(x.rg - old_rg);
    double dot_ys = arma::dot(y[k], s[k]);
    double norm_s = arma::norm(s[k], 2);
    double norm_y = arma::norm(y[k], 2);
    double curvature_eps = std::sqrt(std::numeric_limits<double>::epsilon());

    if(!std::isfinite(dot_ys) ||
       !std::isfinite(norm_s) ||
       !std::isfinite(norm_y) ||
       dot_ys <= curvature_eps*norm_s*norm_y) {

      s[k].zeros(p1);
      y[k].zeros(p1);
      p[k] = 0.0;
      x.dir = -x.rg;
      continue;

    }

    p[k] = 1/dot_ys;

    if(!std::isfinite(p[k])) {

      s[k].zeros(p1);
      y[k].zeros(p1);
      p[k] = 0.0;
      x.dir = -x.rg;
      continue;

    }

    for(int i=k; i > (k-m-1); --i) {

      alpha[i] = p[i]*arma::dot(s[i], q);
      q -= alpha[i] * y[i];

    }

    if(!q.is_finite()) {

      for(int i=0; i <= k; ++i) {
        s[i].zeros(p1);
        y[i].zeros(p1);
        p[i] = 0.0;
      }

      x.dir = -x.rg;
      continue;

    }

    double dot_yy = arma::dot(y[k], y[k]);
    double gamma = dot_ys / dot_yy;

    if(!std::isfinite(gamma) || gamma <= arma::datum::eps) {

      s[k].zeros(p1);
      y[k].zeros(p1);
      p[k] = 0.0;
      x.dir = -x.rg;
      continue;

    }

    arma::mat H0 = gamma*B;
    arma::vec z = H0 * q;

    for(int i=(k-m); i < (k+1); ++i) {

      beta[i] = p[i]*arma::dot(y[i], z);
      z += s[i] * (alpha[i] - beta[i]);

    }

    if(!z.is_finite()) {

      for(int i=0; i <= k; ++i) {
        s[i].zeros(p1);
        y[i].zeros(p1);
        p[i] = 0.0;
      }

      x.dir = -x.rg;

    } else {

      x.dir = -z;

    }

    dir_dot_grad = arma::dot(x.dir, x.rg);

    if(!x.dir.is_finite() ||
       !std::isfinite(dir_dot_grad) ||
       dir_dot_grad >= 0.0) {

      for(int i=0; i <= k; ++i) {
        s[i].zeros(p1);
        y[i].zeros(p1);
        p[i] = 0.0;
      }

      x.dir = -x.rg;

    }

    // Check convergence:
    // x.inprod = arma::accu(x.dir % x.dir);
    // x.inprod = arma::accu(x.rg % x.rg);
    x.inprod = arma::dot(-x.dir, x.rg);
    x.ng = std::sqrt(x.inprod);
    x.max_rg = arma::abs(x.rg).max();

    // Print info:
    if((x.print & x.rstarts) == 1L) {

      if (x.iterations % x.print_interval == 0) {
        Rprintf("iter = %d  f = %.8f  ng = %.8f\r",
                x.iterations, x.f, x.ng);
        R_FlushConsole();
        R_ProcessEvents();
      }

    }

    if(!std::isfinite(x.ng)) {
      x.convergence = false;
      break;
    }
    // if ((x.ng/std::sqrt(x.rg.n_elem)) < x.eps && x.max_rg < x.eps) {
    if ((x.ng/std::sqrt(x.rg.n_elem)) < x.eps) {
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
