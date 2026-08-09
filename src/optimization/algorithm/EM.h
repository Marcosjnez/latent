/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 09/08/2026
 */

// Expectation-Maximization:

optim_result moptim(arguments_optim x,
                    std::vector<transformations*>& xtransforms,
                    std::vector<manifolds*>& xmanifolds,
                    std::vector<estimators*>& xestimators) {

  // Optimization algorithm within the M step

  if(x.mopt == "grad") {

    return gd(x, xtransforms, xmanifolds, xestimators);

  } else if(x.mopt == "lbfgs") {

    return lbfgs(x, xtransforms, xmanifolds, xestimators);

  } else if(x.mopt == "newton") {

    return ntr(x, xtransforms, xmanifolds, xestimators);

  } else {

    Rf_error(
      "Available optimization routines for the M step: "
      "grad, lbfgs, newton"
    );

  }

}

optim_result em(arguments_optim x,
                std::vector<transformations*>& xtransforms,
                std::vector<manifolds*>& xmanifolds,
                std::vector<estimators*>& xestimators) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  // Ensure the starting parameters satisfy the manifold constraints:
  final_manifold.param(x, xmanifolds);
  final_manifold.retr(x, xmanifolds);
  final_manifold.param(x, xmanifolds);

  // Evaluate the initial observed-data objective:
  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.observed_F(x, xestimators);

  double old_f = x.f;
  double diff = arma::datum::inf;

  x.iterations = 0L;
  x.convergence = false;

  do {

    ++x.iterations;

    // ------------------------------------------------------------
    // E step
    // ------------------------------------------------------------

    // param() above calculated the current posterior probabilities.
    // Freeze them throughout the following M step:
    final_estimator.E(x, xestimators);

    // ------------------------------------------------------------
    // M step
    // ------------------------------------------------------------

    arguments_optim x_mstep = x;

    // Use separate convergence controls for the numerical M step:
    x_mstep.maxit = x.mstep_maxit;
    x_mstep.eps = x.mstep_eps;

    x_mstep.iterations = 0L;
    x_mstep.convergence = false;
    // Avoid printing the inner L-BFGS iterations:
    x_mstep.print = false;

    // lcaEM::F() and lcaEM::G() contain the Q-function and its
    // gradient, so the ordinary L-BFGS optimizer can be reused
    // without modification:
    // optim_result mstep = lbfgs(x_mstep, xtransforms, xmanifolds, xestimators);
    optim_result mstep = moptim(x_mstep, xtransforms, xmanifolds, xestimators);

    // Retain the M-step parameter estimates:
    x.parameters = std::get<0>(mstep);
    x.transparameters = std::get<1>(mstep);
    x.dir = std::get<8>(mstep);

    // ------------------------------------------------------------
    // Observed-data likelihood
    // ------------------------------------------------------------

    // Recompute all transformed parameters with the updated estimates:
    final_transform.transform(x, xtransforms);

    // Recompute observed likelihood and posterior probabilities:
    final_estimator.param(x, xestimators);
    final_estimator.observed_F(x, xestimators);

    double new_f = x.f;

    // Convergence of the outer EM algorithm:
    diff =
      std::abs(new_f - old_f) /
        (1.0 + std::abs(old_f));

    x.ng = diff;

    // Print information:
    if(x.print) {

      if(x.iterations % x.print_interval == 0L) {

        Rprintf("EM iter = %d  f = %.8f  diff = %.8f\r",
                x.iterations, x.f, diff);

        R_FlushConsole();
        R_ProcessEvents();

      }

    }

    if(!std::isfinite(new_f)) {

      x.convergence = false;
      break;

    }

    if(diff < x.eps) {

      x.convergence = true;
      break;

    }

    old_f = new_f;

  } while(x.iterations < x.maxit);

  // ------------------------------------------------------------
  // Final observed-data objective and gradient
  // ------------------------------------------------------------

  final_transform.transform(x, xtransforms);

  final_estimator.param(x, xestimators);
  final_estimator.observed_F(x, xestimators);
  final_estimator.observed_G(x, xestimators);

  // Propagate the observed-data gradient through all transformations:
  final_transform.update_grad(x, xtransforms);

  // Project onto the product manifold:
  final_manifold.param(x, xmanifolds);
  final_manifold.proj(x, xmanifolds);

  optim_result result =
    std::make_tuple(x.parameters,
                    x.transparameters,
                    x.f,
                    x.iterations,
                    x.convergence,
                    x.ng,
                    x.rg,
                    x.g,
                    x.dir);

  return result;

}
