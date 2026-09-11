/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

#ifndef LATENT_NEWTON_H
#define LATENT_NEWTON_H

// Riemannian Newton-CG with trust-region or line-search step selection.
optim_result newton(arguments_optim x,
                     std::vector<transformations*>& xtransforms,
                     std::vector<manifolds*>& xmanifolds,
                     std::vector<estimators*>& xestimators) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;
  std::unique_ptr<stepsize> step(choose_stepsize(x.step, true));

  final_manifold.param(x, xmanifolds);
  final_manifold.retr(x, xmanifolds);
  final_manifold.param(x, xmanifolds);
  final_transform.transform(x, xtransforms);
  final_estimator.param(x, xestimators);
  final_estimator.F(x, xestimators);
  final_estimator.G(x, xestimators);
  final_transform.update_grad(x, xtransforms);
  final_manifold.proj(x, xmanifolds);

  x.inprod = arma::accu(x.rg % x.rg);
  x.ng = std::sqrt(x.inprod);
  x.iterations = 0L;
  x.convergence = false;

  while(x.iterations < x.maxit) {

    if(!std::isfinite(x.f) || !x.rg.is_finite() || !std::isfinite(x.ng)) break;
    if(x.rg.n_elem == 0L || x.ng/std::sqrt(x.rg.n_elem) < x.eps) {
      x.convergence = true;
      break;
    }

    if(!step->uses_trust_region()) {

      tcg(x, xtransforms, xmanifolds, xestimators);
      double slope = arma::dot(x.dir, x.rg);
      if(!x.dir.is_finite() || !std::isfinite(slope) || slope >= 0.0) {
        x.dir = -x.rg;
      }

    }

    step->update(x, xtransforms, xmanifolds, xestimators);

    // The selected point may be the old point after a rejected TRUST proposal.
    final_estimator.G(x, xestimators);
    final_transform.update_grad(x, xtransforms);
    final_manifold.param(x, xmanifolds);
    final_manifold.proj(x, xmanifolds);
    x.inprod = arma::accu(x.rg % x.rg);
    x.ng = std::sqrt(x.inprod);
    x.max_rg = arma::abs(x.rg).max();
    ++x.iterations;

    if(x.print && x.rstarts == 1L && x.print_interval > 0L &&
       x.iterations % x.print_interval == 0L) {
      Rprintf("iter = %d  f = %.8f  ng = %.8f\r", x.iterations, x.f, x.ng);
      R_FlushConsole();
      R_ProcessEvents();
    }

    // Retain Newton's gradient-based stopping rule, also on the last iteration.
    if(std::isfinite(x.f) && std::isfinite(x.ng) &&
       x.ng/std::sqrt(x.rg.n_elem) < x.eps) {
      x.convergence = true;
      break;
    }

  }

  return std::make_tuple(x.parameters, x.transparameters, x.f, x.iterations,
                         x.convergence, x.ng, x.rg, x.g, x.dir);

}

#endif
