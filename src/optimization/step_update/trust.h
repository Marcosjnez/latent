/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

#ifndef LATENT_TRUST_STEP_H
#define LATENT_TRUST_STEP_H

#include "../direction/tcg.h"

// One trust-region proposal. The TRUST strategy retains rad between calls.
void trust(arguments_optim& x,
           std::vector<transformations*>& xtransforms,
           std::vector<manifolds*>& xmanifolds,
           std::vector<estimators*>& xestimators,
           double& rad, bool& att_bnd) {

  product_manifold final_manifold;
  product_transform final_transform;
  product_estimator final_estimator;

  const double max_rad = 10.0;
  const arma::vec fac_rad = {0.25, 2.0};
  const arma::vec crit_goa = {0.20, 0.25, 0.75};
  const arma::vec c = {1.0, 0.01};

  tcg(x, xtransforms, xmanifolds, xestimators, att_bnd, c, rad);

  if(!x.dir.is_finite()) {
    rad = fac_rad[0]*rad;
    x.dir.zeros();
    x.dparameters = x.dir;
    x.df = 0.0;
    return;
  }

  final_transform.dtransform(x, xtransforms);
  final_estimator.dG(x, xestimators);
  final_transform.update_dgrad(x, xtransforms);
  final_manifold.param(x, xmanifolds);
  final_manifold.proj(x, xmanifolds);
  final_manifold.hess(x, xmanifolds);

  double preddiff = -arma::accu(x.dparameters % (x.rg + 0.5*x.dH));

  if(!std::isfinite(preddiff)) {
    rad = fac_rad[0]*rad;
    x.df = 0.0;
    return;
  }

  arguments_optim new_x = x;
  new_x.parameters += x.dir;
  final_manifold.param(new_x, xmanifolds);
  final_manifold.retr(new_x, xmanifolds);
  final_manifold.param(new_x, xmanifolds);
  final_transform.transform(new_x, xtransforms);
  final_estimator.param(new_x, xestimators);
  final_estimator.F(new_x, xestimators);

  x.df = x.f - new_x.f;
  new_x.df = x.df;
  double goa = (std::abs(preddiff) <= arma::datum::eps) ?
    1.0 : x.df/preddiff;

  // Invalid proposals must not be accepted through the tiny-prediction case.
  if(!new_x.parameters.is_finite() || !new_x.transparameters.is_finite() ||
     !std::isfinite(new_x.f) || !std::isfinite(goa)) {
    goa = -arma::datum::inf;
  }

  if(goa < crit_goa[1]) {
    rad = fac_rad[0]*rad;
  } else if(goa > crit_goa[2] && att_bnd) {
    rad = std::min(fac_rad[1]*rad, max_rad);
  }

  if(goa > crit_goa[0]) {

    x = new_x;

  } else {

    // Trial evaluation changed the shared caches, even though x was retained.
    final_manifold.param(x, xmanifolds);
    final_transform.transform(x, xtransforms);
    final_estimator.param(x, xestimators);
    final_estimator.F(x, xestimators);
    x.df = 0.0;

  }

}

#endif
