/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "structures.h"
// #include "auxiliary_manifolds.h"
// #include "manifold.h"
// #include "auxiliary_criteria.h"
// #include "criteria.h"

// Line-search algorithm satisfying the armijo condition:

void armijo(arguments_efa& x, efa_manifold *manifold, efa_criterion *criterion,
            double ss_fac, double ss_min) {

  x.ss = std::max(ss_min, x.ss * ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  int iteration = 0;
  arma::vec parameters = x.parameters;
  x.inprod = arma::accu(x.dir % x.rg);

  do{

    ++iteration;
    x.parameters = parameters + x.ss*x.dir;
    // Projection onto the manifold
    manifold->retr(x); // update x.T
    // Parameterization
    manifold->param(x); // update x.L, x.Phi and x.Inv_T
    criterion->F(x);
    double df = x.f - f0;
    if (df < x.c1 * x.ss * x.inprod || // armijo condition
        x.ss < 1e-09) break;
    x.ss *= x.rho;

  } while (iteration <= x.armijo_maxit);

}

void strong_wolfe(arguments_efa& x, efa_manifold *manifold, efa_criterion *criterion) {

  x.inprod = arma::accu(x.dir % x.rg);
  double y_prev, alpha_prev, alpha_lo, alpha_hi;
  arguments_efa x2 = x;
  arguments_efa xlo = x;
  double beta = 10e-04, sigma = 0.1;
  int inner_iter = 0;

  // bracket phase
  do{

    ++inner_iter;

    x2.parameters = x.parameters + x.ss*x.dir;
    // Parameterization
    manifold->param(x);
    criterion->F(x2);
    if((x2.f > (x.f + beta*x.ss*x.inprod)) | (inner_iter > 1 & (x2.f >= y_prev))) {
      alpha_lo = alpha_prev;
      alpha_hi = x.ss;
      break;
    }

    // update gradient
    criterion->G(x2);
    manifold->grad(x2);
    // criterion->G(x2);
    manifold->proj(x2);
    x2.inprod = arma::accu(x.dir % x2.rg);

    if(abs(x2.inprod) <= (-sigma * x.inprod)) {
      return;
    } else if(x2.inprod >= 0) {
      alpha_lo = x.ss;
      alpha_hi = alpha_prev;
      break;
    }

    y_prev = x2.f;
    alpha_prev = x.ss;
    x.ss *= 2;

  } while (true);

  // zoom phase
  for(int i=0; i < 5; ++i) {

    xlo.parameters = x.parameters + alpha_lo*x.dir;
    criterion->F(xlo);
    x.ss = 0.5*(alpha_lo + alpha_hi);
    x2.parameters = x.parameters + x.ss*x.dir;
    // Parameterization
    manifold->param(x); // update x.L, x.Phi and x.Inv_T
    criterion->F(x);

    if((x2.f > (x.f + beta * x.ss * x.inprod)) | (x2.f >= xlo.f)) {
      alpha_hi = x.ss;
    } else {
      // update gradient
      criterion->G(x);
      manifold->grad(x);
      manifold->proj(x2);
      x2.inprod = arma::accu(x.dir % x2.rg);
      if(abs(x2.inprod) <= (-sigma*x.inprod)) {
        return;
      } else if((x2.inprod*(alpha_hi - alpha_lo)) >= 0) {
        alpha_hi = alpha_lo;
      }
      alpha_lo = x.ss;
    }
  }

  return;

}

// Conjugate-gradient method to solve the Riemannian Newton equation:

void tcg(arguments_efa x, efa_manifold *manifold, efa_criterion *criterion,
         arma::vec& dir, bool& att_bnd, double ng, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  dir.zeros();
  arma::vec dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dparameters = -x.rg; // Initial search direction
  arma::vec r = x.dparameters; // Initial residual
  double rr = ng * ng;
  double tol = ng * std::min(pow(ng, c[0]), c[1]);

  int iter = 0;

  do{

    // Differential of L and P
    manifold->dparam(x);

    // Differential of the gradient of L and P
    criterion->dG(x);

    // Differential of g
    manifold->dgrad(x);

    // Riemannian hessian
    manifold->hess(x);

    dHd = arma::accu(x.dparameters % x.dH);

    if(dHd <= 0) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters), 2 * arma::accu(dir % x.dparameters),
                      arma::accu(dir % dir) - rad * rad); // Solve equation 39
      dir = dir + tau * x.dparameters;
      att_bnd = true;

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = dir;
    dir = dir + alpha * x.dparameters; // update proposal

    if (sqrt(arma::accu(dir % dir)) >= rad) {

      tau = root_quad(arma::accu(x.dparameters % x.dparameters), 2 * arma::accu(dir0 % x.dparameters),
                      arma::accu(dir0 % dir0) - rad * rad); // Solve equation 39
      dir = dir0 + tau * x.dparameters;
      att_bnd = true;

      break;

    }

    r = r - alpha * x.dH; // update gradient
    rr = arma::accu(r % r);

    if (sqrt(rr) < tol) {

      att_bnd = false;
      break;

    }

    beta = rr / rr0;
    x.dparameters = r + beta * x.dparameters;
    iter = iter + 1;

  } while (iter < 5);

}

// Newton Trust-region algorithm:

efa_NTR ntr(arguments_efa x, efa_manifold *manifold, efa_criterion *criterion) {

  /*
   * Riemannian trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T

  // Objective
  criterion->F(x); // update x.f, x.L2, x.IgCL2N, x.term, x.f1 and x.f2

  // Rcpp::Rcout << "x.f = " << x.f << std::endl;

  // Gradient wrt L
  criterion->G(x); // update x.gL, x.gP, x.f1, x.f2 and x.LoL2

  // Rcpp::Rcout << "x.gL = " << x.gL << std::endl;

  // Gradient wtr T
  manifold->grad(x); // update x.g

  // Riemannian gradient
  manifold->proj(x); // update x.rg and x.A

  // Differential of the gradient of L and P
  // criterion->dG(x); // update dgL and dgP

  // Rcpp::Rcout << "x.dgL = " << x.dgL << std::endl;

  double ng = sqrt(arma::accu(x.rg % x.rg));

  double max_rad = 10;

  arma::vec fac_rad(2);
  fac_rad[0] = 0.25;
  fac_rad[1] = 2;

  arma::vec crit_goa(3);
  crit_goa[0] = 0.2;
  crit_goa[1] = 0.25;
  crit_goa[2] = 0.75;

  arma::vec c(2);
  c[0] = 1;
  c[1] = 0.01;

  double rad = 1;
  bool att_bnd = false;

  x.iteration = 0;
  double goa, preddiff;

  arguments_efa new_x;
  arma::vec dir(x.p);

  do {

    if (ng < x.eps) break;

    ++x.iteration;

    // subsolver
    tcg(x, manifold, criterion, dir, att_bnd, ng, c, rad);
    x.dparameters = dir;
    new_x = x;
    new_x.parameters += dir;

    // Projection onto the manifold
    manifold->retr(new_x); // update x.T

    // Differential of L and P
    manifold->dparam(x); // update x.dL, x.dP and Inv_T_dt

    // Differential of the gradient of L and P
    criterion->dG(x); // update dgL and dgP

    // Differential of g
    manifold->dgrad(x); // update dg

    // Riemannian hessian
    manifold->hess(x); // update dH

    preddiff = - arma::accu(x.dparameters % ( x.rg + 0.5 * x.dH) );

    // Parameterization
    manifold->param(new_x);

    // objective
    criterion->F(new_x);

    if ( std::abs(preddiff) <= arma::datum::eps ) {

      goa = 1;

    } else {

      goa = (x.f - new_x.f) / preddiff;

    }
    if (goa < crit_goa[1]) {

      rad = fac_rad[0] * rad;

    } else if (goa > crit_goa[2] && att_bnd) {

      rad = std::min(fac_rad[1] * rad, max_rad);

    }

    // accepted iteration
    if (goa > crit_goa[0]) {

      x = new_x;

      // update gradient
      criterion->G(x);
      manifold->grad(x);

      // Riemannian gradient
      manifold->proj(x);

      ng = sqrt(arma::accu(x.rg % x.rg));

    }

  } while (x.iteration <= x.maxit);


  x.convergence = true;
  if(x.iteration > x.maxit) {

    x.convergence = false;

  }

  criterion->outcomes(x);

  NTR result = std::make_tuple(x.lambda, x.uniquenesses, x.Rhat,
                               x.f, x.iteration, x.convergence);

  return result;

}

// Gradient descent algorithm:

efa_NTR gd(arguments_efa x, efa_manifold *manifold, efa_criterion *criterion) {

  double ss_fac = 2, ss_min = 0.1;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);
  // update gradient
  criterion->G(x);
  manifold->grad(x);
  // criterion->G(x);
  // Riemannian gradient
  manifold->proj(x);
  // x.ss = 1;

  do{

    // x.ss *= 2;
    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    if (x.ng < x.eps) break;

    ++x.iterations;

    if(x.search == "back") {
      armijo(x, manifold, criterion, ss_fac, ss_min);
    } else if(x.search == "wolfe") {
      strong_wolfe(x, manifold, criterion);
      x.parameters += x.ss * x.dir; // For strong_wolve
      manifold->param(x); // update x.L, x.Phi and x.Inv_T
      criterion->F(x); // For strong_wolve
    }

    // update gradient
    criterion->G(x);
    manifold->grad(x);
    // criterion->G(x);
    // Riemannian gradient
    manifold->proj(x);

  } while (x.iterations < x.maxit);

  x.convergence = true;
  if(x.iterations > x.maxit) {

    x.convergence = false;

  }

  criterion->outcomes(x);

  efa_NTR result = std::make_tuple(x.lambda, x.uniquenesses, x.Rhat,
                                   x.f, x.iterations, x.convergence);

  return result;

}

// L-BFGS algorithm:

efa_NTR lbfgs(arguments_efa x, efa_manifold *manifold, efa_criterion *criterion) {

  double ss_fac = 2, ss_min = 0.1;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);
  // update the gradient
  criterion->G(x);  // Update the gradient wrt x.L and x.Phi
  manifold->grad(x);  // Update the gradient wrt x.T
  // criterion->G(x);
  // Riemannian gradient
  manifold->proj(x);  // Update the Riemannian gradient x.rg
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;
  int p1 = x.parameters.size();
  arma::mat B(p1, p1, arma::fill::eye);

  std::vector<arma::vec> s(x.maxit), y(x.maxit);
  std::vector<double> p(x.maxit), alpha(x.maxit), beta(x.maxit);

  x.convergence = false;

  do{

    // x.ss *= 2;

    int k = x.iterations;
    arma::uvec seq(2);
    seq[0] = x.M; seq[1] = k;
    int min = seq.min();
    arma::vec max(2);
    max[0] = min; max[1] = 0;
    int m = max.max();

    arma::vec old_parameters = x.parameters;
    arma::vec old_rg = x.rg;

    armijo(x, manifold, criterion, ss_fac, ss_min);

    // update gradient
    criterion->G(x);
    manifold->grad(x);
    // criterion->G(x);
    // Riemannian gradient
    manifold->proj(x);

    arma::vec q = arma::vectorise(x.rg);
    s[k] = arma::vectorise(x.parameters - old_parameters);
    y[k] = arma::vectorise(x.rg - old_rg);
    p[k] = 1/arma::accu(y[k] % s[k]);

    for(int i=k; i > (k-m-1); --i) {

      alpha[i] = p[i]*arma::accu(s[i] % q);
      q -= alpha[i] * y[i];

    }

    double gamma = arma::accu(s[k] % y[k]) / arma::accu(y[k] % y[k]);
    arma::mat H0 = gamma*B;
    arma::vec z = H0 * q;

    for(int i=(k-m); i < (k+1); ++i) {

      beta[i] = p[i]*arma::accu(y[i] % z);
      z += s[i] * (alpha[i] - beta[i]);

    }

    // z.reshape(p1, p2);
    z = arma::vectorise(z);

    x.dir = -z;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iterations;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iterations < x.maxit);

  criterion->outcomes(x);

  efa_NTR result = std::make_tuple(x.lambda, x.uniquenesses, x.Rhat,
                                   x.f, x.iterations, x.convergence);

  return result;

}

