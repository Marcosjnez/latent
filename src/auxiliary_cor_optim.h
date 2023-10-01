// Conjugate-gradient method to solve the Riemannian Newton equation in the Trust-Region algorithm:
// Line-search satisfying the armijo condition:

void cor_armijo(arguments_cor& x, cor_manifold *manifold,
            cor_criterion *criterion,
            double ss_fac, double ss_min, double max_iter,
            double c1, double c2, double eps) {

  x.ss = std::max(ss_min, x.ss * ss_fac);
  // x.ss = x.ss*2;
  double f0 = x.f;
  int iteration = 0;
  arma::mat X = x.T;
  x.inprod = arma::accu(x.dir % x.rg);

  do{

    ++iteration;
    x.T = X + x.ss*x.dir;
    // Projection onto the manifold
    manifold->retr(x); // update x.T
    // Parameterization
    manifold->param(x); // update x.L, x.Phi and x.Inv_T
    criterion->F(x);
    double df = x.f - f0;
    if (df < c1 * x.ss * x.inprod || x.ss < 1e-09) // armijo condition
      break;
    x.ss *= c2;

  } while (iteration <= max_iter);

  bool convergence = true;
  if(iteration > max_iter) {

    convergence = false;

  }

}

void cor_tcg(arguments_cor x, cor_manifold *manifold, cor_criterion *criterion,
         arma::mat& dir, bool& att_bnd, double ng, arma::vec c, double rad) {

  /*
   * Truncated conjugate gradient sub-solver for the trust-region sub-problem
   * From Liu (Algorithm 4; 2020)
   */

  dir.zeros();
  arma::mat dir0;

  double alpha, rr0, tau, beta, dHd;
  x.dT = -x.rg; // Initial search direction
  arma::mat r = x.dT; // Initial residual
  double rr = ng * ng;
  double tol = ng * std::min(pow(ng, c[0]), c[1]);

  int iter = 0;

  do{

    // Differential of L and P
    manifold->dcor(x);

    // Differential of the gradient of L and P
    criterion->dgcor(x);

    // Differential of g
    manifold->dgrad(x);

    // Riemannian hessian
    manifold->hess(x);

    dHd = arma::accu(x.dT % x.dH);

    if(dHd <= 0) {

      tau = root_quad(arma::accu(x.dT % x.dT), 2 * arma::accu(dir % x.dT),
                      arma::accu(dir % dir) - rad * rad); // Solve equation 39
      dir = dir + tau * x.dT;
      att_bnd = true;

      break;

    }

    rr0 = rr;
    alpha = rr0 / dHd;
    dir0 = dir;
    dir = dir + alpha * x.dT; // update proposal

    if (sqrt(arma::accu(dir % dir)) >= rad) {

      tau = root_quad(arma::accu(x.dT % x.dT), 2 * arma::accu(dir0 % x.dT),
                      arma::accu(dir0 % dir0) - rad * rad); // Solve equation 39
      dir = dir0 + tau * x.dT;
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
    x.dT = r + beta * x.dT;
    iter = iter + 1;

  } while (iter < 5);

}

// Newton Trust-region algorithm:
typedef std::tuple<arma::mat, arma::mat, double, int, bool> cor_NTR;
cor_NTR cor_ntr(arguments_cor x, cor_manifold *manifold, cor_criterion *criterion) {

  /*
   * Riemannian trust-region algorithm
   * From Liu (Algorithm 2; 2020)
   */

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T

  // Objective
  criterion->F(x); // update x.f, x.L2, x.IgCL2N, x.term, x.f1 and x.f2

  // Rcpp::Rcout << "x.f = " << x.f << std::endl;

  // Gradient wrt L and P
  criterion->gcor(x); // update x.gL, x.gP, x.f1, x.f2 and x.LoL2

  // Rcpp::Rcout << "x.gL = " << x.gL << std::endl;

  // Gradient wtr T
  manifold->grad(x); // update x.g

  // Riemannian gradient
  manifold->proj(x); // update x.rg and x.A

  // Differential of the gradient of L and P
  // criterion->dgLP(x); // update dgL and dgP

  // Rcpp::Rcout << "x.dgL = " << x.dgL << std::endl;

  x.ng = sqrt(arma::accu(x.rg % x.rg));

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

  arguments_cor new_x;
  arma::mat dir(x.q, x.q);

  x.convergence = false;

  do{

    // subsolver
    cor_tcg(x, manifold, criterion, dir, att_bnd, x.ng, c, rad);
    x.dT = dir;
    new_x = x;
    new_x.T += dir;

    // Projection onto the manifold
    manifold->retr(new_x); // update x.T

    // Differential of L and P
    manifold->dcor(x); // update x.dL, x.dP and Inv_T_dt

    // Differential of the gradient of L and P
    criterion->dgcor(x); // update dgL and dgP

    // Differential of g
    manifold->dgrad(x); // update dg

    // Riemannian hessian
    manifold->hess(x); // update dH

    preddiff = - arma::accu(x.dT % ( x.rg + 0.5 * x.dH) );

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

      // update the gradient
      criterion->gcor(x);
      manifold->grad(x);

      // Riemannian gradient
      manifold->proj(x);

      x.ng = sqrt(arma::accu(x.rg % x.rg));

    }

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);

  cor_NTR result = std::make_tuple(x.correlation, x.T, x.f, x.iteration, x.convergence);

  return result;

}

// Gradient descent algorithm:
cor_NTR cor_gd(arguments_cor x, cor_manifold *manifold, cor_criterion *criterion) {

  x.iteration = 0;
  double ss_fac = 2, ss_min = 0.1, c1 = 0.5, c2 = 0.5;

  // Parameterization
  manifold->param(x); // update x.L, x.Phi and x.Inv_T
  criterion->F(x);
  // update gradient
  criterion->gcor(x);
  manifold->grad(x);
  // Riemannian gradient
  manifold->proj(x);
  x.dir = -x.rg;
  x.inprod = arma::accu(-x.dir % x.rg);
  x.ng = sqrt(x.inprod);
  // x.ss = 1;

  x.convergence = false;

  do{

    // x.ss *= 2;

    cor_armijo(x, manifold, criterion, ss_fac, ss_min,
               10, c1, c2, x.eps);

    // update gradient
    criterion->gcor(x);
    manifold->grad(x);
    // Riemannian gradient
    manifold->proj(x);
    x.dir = -x.rg;
    x.inprod = arma::accu(-x.dir % x.rg);
    x.ng = sqrt(x.inprod);

    ++x.iteration;
    if (x.ng < x.eps) {
      x.convergence = true;
      break;
    }

  } while (x.iteration < x.maxit);

  cor_NTR result = std::make_tuple(x.correlation, x.T, x.f, x.iteration, x.convergence);

  return result;

}
