/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 05/09/2026
 */

#include <memory>
#include "../component_cleanup.h"

#include "step_update/armijo.h"
#include "step_update/wolfe.h"
#include "step_update/tcg.h"

typedef std::tuple<arma::vec, arma::vec, double, int, bool, double, arma::mat,
                   arma::vec, arma::vec> optim_result;

#include "algorithm/gradient_descent.h"
#include "algorithm/lbfgs.h"
#include "algorithm/newtonTR.h"

class optim {

public:

  virtual ~optim() = default;

  virtual optim_result optimize(arguments_optim x,
                                std::vector<transformations*>& xtransforms,
                                std::vector<manifolds*>& xmanifolds,
                                std::vector<estimators*>& xestimators) = 0;

};

// Define EM after optim because EM uses optim in the M step:
#include "algorithm/EM.h"

// Riemannian Gradient Descent:

class RGD:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return gd(x, xtransforms, xmanifolds, xestimators);

  }

};

// L-BFGS:

class LBFGS:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return lbfgs(x, xtransforms, xmanifolds, xestimators);

  }

};

// Riemannian Newton Trust-Region:

class RNTR:public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return ntr(x, xtransforms, xmanifolds, xestimators);

  }

};

// Expectation-Maximization:

class EM: public optim {

public:

  optim_result optimize(arguments_optim x,
                        std::vector<transformations*>& xtransforms,
                        std::vector<manifolds*>& xmanifolds,
                        std::vector<estimators*>& xestimators) {

    return em(x, xtransforms, xmanifolds, xestimators);

  }

};

optim* choose_optim(arguments_optim& x, Rcpp::List control_optimizer,
                     const arma::vec& parameters,
                     const arma::vec& transparameters,
                     bool random_direction = true) {

  // Store the indices relating parameters and transformed parameters:
  arma::uvec transparam2param = control_optimizer["transparam2param"];
  x.transparam2param = transparam2param;

  // Store the parameters:
  x.nparam = parameters.n_elem;
  x.parameters = parameters;
  x.dir.set_size(x.nparam); x.dir.zeros();

  // Store the transformed parameters:
  x.ntransparam = transparameters.n_elem;
  x.transparameters = transparameters;

  if(x.transparam2param.n_elem != parameters.n_elem ||
     (x.transparam2param.n_elem > 0L &&
      x.transparam2param.max() >= transparameters.n_elem)) {
    Rcpp::stop("The parameter vectors do not match transparam2param.");
  }

  x.transparameters(x.transparam2param) = x.parameters;
  x.transparameters_init = x.transparameters;

  // Initial values for dparameters:
  if(control_optimizer.containsElementNamed("dparameters")) {
    arma::vec dparameters = control_optimizer["dparameters"];
    x.dparameters = dparameters;
  } else if(random_direction) {
    x.dparameters = arma::randu(x.nparam);
  } else {
    x.dparameters.zeros(x.nparam);
  }
  if(x.dparameters.n_elem != parameters.n_elem) {
    Rcpp::stop("dparameters must have one element per free parameter.");
  }

  x.dtransparameters.set_size(x.ntransparam);
  x.dtransparameters.zeros();
  x.dtransparameters(x.transparam2param) = x.dparameters;

  // Initialize objects:
  x.grad.set_size(x.ntransparam);
  x.dgrad.set_size(x.ntransparam);
  x.g.set_size(x.nparam);
  x.rg.set_size(x.nparam);
  x.dH.set_size(x.nparam);
  // x.hess.set_size(x.ntransparam, x.ntransparam);

  // Pass optimization parameters to the x structure:
  if(control_optimizer.containsElementNamed("opt")) {
    std::string opt = control_optimizer["opt"];
    x.optimizer = opt;
  }
  if(control_optimizer.containsElementNamed("maxit")) {
    int maxit = control_optimizer["maxit"];
    x.maxit = maxit;
  }
  if(control_optimizer.containsElementNamed("step_maxit")) {
    int step_maxit = control_optimizer["step_maxit"];
    x.step_maxit = step_maxit;
  }
  if(control_optimizer.containsElementNamed("step_eps")) {
    double step_eps = control_optimizer["step_eps"];
    x.step_eps = step_eps;
  }
  if(control_optimizer.containsElementNamed("df_eps")) {
    double df_eps = control_optimizer["df_eps"];
    x.df_eps = df_eps;
  }
  if(control_optimizer.containsElementNamed("M")) {
    int M = control_optimizer["M"];
    x.M = M;
  }
  if(control_optimizer.containsElementNamed("rstarts")) {
    int rstarts = control_optimizer["rstarts"];
    x.rstarts = rstarts;
  }
  if(control_optimizer.containsElementNamed("cores")) {
    int cores = control_optimizer["cores"];
    x.cores = cores;
  }
  if(control_optimizer.containsElementNamed("c2")) {
    double c2 = control_optimizer["c2"];
    x.c2 = c2;
  }
  if(control_optimizer.containsElementNamed("ss")) {
    double ss = control_optimizer["ss"];
    x.ss = ss;
  }
  if(control_optimizer.containsElementNamed("ss_fac")) {
    double ss_fac = control_optimizer["ss_fac"];
    x.ss_fac = ss_fac;
  }
  if(control_optimizer.containsElementNamed("ss_min")) {
    double ss_min = control_optimizer["ss_min"];
    x.ss_min = ss_min;
  }
  if(control_optimizer.containsElementNamed("eps")) {
    double eps = control_optimizer["eps"];
    x.eps = eps;
  }
  if(control_optimizer.containsElementNamed("tcg_maxit")) {
    double tcg_maxit = control_optimizer["tcg_maxit"];
    x.tcg_maxit = tcg_maxit;
  }
  if(control_optimizer.containsElementNamed("print")) {
    bool print = control_optimizer["print"];
    x.print = print;
  }
  if(control_optimizer.containsElementNamed("print_interval")) {
    int print_interval = control_optimizer["print_interval"];
    x.print_interval = print_interval;
  }
  if(control_optimizer.containsElementNamed("pick")) {
    int pick = control_optimizer["pick"];
    x.pick = pick; // Pick the "pick" number of rstarts with minimum objective
  }
  if(control_optimizer.containsElementNamed("idx_transforms")) {
    // Compute the jacobians and update the vcov of the parameters that are in
    // the control_transform structures indexed by idx_transforms:
    arma::uvec idx_transforms = control_optimizer["idx_transforms"];
    x.idx_transforms = idx_transforms;
  }
  if(control_optimizer.containsElementNamed("mstep_maxit")) {
    int mstep_maxit = control_optimizer["mstep_maxit"];
    x.mstep_maxit = mstep_maxit;
  }
  if(control_optimizer.containsElementNamed("mstep_eps")) {
    double mstep_eps = control_optimizer["mstep_eps"];
    x.mstep_eps = mstep_eps;
  }
  if(control_optimizer.containsElementNamed("mopt")) {
    std::string mopt = control_optimizer["mopt"];
    x.mopt = mopt;
  }

  // Select the optimization algorithm and set defaults:

  optim* algorithm;
  if(x.optimizer == "grad") {

    if(control_optimizer.containsElementNamed("c1")) {

      double c1 = control_optimizer["c1"];
      x.c1 = c1;

    } else {

      x.c1 = 0.5;

    }

    algorithm = new RGD();

  } else if(x.optimizer == "lbfgs") {

    if(control_optimizer.containsElementNamed("c1")) {

      double c1 = control_optimizer["c1"];
      x.c1 = c1;

    } else {

      // x.c1 = 10e-04;
      x.c1 = 0.5;

    }

    algorithm = new LBFGS();

  } else if(x.optimizer == "newton") {

    algorithm = new RNTR();

  // } else if(x.optimizer == "em") {
  //
  //   algorithm = new EM();
  //
  // } else if(x.optimizer == "em-lbfgs") {
  //
  //   algorithm = new EM();

  } else if(x.optimizer == "em") {

    algorithm = new EM();

  } else {

    Rf_error("Available optimization routines: \n grad, lbfgs, newton, em");

  }

  return algorithm;

}

// Starting values for model fitting:
optim* choose_optim(arguments_optim& x, Rcpp::List control_optimizer) {

  std::vector<arma::vec> parameters = control_optimizer["parameters"];
  std::vector<arma::vec> transparameters = control_optimizer["transparameters"];

  return choose_optim(x, control_optimizer, parameters[0], transparameters[0]);

}

// Fitted values for post-estimation computations:
optim* choose_optim(arguments_optim& x, Rcpp::S4 fit,
                     bool random_direction = true) {

  Rcpp::List modelInfo = fit.slot("modelInfo");
  Rcpp::List control_optimizer = modelInfo["control_optimizer"];
  Rcpp::List control_transform = modelInfo["control_transform"];
  Rcpp::List Optim = fit.slot("Optim");

  arma::vec parameters = Optim["parameters"];
  arma::vec transparameters = Optim["transparameters"];
  Rcpp::List initial = control_optimizer["transparameters"];
  arma::vec transparameters_init = initial[0];

  if(transparameters.n_elem != transparameters_init.n_elem) {
    Rcpp::stop("Optim$transparameters does not match the model coordinates.");
  }

  // Recompute outputs from their original seeds, not their fitted values.
  // This is essential for transformations that accumulate log-likelihoods.
  // All other values, including fixed inputs, are taken from Optim.
  for(int i=0; i < control_transform.size(); ++i) {
    Rcpp::List setup = control_transform[i];
    std::vector<arma::uvec> indices_out = setup["indices_out"];
    for(const arma::uvec& indices : indices_out) {
      transparameters.elem(indices) = transparameters_init.elem(indices);
    }
  }

  if(!control_optimizer.containsElementNamed("idx_transforms")) {
    x.idx_transforms.set_size(x.ntransforms);
    for(int i=0; i < x.ntransforms; ++i) x.idx_transforms[i] = i;
  }

  return choose_optim(x, control_optimizer, parameters, transparameters,
                       random_direction);

}
