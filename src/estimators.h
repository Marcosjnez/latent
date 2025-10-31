/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 30/10/2025
 */

class estimators {

public:

  double f, loglik;
  arma::vec parameters, transparameters, dparameters;
  arma::mat g, dg;
  arma::mat grad, dgrad;
  arma::mat hess;
  arma::mat dparam_dS, modhessian;
  arma::mat posterior, latentloglik;
  arma::vec uniquenesses;
  // arma::mat lambda, phi, psi, model, residuals;
  int nhessian, nS, p, q;
  arma::vec weights, latentpars, loglatentpars; // P(X = c) // classes_hat
  std::vector<arma::mat> conditionals; // conditionals_hat

  std::vector<arma::uvec> indices;
  arma::mat freqs;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void param(arguments_optim& x) = 0;

  virtual void F(arguments_optim& x) = 0;

  virtual void G(arguments_optim& x) = 0;

  virtual void dG(arguments_optim& x) = 0;

  virtual void E(arguments_optim& x) = 0;

  virtual void M(arguments_optim& x) = 0;

  virtual void H(arguments_optim& x) = 0;

  virtual void outcomes(arguments_optim& x) = 0;

};

#include "estimators/efa/ml_efa.h"
#include "estimators/efa/uls_efa.h"
#include "estimators/efa/dwls_lt.h"

#include "estimators/rotation/cf.h"
#include "estimators/rotation/oblimin.h"
#include "estimators/rotation/geomin.h"
#include "estimators/rotation/varimax.h"
#include "estimators/rotation/varimin.h"
#include "estimators/rotation/target.h"
#include "estimators/rotation/xtarget.h"
#include "estimators/rotation/lclf.h"

#include "estimators/lca/lca.h"
#include "estimators/lca/bayesconst2.h"
#include "estimators/lca/bayesconst3.h"
#include "estimators/lca/bayesconst1.h"
#include "estimators/cfa/logdetmat.h"

#include "estimators/lreg/lreg.h"

#include "estimators/cfa/cfa_dwls.h"
#include "estimators/cfa/cfa_ml.h"

#include "estimators/correlation/polycor.h"

using EstimatorFactory = std::function<estimators*(const Rcpp::List&)>;

static const std::unordered_map<std::string, EstimatorFactory> estimator_factories = {
  { "ml_efa",                      choose_ml_efa                    },
  { "uls_efa",                     choose_uls_efa                   },
  { "dwls_lt",                     choose_dwls_lt                   },
  { "cfa_dwls",                    choose_cfa_dwls                  },
  { "cfa_ml",                      choose_cfa_ml                    },
  { "cf",                          choose_cf                        },
  { "oblimin",                     choose_oblimin                   },
  { "geomin",                      choose_geomin                    },
  { "varimax",                     choose_varimax                   },
  { "varimin",                     choose_varimin                   },
  { "target",                      choose_target                    },
  { "xtarget",                     choose_xtarget                   },
  { "lclf",                        choose_lclf                      },
  { "lca",                         choose_lca                       },
  { "bayesconst1",                 choose_bayesconst1               },
  { "bayesconst2",                 choose_bayesconst2               },
  { "bayesconst3",                 choose_bayesconst3               },
  { "logdetmat",                   choose_logdetmat                 },
  { "lreg",                        choose_lreg                      },
  { "polycor",                     choose_polycor                   }
};

estimators* choose_estimator(const Rcpp::List& estimator_setup) {
  const std::string name = Rcpp::as<std::string>(estimator_setup["estimator"]);
  auto it = estimator_factories.find(name);
  if (it == estimator_factories.end()) {
    Rcpp::stop("Unknown estimator: " + name);
  }
  return it->second(estimator_setup);
}

// Product Estimator

class product_estimator {

public:

  void param(arguments_optim& x, std::vector<estimators*>& xestimators) {

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->param(x);

    }

  }

  void F(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.f = 0;

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->F(x);

    }

  }

  void G(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.grad.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->G(x);

    }

  }

  void dG(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.dgrad.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->dG(x);

    }

  }

  void H(arguments_optim& x, std::vector<estimators*>& xestimators) {

    int ntrans = x.transparameters.n_elem;
    x.hess.set_size(ntrans, ntrans);
    x.hess.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->H(x);

    }

  }

  void E(arguments_optim& x, std::vector<estimators*>& xestimators) {

    // Update the log-likelihood and the estimated posterior using the
    // parameter estimates:

    // x.posterior.zeros();
    x.loglik = 0.00;

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->E(x);

    }

  }

  void M(arguments_optim& x, std::vector<estimators*>& xestimators) {

    // // Update the parameter estimates using the estimated posterior:
    //
    // for(int i=0; i < x.nestimators; ++i) {
    //
    //   xestimators[i]->posterior = x.posterior;
    //   arma::uvec indices = xestimators[i]->indices[0];
    //   xestimators[i]->transparameters = x.transparameters.elem(indices);
    //   xestimators[i]->M();
    //
    //   x.transparameters.elem(indices) = xestimators[i]->transparameters;
    //
    // }

  }

  void outcomes(arguments_optim& x, std::vector<estimators*>& xestimators) {

    std::get<0>(x.outputs_estimator).resize(x.nestimators);
    std::get<1>(x.outputs_estimator).resize(x.nestimators);
    std::get<2>(x.outputs_estimator).resize(x.nestimators);
    std::get<3>(x.outputs_estimator).resize(x.nestimators);
    std::get<4>(x.outputs_estimator).resize(x.nestimators);
    std::get<5>(x.outputs_estimator).resize(x.nestimators);

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->outcomes(x);

      std::get<0>(x.outputs_estimator)[i] = xestimators[i]->doubles;
      std::get<1>(x.outputs_estimator)[i] = xestimators[i]->vectors;
      std::get<2>(x.outputs_estimator)[i] = xestimators[i]->matrices;
      std::get<3>(x.outputs_estimator)[i] = xestimators[i]->cubes;
      std::get<4>(x.outputs_estimator)[i] = xestimators[i]->list_vectors;
      std::get<5>(x.outputs_estimator)[i] = xestimators[i]->list_matrices;

    }

  }

};

