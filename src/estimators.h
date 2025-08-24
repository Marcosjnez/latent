/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2025
 */

class estimators {

public:

  double f, loglik;
  arma::vec parameters, transparameters, dparameters;
  arma::mat g, dg;
  arma::mat grad, dgrad;
  arma::mat hess, dparam_dS, modhessian;
  arma::mat posterior, latentloglik;
  arma::vec uniquenesses;
  // arma::mat lambda, phi, psi, model, residuals;
  int nhessian, nS, p, q;
  arma::vec weights, latentpars, loglatentpars; // P(X = c) // classes_hat
  std::vector<arma::mat> conditionals; // conditionals_hat

  std::vector<arma::uvec> indices;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void param() = 0;

  virtual void F() = 0;

  virtual void G() = 0;

  virtual void dG() = 0;

  virtual void E() = 0;

  virtual void M() = 0;

  virtual void H() = 0;

  virtual void outcomes() = 0;

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
#include "estimators/lca/bayesconst1.h"
#include "estimators/lca/bayesconst2.h"
#include "estimators/lca/bayesconst3.h"

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

    // x.latentloglik.zeros();
    // x.latentpars = xestimators[1]->latentpars;

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices[0];
      xestimators[i]->transparameters = x.transparameters.elem(indices);
      // xestimators[i]->latentpars = x.latentpars;
      // xestimators[i]->latentloglik = x.latentloglik;

      xestimators[i]->param();

      // arma::vec v = xestimators[i]->latentpars;
      // for (arma::uword j = 0; j < v.n_elem; ++j) {
      //   Rprintf("%g \n", v[j]);
      // }

      // x.latentloglik += xestimators[i]->latentloglik;
      // x.latentpars = xestimators[i]->latentpars; // ADD INDICES
      // x.loglatentpars = xestimators[i]->loglatentpars; // ADD INDICES

    }

  }

  void F(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.f = 0;

    for(int i=0; i < x.nestimators; ++i) {

      // xestimators[i]->latentloglik = x.latentloglik;
      // xestimators[i]->latentpars = x.latentpars; // ADD INDICES

      xestimators[i]->F();
      x.f += xestimators[i]->f;

    }

  }

  void G(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.grad.set_size(x.transparameters.n_elem); x.grad.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->G();

      arma::uvec indices = xestimators[i]->indices[0];
      x.grad.elem(indices) += xestimators[i]->grad;

    }

  }

  void dG(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.dg.set_size(x.transparameters.n_elem); x.dg.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices[0];
      xestimators[i]->dparameters = x.dparameters.elem(indices);
      xestimators[i]->g = x.grad.elem(indices); // Maybe delete PROBLEMS
      xestimators[i]->dG();
      x.dg.elem(indices) += xestimators[i]->dg;

    }

  }

  void H(arguments_optim& x, std::vector<estimators*>& xestimators) {

    int nhessian = x.transparameters.n_elem;
    x.hess.set_size(nhessian, nhessian); x.hess.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->H();
      arma::uvec indices = xestimators[i]->indices[0];
      x.hess(indices, indices) += xestimators[i]->hess;

    }

  }

  void E(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.latentloglik.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      // Rprintf("Posterior:\n");
      // Rprintf("%zu ", x.posterior.n_rows);
      // Rprintf("%zu ", x.posterior.n_cols);
      // Rprintf("\n\n");

        arma::mat freqs = x.posterior; // S x nclasses
        x.weights = xestimators[0]->weights;
        freqs.each_col() %= x.weights;
        arma::vec post_total = arma::sum(freqs, 0).t();
        x.latentpars = post_total / arma::accu(post_total);
        x.loglatentpars = trunc_log(x.latentpars);

        xestimators[i]->latentpars = x.latentpars;
        xestimators[i]->loglatentpars = x.loglatentpars;
        xestimators[i]->posterior = x.posterior;
        xestimators[i]->E();

        // Rprintf("Indices:\n");
        // for (arma::uword i = 0; i < indices.n_elem; ++i) {
        //   Rprintf("%u ", indices[i]);
        // }
        // Rprintf("\n\n");
        //
        // Rprintf("xestimators[i]->transparameters.n_elem:\n");
        // Rprintf("%zu ", xestimators[i]->transparameters.n_elem);
        // Rprintf("\n\n");
        //
        // Rf_error("812");

        arma::uvec indices = xestimators[i]->indices[0];
        x.transparameters.elem(indices) = xestimators[i]->transparameters;
        x.latentloglik += xestimators[i]->latentloglik;

    }

  }

  void M(arguments_optim& x, std::vector<estimators*>& xestimators) {

    // x.loglik = 0.00;
    //
    // for(int i=0; i < x.nestimators; ++i) {
    //
    //   xestimators[i]->latentloglik = x.latentloglik;
    //
    //   xestimators[i]->M();
    //   x.loglik += xestimators[i]->loglik;
    //   x.posterior = xestimators[i]->posterior;
    //
    // }
    //
    // x.loglik /= x.nestimators;

    // For each response pattern, compute its joint and marginal probabilities:
    int S = x.latentloglik.n_rows;
    arma::vec mid(S);
    x.latentpars = xestimators[0]->latentpars;
    x.loglatentpars = xestimators[0]->loglatentpars;
    x.weights = xestimators[0]->weights;
    // for(int s=0; s < S; ++s) {
    //   x.posterior.row(s) = arma::trunc_exp(x.latentloglik.row(s) + x.latentpars.t()); // P(data |X = c) P(X = c)
    //   mid[s] = arma::accu(x.posterior.row(s)); // P(data)
    // }

    x.posterior = x.latentloglik;
    x.posterior.each_row() += x.loglatentpars.t();
    x.posterior = arma::trunc_exp(x.posterior);
    mid = arma::sum(x.posterior, 1);

    x.posterior.each_col() /= mid; // P(X = c | data) = P(data | X = c) P(X = c) / P(data)

    // Rprintf("n:\n");
    // Rprintf("%zu ", x.n.n_elem);
    // Rprintf("%zu ", x.n.n_cols);
    // Rprintf("\n\n");
    //
    // Rprintf("mid:\n");
    // Rprintf("%zu ", mid.n_rows);
    // Rprintf("%zu ", mid.n_cols);
    // Rprintf("\n\n");

    arma::vec logliks = x.weights % arma::trunc_log(mid);
    x.loglik = -arma::accu(logliks);

  }

  void outcomes(arguments_optim& x, std::vector<estimators*>& xestimators) {

    std::get<0>(x.outputs_estimator).resize(x.nestimators);
    std::get<1>(x.outputs_estimator).resize(x.nestimators);
    std::get<2>(x.outputs_estimator).resize(x.nestimators);
    std::get<3>(x.outputs_estimator).resize(x.nestimators);
    std::get<4>(x.outputs_estimator).resize(x.nestimators);
    std::get<5>(x.outputs_estimator).resize(x.nestimators);

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->outcomes();

      std::get<0>(x.outputs_estimator)[i] = xestimators[i]->doubles;
      std::get<1>(x.outputs_estimator)[i] = xestimators[i]->vectors;
      std::get<2>(x.outputs_estimator)[i] = xestimators[i]->matrices;
      std::get<3>(x.outputs_estimator)[i] = xestimators[i]->cubes;
      std::get<4>(x.outputs_estimator)[i] = xestimators[i]->list_vectors;
      std::get<5>(x.outputs_estimator)[i] = xestimators[i]->list_matrices;

    }

  }

};

