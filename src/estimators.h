/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/06/2025
 */

arma::vec row_sum_ignore_nan(const arma::mat& X) {
  arma::vec result(X.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      double val = X(i, j);
      if (arma::is_finite(val)) {
        result(i) += val;
      }
    }
  }
  return result;
}

public:

  double f, loglik;
  arma::vec parameters, transparameters, dparameters;
  arma::mat g, dg;
  arma::mat grad, dgrad;
  arma::mat hessian, dparam_dS, modhessian;
  arma::mat posterior, latentloglik;
  arma::vec uniquenesses;
  // arma::mat lambda, phi, psi, model, residuals;
  int nhessian, nS;
  arma::vec n, latentpars, loglatentpars; // P(X = c) // classes_hat
  std::vector<arma::mat> conditionals; // conditionals_hat

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::mat>> list_matrices;

  arma::uvec indices;
  int q;
  bool evalEM;

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

#include "estimators/lca/lca_multinomial.h"
#include "estimators/lca/lca_gaussian.h"
#include "estimators/lca/latentloglik_combination.h"
#include "estimators/lca/lca.h"

#include "estimators/lreg/lreg.h"

#include "estimators/cfa/cfa_dwls.h"
#include "estimators/cfa/cfa_ml.h"

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
  { "lca_multinomial",             choose_lca_multinomial           },
  { "lca_gaussian",                choose_lca_gaussian              },
  { "latentloglik_combination",    choose_latentloglik_combination  },
  { "lca",                         choose_lca                       },
  { "lreg",                        choose_lreg                      }
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

    x.latentloglik.zeros();
    // x.latentpars = xestimators[1]->latentpars;

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->transparameters = x.transparameters.elem(indices);
      // xestimators[i]->latentpars = x.latentpars;
      xestimators[i]->latentloglik = x.latentloglik;

      xestimators[i]->param();

      // arma::vec v = xestimators[i]->latentpars;
      // for (arma::uword j = 0; j < v.n_elem; ++j) {
      //   Rprintf("%g \n", v[j]);
      // }

      x.latentloglik += xestimators[i]->latentloglik;
      x.latentpars = xestimators[i]->latentpars; // ADD INDICES
      x.loglatentpars = xestimators[i]->loglatentpars; // ADD INDICES

    }

  }

  void F(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.f = 0;

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->latentloglik = x.latentloglik;
      xestimators[i]->latentpars = x.latentpars; // ADD INDICES

      xestimators[i]->F();
      x.f += xestimators[i]->f;

    }

  }

  void G(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.grad.set_size(x.transparameters.n_elem); x.grad.zeros();

    // Rprintf("%zu ", x.grad.n_elem);
    // Rprintf("\n\n");

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->latentloglik = x.latentloglik;
      xestimators[i]->latentpars = x.latentpars; // ADD INDICES
      xestimators[i]->loglatentpars = x.loglatentpars; // ADD INDICES

      xestimators[i]->G();

      arma::uvec indices = xestimators[i]->indices;

      // Rprintf("Indices:\n");
      // for (arma::uword i = 0; i < indices.n_elem; ++i) {
      //   Rprintf("%u ", indices[i]);
      // }
      // Rprintf("\n\n");
      //
      // Rprintf("xestimators[i]->grad.n_elem:\n");
      // Rprintf("%zu ", xestimators[i]->grad.n_elem);
      // Rprintf("\n\n");
      //
      // Rprintf("x.grad.n_elem:\n");
      // Rprintf("%zu ", x.grad.n_elem);
      // Rprintf("\n\n");
      //
      // Rf_error("921");
      x.grad.elem(indices) = xestimators[i]->grad;
      // Rf_error("922");

      // for (arma::uword j = 0; j < x.grad.n_elem; ++j) {
      //   Rprintf("%grad \n", x.grad[j]);
      // }

    }

  }

  void dG(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.dg.set_size(x.transparameters.n_elem); x.dg.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->dparameters = x.dparameters.elem(indices);
      xestimators[i]->g = x.grad.elem(indices); // Maybe delete PROBLEMS
      xestimators[i]->dG();
      x.dg.elem(indices) += xestimators[i]->dg;

    }

  }

  void H(arguments_optim& x, std::vector<estimators*>& xestimators) {

    int nhessian = x.transparameters.n_elem;
    x.hessian.set_size(nhessian, nhessian); x.hessian.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      xestimators[i]->H();
      arma::uvec indices = xestimators[i]->indices;
      x.hessian(indices, indices) += xestimators[i]->hessian;

    }

  }

  void E(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.latentloglik.zeros();

    for(int i=0; i < x.nestimators; ++i) {

      // Rprintf("Posterior:\n");
      // Rprintf("%zu ", x.posterior.n_rows);
      // Rprintf("%zu ", x.posterior.n_cols);
      // Rprintf("\n\n");

      if(xestimators[i]->evalEM) {

        arma::mat freqs = x.posterior; // S x nclasses
        x.n = xestimators[0]->n;
        freqs.each_col() %= x.n;
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

        arma::uvec indices = xestimators[i]->indices;
        x.transparameters.elem(indices) = xestimators[i]->transparameters;
        x.latentloglik += xestimators[i]->latentloglik;

      }

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
    x.n = xestimators[0]->n;
    // for(int s=0; s < S; ++s) {
    //   x.posterior.row(s) = arma::trunc_exp(x.latentloglik.row(s) + x.latentpars.t()); // P(data |X = c) P(X = c)
    //   mid[s] = arma::accu(x.posterior.row(s)); // P(data)
    // }

    x.posterior = x.latentloglik;
    x.posterior.each_row() += x.loglatentpars.t();
    x.posterior = arma::trunc_exp(x.posterior);
    // mid = row_sum_ignore_nan(x.posterior);
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

    arma::vec logliks = x.n % arma::trunc_log(mid);
    x.loglik = -sum_ignore_na(logliks);
    x.loglik = -arma::accu(logliks);

  }

  void outcomes(arguments_optim& x, std::vector<estimators*>& xestimators) {

    x.modhessian.resize(x.nestimators);
    x.dparam_dS.resize(x.nestimators);
    // x.doubles.resize(x.nestimators);
    // x.vectors.resize(x.nestimators);
    // x.matrices.resize(x.nestimators);
    // x.list_matrices.resize(x.nestimators);
    std::get<0>(x.outputs_estimator).resize(x.nestimators);
    std::get<1>(x.outputs_estimator).resize(x.nestimators);
    std::get<2>(x.outputs_estimator).resize(x.nestimators);
    std::get<3>(x.outputs_estimator).resize(x.nestimators);
    std::get<4>(x.outputs_estimator).resize(x.nestimators);

    for(int i=0; i < x.nestimators; ++i) {

      arma::uvec indices = xestimators[i]->indices;
      xestimators[i]->transparameters = x.transparameters.elem(indices);
      xestimators[i]->outcomes();
      x.modhessian[i] = xestimators[i]->modhessian;
      x.dparam_dS[i] = xestimators[i]->dparam_dS;
      // x.doubles[i] = xestimators[i]->doubles;
      // x.vectors[i] = xestimators[i]->vectors;
      // x.matrices[i] = xestimators[i]->matrices;
      // x.list_matrices[i] = xestimators[i]->list_matrices;
      std::get<0>(x.outputs_estimator)[i] = xestimators[i]->doubles;
      std::get<1>(x.outputs_estimator)[i] = xestimators[i]->vectors;
      std::get<2>(x.outputs_estimator)[i] = xestimators[i]->matrices;
      std::get<3>(x.outputs_estimator)[i] = xestimators[i]->cubes;
      std::get<4>(x.outputs_estimator)[i] = xestimators[i]->list_matrices;

    }

  }

};

