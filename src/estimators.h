/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/05/2025
 */

class estimators {

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

#include "estimators/lreg/lreg.h"

#include "estimators/cfa/cfa_dwls.h"
#include "estimators/cfa/cfa_ml.h"

// Choose the estimator:

estimators* choose_estimator(Rcpp::List estimator_setup, estimators* xestimator) {

  estimators* criterion;
  std::string estimator = estimator_setup["estimator"];

  if (estimator == "ml_efa") {

    criterion = choose_ml_efa(estimator_setup);

  } else if (estimator == "uls_efa") {

    criterion = choose_uls_efa(estimator_setup);

  } else if (estimator == "dwls_lt") {

    criterion = choose_dwls_lt(estimator_setup);

  } else if(estimator == "cfa_dwls") {

    criterion = choose_cfa_dwls(estimator_setup);

  } else if(estimator == "cfa_ml") {

    criterion = choose_cfa_ml(estimator_setup);

  } else if (estimator == "cf") {

    criterion = choose_cf(estimator_setup);

  } else if (estimator == "oblimin") {

    criterion = choose_oblimin(estimator_setup);

  } else if (estimator == "geomin") {

    criterion = choose_geomin(estimator_setup);

  } else if (estimator == "varimax") {

    criterion = choose_varimax(estimator_setup);

  } else if (estimator == "varimin") {

    criterion = choose_varimin(estimator_setup);

  } else if (estimator == "target") {

    criterion = choose_target(estimator_setup);

  } else if (estimator == "xtarget") {

    criterion = choose_xtarget(estimator_setup);

  } else if (estimator == "lclf") {

    criterion = choose_lclf(estimator_setup);

  } else if (estimator == "lca_multinomial") {

    criterion = choose_lca_multinomial(estimator_setup);

  } else if (estimator == "lca_gaussian") {

    criterion = choose_lca_gaussian(estimator_setup);

  } else if (estimator == "latentloglik_combination") {

    criterion = choose_latentloglik_combination(estimator_setup);

  } else if (estimator == "lreg") {

    criterion = choose_lreg(estimator_setup);

  } else {

    Rcpp::stop("Unkown estimator");

  }

  return criterion;

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
      std::get<3>(x.outputs_estimator)[i] = xestimators[i]->list_matrices;

    }

  }

};

