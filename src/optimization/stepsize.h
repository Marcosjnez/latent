/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/09/2026
 */

#ifndef LATENT_STEPSIZE_H
#define LATENT_STEPSIZE_H

#include <memory>

#include "step_update/armijo.h"
#include "step_update/wolfe.h"
#include "step_update/tcg.h"

// Step-size strategies. Line searches update the parameters; a trust-region
// subsolver only proposes a direction, which Newton subsequently evaluates.

class stepsize {

public:

  virtual ~stepsize() = default;

  virtual void update(arguments_optim& x,
                      std::vector<transformations*>& xtransforms,
                      std::vector<manifolds*>& xmanifolds,
                      std::vector<estimators*>& xestimators) = 0;

  virtual void update(arguments_optim&,
                      std::vector<transformations*>&,
                      std::vector<manifolds*>&,
                      std::vector<estimators*>&,
                      bool&, const arma::vec&, double) {

    Rf_error("A Newton trust-region step requires control$step = 'tcg'.");

  }

};

class ARMIJO: public stepsize {

public:

  using stepsize::update;

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators) override {

                armijo(x, xtransforms, xmanifolds, xestimators);

              }

};

class WOLFE: public stepsize {

public:

  using stepsize::update;

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators) override {

                wolfe(x, xtransforms, xmanifolds, xestimators);

              }

};

class TCG: public stepsize {

public:

  void update(arguments_optim&,
              std::vector<transformations*>&,
              std::vector<manifolds*>&,
              std::vector<estimators*>&) override {

                Rf_error("control$step = 'tcg' requires the Newton trust-region optimizer.");

              }

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators,
              bool& att_bnd, const arma::vec& c, double rad) override {

                tcg(x, xtransforms, xmanifolds, xestimators, att_bnd, c, rad);

              }

};

inline void validate_stepsize(const std::string& step, bool trust_region) {

  if(step != "armijo" && step != "wolfe" && step != "tcg") {
    Rf_error("control$step must be one of 'armijo', 'wolfe', or 'tcg'.");
  }

  if(trust_region && step != "tcg") {
    Rf_error("The Newton trust-region optimizer requires control$step = 'tcg'. "
               "Use 'armijo' or 'wolfe' with opt = 'grad' or 'lbfgs', "
               "or with the corresponding EM M-step optimizer.");
  }

  if(!trust_region && step == "tcg") {
    Rf_error("control$step = 'tcg' requires opt = 'newton', "
               "or opt = 'em' with mopt = 'newton'.");
  }

}

inline stepsize* choose_stepsize(const std::string& step,
                                 bool trust_region = false) {

  validate_stepsize(step, trust_region);

  if(step == "armijo") {
    return new ARMIJO();
  }

  if(step == "wolfe") {
    return new WOLFE();
  }

  return new TCG();

}

#endif
