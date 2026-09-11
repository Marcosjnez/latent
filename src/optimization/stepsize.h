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
#include "step_update/trust.h"

// Step-selection strategies. Every update leaves x at the selected point.
class stepsize {

public:

  virtual ~stepsize() = default;
  virtual bool uses_trust_region() const { return false; }

  virtual void update(arguments_optim& x,
                      std::vector<transformations*>& xtransforms,
                      std::vector<manifolds*>& xmanifolds,
                      std::vector<estimators*>& xestimators) = 0;

};

class ARMIJO: public stepsize {

public:

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators) override {

    armijo(x, xtransforms, xmanifolds, xestimators);

  }

};

class WOLFE: public stepsize {

public:

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators) override {

    wolfe(x, xtransforms, xmanifolds, xestimators);

  }

};

class TRUST: public stepsize {

public:

  double rad = 1.0;
  bool att_bnd = false;

  bool uses_trust_region() const override { return true; }

  void update(arguments_optim& x,
              std::vector<transformations*>& xtransforms,
              std::vector<manifolds*>& xmanifolds,
              std::vector<estimators*>& xestimators) override {

    trust(x, xtransforms, xmanifolds, xestimators, rad, att_bnd);

  }

};

inline void validate_stepsize(const std::string& step, bool newton_optimizer) {

  if(step != "armijo" && step != "wolfe" && step != "trust" && step != "tcg") {
    Rf_error("control$step must be one of 'armijo', 'wolfe', or 'trust'.");
  }

  if(!newton_optimizer && (step == "trust" || step == "tcg")) {
    Rf_error("control$step = 'trust' requires opt = 'newton', "
             "or opt = 'em' with mopt = 'newton'.");
  }

}

inline stepsize* choose_stepsize(const std::string& step,
                                 bool newton_optimizer = false) {

  validate_stepsize(step, newton_optimizer);

  if(step == "armijo") return new ARMIJO();
  if(step == "wolfe") return new WOLFE();
  return new TRUST(); // Includes the legacy 'tcg' spelling.

}

#endif
