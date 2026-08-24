/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 24/08/2026
 */

class estimators {

public:

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  std::vector<std::string> names_doubles, names_vectors, names_matrices,
  names_cubes, names_list_vectors, names_list_matrices;

  virtual void param(arguments_optim& x) = 0;
  virtual void F(arguments_optim& x) = 0;
  virtual void G(arguments_optim& x) = 0;
  virtual void dG(arguments_optim& x) = 0;
  virtual void E(arguments_optim& x) {}

  virtual void observed_F(arguments_optim& x) {
    F(x);
  }

  virtual void observed_G(arguments_optim& x) {
    G(x);
  }

  virtual void outcomes(arguments_optim& x) = 0;

};

#include "rotation/cf.h"
#include "rotation/oblimin.h"
#include "rotation/geomin.h"
#include "rotation/varimax.h"
#include "rotation/varimin.h"
#include "rotation/target.h"
#include "rotation/xtarget.h"
#include "rotation/lclf.h"

#include "lca/lca.h"
#include "lca/lcaEM.h"
#include "lca/bayesconst1.h"
#include "lca/bayesconst2.h"
#include "lca/bayesconst3.h"
#include "lca/bayesconst4.h"
#include "cfa/logdetmat.h"
#include "cfa/logdetR.h"

#include "lreg/lreg.h"

#include "cfa/cfa_dwls.h"
#include "cfa/cfa_dwls_poly.h"
#include "cfa/cfa_ml.h"
#include "cfa/cfa_fml.h"

#include "loglik/gaussian_loglik.h"
#include "loglik/poisson_loglik.h"
#include "penalties/ridge.h"
#include "correlation/polycor.h"

using EstimatorFactory = std::function<estimators*(const Rcpp::List&)>;

static const std::unordered_map<std::string, EstimatorFactory> estimator_factories = {
  { "cfa_dwls",                    choose_cfa_dwls                  },
  { "cfa_dwls_poly",               choose_cfa_dwls_poly             },
  { "cfa_ml",                      choose_cfa_ml                    },
  { "cfa_fml",                     choose_cfa_fml                   },

  { "cf",                          choose_cf                        },
  { "oblimin",                     choose_oblimin                   },
  { "geomin",                      choose_geomin                    },
  { "varimax",                     choose_varimax                   },
  { "varimin",                     choose_varimin                   },
  { "target",                      choose_target                    },
  { "xtarget",                     choose_xtarget                   },
  { "lclf",                        choose_lclf                      },
  { "lca",                         choose_lca                       },
  { "lcaEM",                       choose_lcaEM                     },
  { "bayesconst1",                 choose_bayesconst1               },
  { "bayesconst2",                 choose_bayesconst2               },
  { "bayesconst3",                 choose_bayesconst3               },
  { "bayesconst4",                 choose_bayesconst4               },
  { "logdetmat",                   choose_logdetmat                 },
  { "logdetR",                     choose_logdetR                   },
  { "lreg",                        choose_lreg                      },
  { "gaussian_loglik",             choose_gaussian_loglik           },
  { "poisson_loglik",              choose_poisson_loglik            },
  { "ridge",                       choose_ridge                     },
  { "polycor",                     choose_polycor                   }
};

estimators* choose_estimator(const Rcpp::List& estimator_setup) {

  const std::string name = Rcpp::as<std::string>(estimator_setup["estimator"]);
  auto it = estimator_factories.find(name);

  if(it == estimator_factories.end()) {
    Rcpp::stop("Unknown estimator: " + name);
  }

  return it->second(estimator_setup);

}

class product_estimator {

public:

  void param(arguments_optim& x, std::vector<estimators*>& xestimators) {
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->param(x);
  }

  void F(arguments_optim& x, std::vector<estimators*>& xestimators) {
    x.f = 0;
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->F(x);
  }

  void G(arguments_optim& x, std::vector<estimators*>& xestimators) {
    x.grad.zeros();
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->G(x);
  }

  void dG(arguments_optim& x, std::vector<estimators*>& xestimators) {
    x.dgrad.zeros();
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->dG(x);
  }

  void E(arguments_optim& x, std::vector<estimators*>& xestimators) {
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->E(x);
  }

  void observed_F(arguments_optim& x,
                  std::vector<estimators*>& xestimators) {
    x.f = 0.0;
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->observed_F(x);
  }

  void observed_G(arguments_optim& x,
                  std::vector<estimators*>& xestimators) {
    x.grad.zeros();
    for(int i = 0; i < x.nestimators; ++i) xestimators[i]->observed_G(x);
  }

  void outcomes(arguments_optim& x,
                std::vector<estimators*>& xestimators) {

    std::get<0>(x.outputs_estimator).resize(x.nestimators);
    std::get<1>(x.outputs_estimator).resize(x.nestimators);
    std::get<2>(x.outputs_estimator).resize(x.nestimators);
    std::get<3>(x.outputs_estimator).resize(x.nestimators);
    std::get<4>(x.outputs_estimator).resize(x.nestimators);
    std::get<5>(x.outputs_estimator).resize(x.nestimators);
    std::get<6>(x.outputs_estimator).resize(x.nestimators);
    std::get<7>(x.outputs_estimator).resize(x.nestimators);
    std::get<8>(x.outputs_estimator).resize(x.nestimators);
    std::get<9>(x.outputs_estimator).resize(x.nestimators);
    std::get<10>(x.outputs_estimator).resize(x.nestimators);
    std::get<11>(x.outputs_estimator).resize(x.nestimators);

    for(int i = 0; i < x.nestimators; ++i) {

      xestimators[i]->outcomes(x);

      std::get<0>(x.outputs_estimator)[i] = xestimators[i]->doubles;
      std::get<1>(x.outputs_estimator)[i] = xestimators[i]->vectors;
      std::get<2>(x.outputs_estimator)[i] = xestimators[i]->matrices;
      std::get<3>(x.outputs_estimator)[i] = xestimators[i]->cubes;
      std::get<4>(x.outputs_estimator)[i] = xestimators[i]->list_vectors;
      std::get<5>(x.outputs_estimator)[i] = xestimators[i]->list_matrices;
      std::get<6>(x.outputs_estimator)[i] = xestimators[i]->names_doubles;
      std::get<7>(x.outputs_estimator)[i] = xestimators[i]->names_vectors;
      std::get<8>(x.outputs_estimator)[i] = xestimators[i]->names_matrices;
      std::get<9>(x.outputs_estimator)[i] = xestimators[i]->names_cubes;
      std::get<10>(x.outputs_estimator)[i] = xestimators[i]->names_list_vectors;
      std::get<11>(x.outputs_estimator)[i] = xestimators[i]->names_list_matrices;

    }

  }

};
