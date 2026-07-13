/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 13/07/2026
 */

// Manifolds

class manifolds {

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

  virtual void proj(arguments_optim& x) = 0;

  virtual void hess(arguments_optim& x) = 0;

  virtual void retr(arguments_optim& x) = 0;

  virtual void dconstraints(arguments_optim& x) = 0;

  virtual void outcomes(arguments_optim& x) = 0;

};

#include "manifolds/euclidean.h"
#include "manifolds/unit.h"
#include "manifolds/orth.h"
#include "manifolds/oblq.h"
#include "manifolds/poblq.h"
#include "manifolds/simplex.h"

// type alias for factory functions
using ManifoldFactory = std::function<manifolds*(const Rcpp::List&)>;

// dispatch table
static const std::unordered_map<std::string, ManifoldFactory> manifold_factories = {
  { "euclidean", choose_euclidean },
  { "unit",      choose_unit      },
  { "simplex",   choose_simplex   },
  { "orth",      choose_orth      },
  { "oblq",      choose_oblq      },
  { "poblq",     choose_poblq     }
};

manifolds* choose_manifold(const Rcpp::List& manifold_setup) {
  const std::string name = Rcpp::as<std::string>(manifold_setup["manifold"]);
  auto it = manifold_factories.find(name);
  if (it == manifold_factories.end()) {
    Rcpp::stop(
      "Unknown manifold ‘" + name +
        "’. Available: euclidean, unit, simplex, orth, oblq, poblq"
    );
  }
  return it->second(manifold_setup);
}

// Product Manifold:

class product_manifold {

public:

  void param(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    for(int i=0; i < x.nmanifolds; ++i) {

      xmanifolds[i]->param(x);

    }

  }

  void proj(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    x.rg.zeros();

    for(int i=0; i < x.nmanifolds; ++i) {

      xmanifolds[i]->proj(x);

    }

  }

  void hess(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    x.dH.zeros();

    for(int i=0; i < x.nmanifolds; ++i) {

      xmanifolds[i]->hess(x);

    }

  }

  void retr(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    for(int i=0; i < x.nmanifolds; ++i) {

      xmanifolds[i]->retr(x);

    }

    x.transparameters(x.transparam2param) = x.parameters;

  }

  void dconstraints(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    // for(int i=0; i < x.nmanifolds; ++i) {
    //
    // }

  }

  void outcomes(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    std::get<0>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<1>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<2>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<3>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<4>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<5>(x.outputs_manifold).resize(x.nmanifolds);

    std::get<6>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<7>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<8>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<9>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<10>(x.outputs_manifold).resize(x.nmanifolds);
    std::get<11>(x.outputs_manifold).resize(x.nmanifolds);

    for(int i=0; i < x.nmanifolds; ++i) {

      xmanifolds[i]->outcomes(x);

      std::get<0>(x.outputs_manifold)[i] = xmanifolds[i]->doubles;
      std::get<1>(x.outputs_manifold)[i] = xmanifolds[i]->vectors;
      std::get<2>(x.outputs_manifold)[i] = xmanifolds[i]->matrices;
      std::get<3>(x.outputs_manifold)[i] = xmanifolds[i]->cubes;
      std::get<4>(x.outputs_manifold)[i] = xmanifolds[i]->list_vectors;
      std::get<5>(x.outputs_manifold)[i] = xmanifolds[i]->list_matrices;

      std::get<6>(x.outputs_manifold)[i] = xmanifolds[i]->names_doubles;
      std::get<7>(x.outputs_manifold)[i] = xmanifolds[i]->names_vectors;
      std::get<8>(x.outputs_manifold)[i] = xmanifolds[i]->names_matrices;
      std::get<9>(x.outputs_manifold)[i] = xmanifolds[i]->names_cubes;
      std::get<10>(x.outputs_manifold)[i] = xmanifolds[i]->names_list_vectors;
      std::get<11>(x.outputs_manifold)[i] = xmanifolds[i]->names_list_matrices;

    }

  }

};

