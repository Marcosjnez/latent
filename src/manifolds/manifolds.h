/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2026
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

  arma::mat T;

  virtual void param(arguments_optim& x) = 0;

  virtual void proj(arguments_optim& x) = 0;

  virtual void hess(arguments_optim& x) = 0;

  virtual void retr(arguments_optim& x) = 0;

  virtual void tangent_basis(arguments_optim& x) = 0;

  virtual void dconstraints(arguments_optim& x) = 0;

  virtual void outcomes(arguments_optim& x) = 0;

};

inline arma::mat tangent_complement(const arma::vec& x) {

  const arma::uword n = x.n_elem;

  if(n < 2L) {
    return arma::mat(n, 0L);
  }

  const double xnorm = arma::norm(x, 2);

  if(!std::isfinite(xnorm) || xnorm <= 0.00) {
    Rcpp::stop("A tangent-space basis cannot be computed at a zero or non-finite vector.");
  }

  arma::vec u = x/xnorm;
  arma::uword pivot = arma::index_min(arma::abs(u));
  arma::vec e(n, arma::fill::zeros);
  e[pivot] = 1.00;

  arma::vec v = e-u;
  double vnorm2 = arma::dot(v, v);
  arma::mat Q = arma::eye(n, n)-2.00*(v*v.t())/vnorm2;

  arma::uvec keep = arma::regspace<arma::uvec>(0L, n-1L);
  keep.shed_row(pivot);

  return Q.cols(keep);

}

#include "euclidean.h"
#include "unit.h"
#include "orth.h"
#include "orthog.h"
#include "oblq.h"
#include "poblq.h"
#include "poblq_blocks.h"
#include "simplex.h"

// type alias for factory functions
using ManifoldFactory = std::function<manifolds*(const Rcpp::List&)>;

// dispatch table
static const std::unordered_map<std::string, ManifoldFactory> manifold_factories = {
  { "euclidean", choose_euclidean },
  { "unit",      choose_unit      },
  { "simplex",   choose_simplex   },
  { "orth",      choose_orth      },
  { "orthog",    choose_orthog    },
  { "oblq",      choose_oblq      },
  { "poblq",     choose_poblq     },
  { "poblq_blocks", choose_poblq_blocks }
};

manifolds* choose_manifold(const Rcpp::List& manifold_setup) {
  const std::string name = Rcpp::as<std::string>(manifold_setup["manifold"]);
  auto it = manifold_factories.find(name);
  if (it == manifold_factories.end()) {
    Rcpp::stop(
      "Unknown manifold ‘" + name +
        "’. Available: euclidean, unit, simplex, orth, orthog, oblq, poblq, poblq_blocks"
    );
  }
  return it->second(manifold_setup);
}

// Product Manifold:

class product_manifold {

public:

  arma::mat T;

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

  void tangent_basis(arguments_optim& x,
                     std::vector<manifolds*>& xmanifolds) {

    T.set_size(x.nparam, 0L);

    if(x.nmanifolds == 0L) {
      T = arma::eye(x.nparam, x.nparam);
      return;
    }

    for(int i=0L; i < x.nmanifolds; ++i) {

      xmanifolds[i]->tangent_basis(x);

      if(xmanifolds[i]->T.n_rows != static_cast<arma::uword>(x.nparam)) {
        Rcpp::stop("The tangent-space basis of manifold "+
                   std::to_string(i+1L)+
                   " does not have one row per free parameter.");
      }

      T = arma::join_rows(T, xmanifolds[i]->T);

    }

  }

  void dconstraints(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    // Fill-in the constraints:
    for(int i=0L; i < x.nmanifolds ; ++i) {
      xmanifolds[i]->dconstraints(x);
    }

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
