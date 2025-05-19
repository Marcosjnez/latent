/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/05/2025
 */

// Manifolds

class manifolds {

public:

  // Objects that will be exported from manifolds:

  // Random: Provide these in product_manifold:
  arma::vec parameters, dparameters;
  arma::mat g, dg, rg, dH;

  // Fixed: Provide these in choose_manifold:
  arma::uvec indices; // Which parameters are projected onto the manifold (used in product_manifold)
  std::size_t q = indices.n_elem;
  arma::mat PhiTarget;

  double ss;
  arma::vec dir;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void param() = 0;

  virtual void proj() = 0;

  virtual void hess() = 0;

  virtual void retr() = 0;

  virtual void dconstraints() = 0;

  virtual void outcomes() = 0;

};

#include "manifolds/euclidean.h"
#include "manifolds/unit.h"
#include "manifolds/orth.h"
#include "manifolds/oblq.h"
#include "manifolds/poblq.h"
#include "manifolds/simplex.h"

// Choose the manifold:

manifolds* choose_manifold(Rcpp::List manifold_setup, manifolds* xmanifold) {

  manifolds* manifold;
  std::string projection = manifold_setup["manifold"];

  if(projection == "euclidean") {

    manifold = choose_euclidean(manifold_setup);

  } else if(projection == "unit") {

    manifold = choose_unit(manifold_setup);

  } else if(projection == "simplex") {

    manifold = choose_simplex(manifold_setup);

  } else if(projection == "orth") {

    manifold = choose_orth(manifold_setup);

  } else if(projection == "oblq") {

    manifold = choose_oblq(manifold_setup);

  } else if(projection == "poblq") {

    manifold = choose_poblq(manifold_setup);

  } else {

    Rcpp::stop("Available manifolds: \n euclidean, unit, simplex, orth, oblq, and poblq");

  }

  return manifold;

}

// Product Manifold:

class product_manifold {

public:

  void param(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    for(int i=0; i < x.nmanifolds; ++i) {

      arma::uvec indices = xmanifolds[i]->indices;
      xmanifolds[i]->parameters = x.parameters.elem(indices);
      xmanifolds[i]->param();

    }

  }

  void proj(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    x.rg.set_size(x.parameters.n_elem); x.rg.zeros();

    for(int i=0; i < x.nmanifolds; ++i) {

      arma::uvec indices = xmanifolds[i]->indices;
      xmanifolds[i]->g = x.g.elem(indices);
      xmanifolds[i]->proj();
      x.rg.elem(indices) += arma::vectorise(xmanifolds[i]->rg);

    }

  }

  void hess(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    x.dH.set_size(x.parameters.n_elem); x.dH.zeros();

    for(int i=0; i < x.nmanifolds; ++i) {

      arma::uvec indices = xmanifolds[i]->indices;
      xmanifolds[i]->dparameters = x.dparameters.elem(indices);
      xmanifolds[i]->g = x.g.elem(indices);
      xmanifolds[i]->dg = x.dg.elem(indices);
      // xmanifolds[i]->rg = x.rg.elem(indices); // Maybe not necessary
      xmanifolds[i]->hess();
      x.dH.elem(indices) += arma::vectorise(xmanifolds[i]->dH);

    }

  }

  void retr(arguments_optim& x, std::vector<manifolds*>& xmanifolds) {

    for(int i=0; i < x.nmanifolds; ++i) {

      arma::uvec indices = xmanifolds[i]->indices;
      xmanifolds[i]->ss = x.ss;
      xmanifolds[i]->dir = x.dir(indices);
      // xmanifolds[i]->parameters = x.parameters.elem(indices); // Maybe not necessary
      xmanifolds[i]->retr();
      // x.dir(indices) = xmanifolds[i]->dir;
      x.parameters.elem(indices) = xmanifolds[i]->parameters;

    }

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

    for(int i=0; i < x.nmanifolds; ++i) {

      std::get<0>(x.outputs_manifold)[i] = xmanifolds[i]->doubles;
      std::get<1>(x.outputs_manifold)[i] = xmanifolds[i]->vectors;
      std::get<2>(x.outputs_manifold)[i] = xmanifolds[i]->matrices;
      std::get<3>(x.outputs_manifold)[i] = xmanifolds[i]->list_matrices;

    }

  }

};

