/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
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

  virtual void param() = 0;

  virtual void proj() = 0;

  virtual void hess() = 0;

  virtual void retr() = 0;

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

    euclidean* mymanifold = new euclidean();

    arma::uvec indices = manifold_setup["indices"];
    mymanifold->indices = indices;

    manifold = mymanifold;

  } else if(projection == "unit") {

    unit* mymanifold = new unit();

    // Provide these:
    arma::uvec indices = manifold_setup["indices"];
    mymanifold->indices = indices;

    manifold = mymanifold;

  } else if(projection == "simplex") {

    simplex* mymanifold = new simplex();

    // Provide these:
    arma::uvec indices = manifold_setup["indices"];
    mymanifold->indices = indices;

    manifold = mymanifold;

  } else if(projection == "orth") {

    orth* mymanifold = new orth();

    // Provide these:
    arma::uvec indices = manifold_setup["indices"];
    std::size_t q = manifold_setup["q"];
    mymanifold->indices = indices;
    mymanifold->q = q;

    manifold = mymanifold;

  } else if(projection == "oblq") {

    oblq* mymanifold = new oblq();

    // Provide these:
    arma::uvec indices = manifold_setup["indices"];
    std::size_t q = manifold_setup["q"];
    mymanifold->indices = indices;
    mymanifold->q = q;

    manifold = mymanifold;

  } else if(projection == "poblq") {

    poblq* mymanifold = new poblq();

    // Provide these:
    arma::uvec indices = manifold_setup["indices"];
    std::size_t q = manifold_setup["q"];
    // arma::uvec oblq_indices = manifold_setup["oblq_indices"];
    arma::mat PhiTarget = manifold_setup["PhiTarget"];
    mymanifold->indices = indices;
    mymanifold->q = q;
    mymanifold->PhiTarget = PhiTarget;

    PhiTarget.diag() += 10;
    arma::uvec oblq_indices = arma::find(PhiTarget == 1);
    mymanifold->oblq_indices = oblq_indices;

    manifold = mymanifold;

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

};

