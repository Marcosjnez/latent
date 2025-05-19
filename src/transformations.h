/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/05/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::mat jacob, grad, g;
  arma::uvec indices, target_indices;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void transform() = 0;

  virtual void jacobian() = 0;

  virtual void d2jacobian() = 0;

  virtual void dconstraints() = 0;

  virtual void outcomes() = 0;

};

#include "transformations/identity.h"
#include "transformations/softmax.h"
#include "transformations/expo.h"

// Choose the transformation:

transformations* choose_transform(Rcpp::List trans_setup, transformations* xtrans) {

  transformations* trans;
  std::string transform = trans_setup["transform"];

  if(transform == "identity") {

    trans = choose_identity(trans_setup);

  } else if(transform == "softmaxt") {

    trans = choose_softmaxt(trans_setup);

  } else if(transform == "exp") {

    trans = choose_expo(trans_setup);

  } else {

    Rcpp::stop("Available transformations: \n identity, exp, and softmaxt");

  }

  return trans;

}

// Product transformation:

class product_transform {

public:

  void transform(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    for(int i=0; i < x.ntransforms; ++i) {

      arma::uvec indices = xtransformations[i]->indices;
      xtransformations[i]->parameters = x.parameters.elem(indices);
      xtransformations[i]->transform();
      arma::uvec target_indices = xtransformations[i]->target_indices;
      // Rprintf("target_indices:\n");
      // for (arma::uword i = 0; i < target_indices.n_elem; ++i) {
      //   Rprintf("%u ", target_indices[i]);
      // }
      // Rprintf("\n\n");
      //
      // Rprintf("xtransformations[i]->transparameters:\n");
      // for (arma::uword j = 0; j < xtransformations[i]->transparameters.n_elem; ++j) {
      //   Rprintf("%g ", xtransformations[i]->transparameters[j]);
      // }
      // Rprintf("\n\n");
      //
      // Rprintf("x.transparameters:\n");
      // for (arma::uword j = 0; j < x.transparameters.n_elem; ++j) {
      //   Rprintf("%g ", x.transparameters[j]);
      // }
      // Rprintf("\n\n");
      // Rf_error("118");
      x.transparameters.elem(target_indices) = xtransformations[i]->transparameters;
      // Rf_error("128");

    }

  }

  void jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=0; i < x.ntransforms; ++i) {

      arma::uvec target_indices = xtransformations[i]->target_indices;
      xtransformations[i]->grad = x.grad.elem(target_indices);

      // Rprintf("indices:\n");
      // for (arma::uword i = 0; i < indices.n_elem; ++i) {
      //   Rprintf("%u ", indices[i]);
      // }
      // Rprintf("\n\n");
      //
      // Rprintf("x.grad:\n");
      // for (arma::uword j = 0; j < x.grad.n_elem; ++j) {
      //   Rprintf("%g ", x.grad[j]);
      // }
      // Rprintf("\n\n");
      //
      // Rf_error("150");

      xtransformations[i]->jacobian();

      // Rf_error("154");

      arma::uvec indices = xtransformations[i]->indices;
      x.g.elem(indices) += xtransformations[i]->g;

      // Rprintf("x.g:\n");
      // for (arma::uword j = 0; j < x.g.n_elem; ++j) {
      //   Rprintf("%g ", x.g[j]);
      // }
      // Rprintf("\n\n");
      //
      // Rf_error("163");

    }

  }

  void d2jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    for(int i=0; i < x.ntransforms; ++i) {
    }

  }

  void dconstraints(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // for(int i=0; i < x.ntransforms; ++i) {
    //
    // }

  }

  void outcomes(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // x.outputs_transform[0].resize(x.ntransforms);
    // x.outputs_transform[1].resize(x.ntransforms);
    // x.outputs_transform[2].resize(x.ntransforms);
    // x.outputs_transform[3].resize(x.ntransforms);

    std::get<0>(x.outputs_transform).resize(x.ntransforms);
    std::get<1>(x.outputs_transform).resize(x.ntransforms);
    std::get<2>(x.outputs_transform).resize(x.ntransforms);
    std::get<3>(x.outputs_transform).resize(x.ntransforms);

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->outcomes();

      std::get<0>(x.outputs_transform)[i] = xtransformations[i]->doubles;
      std::get<1>(x.outputs_transform)[i] = xtransformations[i]->vectors;
      std::get<2>(x.outputs_transform)[i] = xtransformations[i]->matrices;
      std::get<3>(x.outputs_transform)[i] = xtransformations[i]->list_matrices;

    }

  }

};
