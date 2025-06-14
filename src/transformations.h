/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/06/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::vec grad, g;
  arma::mat jacob;
  arma::uvec indices, target_indices;
  arma::cube jacob2;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void transform() = 0;

  virtual void jacobian() = 0;

  virtual void d2jacobian() = 0;

  virtual void dconstraints() = 0;

  virtual void outcomes() = 0;

};

#include "transformations/identity.h"
#include "transformations/softmax.h"
#include "transformations/exponential.h"
#include "transformations/lgaussian.h"

using TransformFactory =
  std::function< transformations*(const Rcpp::List&) >;

static const std::unordered_map<std::string, TransformFactory> transform_factories = {
  { "identity",    choose_identity    },
  { "softmax",     choose_softmax     },
  { "exponential", choose_exponential },
  { "lgaussian",   choose_lgaussian   }
};

transformations* choose_transform(const Rcpp::List& trans_setup) {
  const std::string name = Rcpp::as<std::string>(trans_setup["transform"]);
  auto it = transform_factories.find(name);
  if (it == transform_factories.end()) {
    Rcpp::stop("Unknown transform: " + name);
  }
  return it->second(trans_setup);
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
      // Rf_error("100");

    }

  }

  void jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // Rf_error("102");
    // x.jacob.set_size(x.transparameters.n_elem, x.parameters.n_elem);
    // x.jacob.zeros();
    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=0; i < x.ntransforms; ++i) {

      arma::uvec target_indices = xtransformations[i]->target_indices;
      xtransformations[i]->grad = x.grad(target_indices);

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
      // Rf_error("124");

      xtransformations[i]->jacobian();

      // Rf_error("128");

      arma::uvec indices = xtransformations[i]->indices;
      // x.jacob(target_indices, indices) += xtransformations[i]->jacob;
      x.g(indices) += xtransformations[i]->g;
      // x.g(indices) += xtransformations[i]->jacob.t() * x.grad(target_indices);

      // Rprintf("x.g:\n");
      // for (arma::uword j = 0; j < x.g.n_elem; ++j) {
      //   Rprintf("%g ", x.g[j]);
      // }
      // Rprintf("\n\n");
      //
      // Rf_error("143");

    }

    // x.g = x.jacob.t() * x.grad;

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
    std::get<4>(x.outputs_transform).resize(x.ntransforms);

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->outcomes();

      std::get<0>(x.outputs_transform)[i] = xtransformations[i]->doubles;
      std::get<1>(x.outputs_transform)[i] = xtransformations[i]->vectors;
      std::get<2>(x.outputs_transform)[i] = xtransformations[i]->matrices;
      std::get<3>(x.outputs_transform)[i] = xtransformations[i]->cubes;
      std::get<4>(x.outputs_transform)[i] = xtransformations[i]->list_matrices;

    }

  }

};
