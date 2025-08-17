/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 17/08/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::vec grad, g;
  arma::mat jacob, h;
  std::vector<arma::uvec> indices_in, indices_out;
  arma::cube jacob2;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
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
#include "transformations/normal.h"
#include "transformations/crossprod.h"
#include "transformations/multinomial.h"

using TransformFactory =
  std::function< transformations*(const Rcpp::List&) >;

static const std::unordered_map<std::string, TransformFactory> transform_factories = {
  { "identity",     choose_identity     },
  { "softmax",     choose_softmax     },
  { "exponential", choose_exponential },
  { "normal",   choose_normal   },
  { "crossprod", choose_crossprod },
  { "multinomial", choose_multinomial }
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

    x.transparameters(x.transparam2param) = x.parameters;

    for(int i=0; i < x.ntransforms; ++i) {

      arma::uvec indices_in = xtransformations[i]->indices_in[0];
      // Rprintf("76 \n");
      xtransformations[i]->transparameters = x.transparameters.elem(indices_in);
      // Rprintf("78\n");
      xtransformations[i]->transform();
      // Rprintf("80\n");
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
      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      // Rprintf("100 \n");
      x.transparameters.elem(indices_out) = xtransformations[i]->transparameters;
      // Rprintf("102 \n");
      // Rf_error("100");
    }

  }

  void jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // x.jacob.set_size(x.transparameters.n_elem, x.parameters.n_elem);
    // x.jacob.zeros();
    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->grad = x.grad(indices_out);

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

      arma::uvec indices_in = xtransformations[i]->indices_in[0];
      // Rprintf("Size: %zu\n", indices_in.n_elem);
      // Rprintf("Size: %zu\n", xtransformations[i]->grad.n_elem);
      x.grad(indices_in) += xtransformations[i]->grad;

      // Rprintf("x.g:\n");
      // for (arma::uword j = 0; j < x.g.n_elem; ++j) {
      //   Rprintf("%g ", x.g[j]);
      // }
      // Rprintf("\n\n");
      //
      // Rf_error("143");

    }

    x.g = x.grad(x.transparam2param);
    // x.g = x.jacob.t() * x.grad;

  }

  void d2jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.h.set_size(x.parameters.n_elem, x.parameters.n_elem); x.g.zeros();

  }

  void dconstraints(arguments_optim& x, std::vector<transformations*>& xtransformations) {

  }

  void outcomes(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    std::get<0>(x.outputs_transform).resize(x.ntransforms);
    std::get<1>(x.outputs_transform).resize(x.ntransforms);
    std::get<2>(x.outputs_transform).resize(x.ntransforms);
    std::get<3>(x.outputs_transform).resize(x.ntransforms);
    std::get<4>(x.outputs_transform).resize(x.ntransforms);
    std::get<5>(x.outputs_transform).resize(x.ntransforms);

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->outcomes();

      std::get<0>(x.outputs_transform)[i] = xtransformations[i]->doubles;
      std::get<1>(x.outputs_transform)[i] = xtransformations[i]->vectors;
      std::get<2>(x.outputs_transform)[i] = xtransformations[i]->matrices;
      std::get<3>(x.outputs_transform)[i] = xtransformations[i]->cubes;
      std::get<4>(x.outputs_transform)[i] = xtransformations[i]->list_vectors;
      std::get<5>(x.outputs_transform)[i] = xtransformations[i]->list_matrices;

    }

  }

};
