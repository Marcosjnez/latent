/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 20/08/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::vec grad_in, grad_out, g;
  arma::mat jacob, sum_djacob, hess_in, hess_out, h;
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
      xtransformations[i]->transparameters = x.transparameters.elem(indices_in);
      xtransformations[i]->transform();
      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      x.transparameters.elem(indices_out) = xtransformations[i]->transparameters;
    }

  }

  void jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // x.jacob.set_size(x.transparameters.n_elem, x.parameters.n_elem);
    // x.jacob.zeros();
    x.g.set_size(x.parameters.n_elem); x.g.zeros();

    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->grad_in = x.grad(indices_out);

      xtransformations[i]->jacobian();

      arma::uvec indices_in = xtransformations[i]->indices_in[0];
      x.grad(indices_in) += xtransformations[i]->grad_out;

    }

    x.g = x.grad(x.transparam2param);
    // x.g = x.jacob.t() * x.grad;

  }

  void d2jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.h.set_size(x.parameters.n_elem, x.parameters.n_elem); x.h.zeros();
    arma::mat jacob(x.transparameters.n_elem, x.transparameters.n_elem, arma::fill::eye);
    arma::mat sum_djacob(x.transparameters.n_elem, x.transparameters.n_elem, arma::fill::zeros);

    // Update the hessian using the chain rule:
    arma::mat hess = x.hess;
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->grad_in = x.grad(indices_out);
      // xtransformations[i]->hess_in = x.hess(indices_out, indices_out);

      // COmpute the jacobian and sum_djacob of each transformation:
      xtransformations[i]->d2jacobian();

      // Update the jacobian and sum_djacob:
      arma::uvec indices_in = xtransformations[i]->indices_in[0];
      jacob(indices_out, indices_in) += xtransformations[i]->jacob;
      sum_djacob(indices_in, indices_in) += xtransformations[i]->sum_djacob;

      x.hess = jacob.t() * x.hess * jacob + sum_djacob;
      // hess(indices_in, indices_in) = x.hess(indices_in, indices_in);
      // hess(indices_out, indices_in) = x.hess(indices_out, indices_in);
      // hess(indices_in, indices_out) = x.hess(indices_in, indices_out);
      // x.hess(indices_in, indices_in) = xtransformations[i]->hess_out;

      // After modifying the hessian, restore the jacobian and sum_djacob:
      // jacob(indices_out, indices_in) -= xtransformations[i]->jacob;
      // sum_djacob(indices_in, indices_in) -= xtransformations[i]->sum_djacob;
      jacob(indices_out, indices_in).zeros();
      sum_djacob(indices_in, indices_in).zeros();

    }

    x.h = x.hess(x.transparam2param, x.transparam2param);
    // x.hess = hess;

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
