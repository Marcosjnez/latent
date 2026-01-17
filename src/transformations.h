/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 30/10/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::vec grad_in, grad_out, g;
  arma::mat jacob, sum_djacob, hess_in, hess_out, h;
  std::vector<arma::uvec> indices_in, indices_out;
  bool constraints;
  arma::mat freqs;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void transform(arguments_optim& x) = 0;

  virtual void update_grad(arguments_optim& x) = 0;

  virtual void dtransform(arguments_optim& x) = 0;

  virtual void update_dgrad(arguments_optim& x) = 0;

  virtual void jacobian(arguments_optim& x) = 0;

  virtual void update_vcov(arguments_optim& x) = 0;

  virtual void dconstraints(arguments_optim& x) = 0;

  virtual void outcomes(arguments_optim& x) = 0;

};

#include "transformations/identity.h"
#include "transformations/softmax.h"
#include "transformations/exponential.h"
#include "transformations/logarithm.h"
#include "transformations/normal.h"
#include "transformations/crossprod.h"
#include "transformations/multinomial.h"
#include "transformations/column_space.h"
#include "transformations/factor_cor.h"
#include "transformations/matrix_inverse.h"
#include "transformations/XY.h"
#include "transformations/XYt.h"

using TransformFactory =
  std::function< transformations*(const Rcpp::List&) >;

static const std::unordered_map<std::string, TransformFactory> transform_factories = {
  { "identity",     choose_identity     },
  { "softmax",     choose_softmax     },
  { "exponential", choose_exponential },
  { "logarithm", choose_logarithm },
  { "normal",   choose_normal   },
  { "crossprod", choose_crossprod },
  { "multinomial", choose_multinomial },
  { "column_space", choose_column_space },
  { "factor_cor", choose_factor_cor },
  { "matrix_inverse", choose_matrix_inverse },
  { "XY", choose_XY },
  { "XYt", choose_XYt }
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

      xtransformations[i]->transform(x);

    }

  }

  void update_grad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Update the gradient after each parameter transformation:
    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // x.g = x.jacob.t() * x.grad;
    // This sequence avoids this multiplication to reduce computing cost

    // x.g.zeros();
    // Iterate top-down:
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      xtransformations[i]->update_grad(x);

    }

    // Gradient:
    x.g = x.grad(x.transparam2param);

  }

  void dtransform(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.dtransparameters.zeros();
    x.dtransparameters(x.transparam2param) = x.dparameters;

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->dtransform(x);

    }

  }

  void update_dgrad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Iterate top-down:
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      xtransformations[i]->update_dgrad(x);

    }

    // Differential of the gradient:
    x.dg = x.dgrad(x.transparam2param);

  }

  void jacobian(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    x.transparameters(x.transparam2param) = x.parameters;

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->jacobian(x);

    }

  }

  void update_vcov(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // This assumes that the jacobian is already computed
    // Approximate the inverse if the hessian is not positive-definite:
    x.inv_h = arma::inv_sympd(x.h, arma::inv_opts::allow_approx);

    // Get the variance-covariance matrix of everything:
    x.vcov.set_size(x.transparameters.n_elem, x.transparameters.n_elem);
    x.vcov.zeros();
    x.vcov(x.transparam2param, x.transparam2param) = x.inv_h;

    for (int i=0; i < x.ntransforms; ++i) {

      // Rprintf("%d \n", i);
      const arma::uvec& indices_in  = xtransformations[i]->indices_in[0];
      const arma::uvec& indices_out = xtransformations[i]->indices_out[0];

      xtransformations[i]->update_vcov(x);

      const arma::mat& J = xtransformations[i]->jacob;
      x.vcov(indices_out, indices_out) = J * x.vcov(indices_in, indices_in) * J.t();

    }

    // Finally, get the standard errors:
    x.se = arma::sqrt(x.vcov.diag());

  }

  void dconstraints(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    int nconstr = 0L;

    // Count the number of sets of constrained parameters:
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      xtransformations[i]->dconstraints(x);
      if(xtransformations[i]->constraints) ++ nconstr;

    }

    // Compute the derivative of the constraints of each transformation:
    x.mat_dconstraints.set_size(x.transparameters.n_elem, nconstr);
    x.mat_dconstraints.zeros();
    int k=0L;
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->dconstraints(x);
      if(xtransformations[i]->constraints) {
        x.mat_dconstraints.submat(indices_out, arma::uvec{ static_cast<arma::uword>(k) }) = xtransformations[i]->dconstr;
        ++k;
      };

    }

  }

  void outcomes(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    std::get<0>(x.outputs_transform).resize(x.ntransforms);
    std::get<1>(x.outputs_transform).resize(x.ntransforms);
    std::get<2>(x.outputs_transform).resize(x.ntransforms);
    std::get<3>(x.outputs_transform).resize(x.ntransforms);
    std::get<4>(x.outputs_transform).resize(x.ntransforms);
    std::get<5>(x.outputs_transform).resize(x.ntransforms);

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->outcomes(x);

      std::get<0>(x.outputs_transform)[i] = xtransformations[i]->doubles;
      std::get<1>(x.outputs_transform)[i] = xtransformations[i]->vectors;
      std::get<2>(x.outputs_transform)[i] = xtransformations[i]->matrices;
      std::get<3>(x.outputs_transform)[i] = xtransformations[i]->cubes;
      std::get<4>(x.outputs_transform)[i] = xtransformations[i]->list_vectors;
      std::get<5>(x.outputs_transform)[i] = xtransformations[i]->list_matrices;

    }

  }

};

