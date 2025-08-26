/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 25/08/2025
 */

// Transformations

class transformations {

public:

  arma::vec parameters, dparameters, transparameters, dconstr;
  arma::vec grad_in, grad_out, g;
  arma::mat jacob, sum_djacob, hess_in, hess_out, h;
  std::vector<arma::uvec> indices_in, indices_out;
  arma::cube jacob2;
  bool constraints;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  virtual void transform() = 0;

  virtual void update_grad() = 0;

  virtual void update_hess() = 0;

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

  void update_grad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Update the gradient after each parameter transformation:
    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // x.g = x.jacob.t() * x.grad;
    // This sequence avoids this multiplication to reduce computing cost

    // Iterate up-down:
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->grad_in = x.grad(indices_out);

      xtransformations[i]->update_grad();

      arma::uvec indices_in = xtransformations[i]->indices_in[0];
      x.grad(indices_in) += xtransformations[i]->grad_out;

    }

    // Gradient:
    x.g = x.grad(x.transparam2param);

  }

  void update_hess(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Faster version of hessian updating:
    for (int i = x.ntransforms - 1; i >= 0; --i) {

      const arma::uvec& indices_in  = xtransformations[i]->indices_in[0];
      const arma::uvec& indices_out = xtransformations[i]->indices_out[0];

      // Compute the jacobian and sum_djacob of each transformation:
      xtransformations[i]->update_hess();

      // Update the hessian: H += H U + U^T H + U^T H U + S
      const arma::mat& J = xtransformations[i]->jacob;
      const arma::mat& S = xtransformations[i]->sum_djacob;

      // HU = H * U:
      arma::mat HU = x.hess.cols(indices_out) * J;

      // H += H U + U^T H:
      x.hess.cols(indices_in) += HU;
      x.hess.rows(indices_in) += HU.t();

      // H += U^T H U + S:
      x.hess(indices_in, indices_in) += J.t() * x.hess(indices_out, indices_out) * J;
      x.hess(indices_in, indices_in) += S;
      x.hess = arma::symmatu(x.hess); // Ensure symmetry

    }

    // Get the hessian and inverse hessian:
    // Warning: if the hessian is not positive-definite, an approximation is done
    x.h = x.hess(x.transparam2param, x.transparam2param);
    x.inv_h = arma::inv_sympd(x.h, arma::inv_opts::allow_approx);

    // Get the variance-covariance matrix of everything:
    x.vcov.set_size(x.transparameters.n_elem, x.transparameters.n_elem);
    x.vcov.zeros();
    x.vcov(x.transparam2param, x.transparam2param) = x.inv_h;

    for (int i=0; i < x.ntransforms; ++i) {

      const arma::uvec& indices_in  = xtransformations[i]->indices_in[0];
      const arma::uvec& indices_out = xtransformations[i]->indices_out[0];

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

      xtransformations[i]->dconstraints();
      if(xtransformations[i]->constraints) ++ nconstr;

    }

    // Compute the derivative of the constraints of each transformation:
    x.mat_dconstraints.set_size(x.transparameters.n_elem, nconstr);
    x.mat_dconstraints.zeros();
    int k=0L;
    for(int i=x.ntransforms-1L; i > -1L ; --i) {

      arma::uvec indices_out = xtransformations[i]->indices_out[0];
      xtransformations[i]->dconstraints();
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


// arma::mat jacob(x.transparameters.n_elem, x.transparameters.n_elem, arma::fill::eye);
// arma::mat sum_djacob(x.transparameters.n_elem, x.transparameters.n_elem, arma::fill::zeros);

// Update the hessian using the chain rule:
// for(int i=x.ntransforms-1L; i > -1L ; --i) {
//
//   // Indices of input parameters:
//   arma::uvec indices_in = xtransformations[i]->indices_in[0];
//
//   // Indices of output parameters:
//   arma::uvec indices_out = xtransformations[i]->indices_out[0];
//
//   // Compute the jacobian and sum_djacob of each transformation:
//   xtransformations[i]->update_hess();
//
//   // Update the overall jacobian and sum_djacob:
//   jacob(indices_out, indices_in) += xtransformations[i]->jacob;
//   sum_djacob(indices_in, indices_in) += xtransformations[i]->sum_djacob;
//
//   x.hess = jacob.t() * x.hess * jacob + sum_djacob;
//
//   // After modifying the hessian, restore the jacobian and sum_djacob
//   // of the completed tranformations so they do not affect next computations:
//   jacob(indices_out, indices_in).zeros();
//   sum_djacob(indices_in, indices_in).zeros();
//
// }
