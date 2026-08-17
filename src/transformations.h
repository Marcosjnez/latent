/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 11/08/2026
 */

// Transformations

class transformations {

public:

  arma::uvec indices_in, indices_out;

  std::vector<double> doubles;
  std::vector<arma::vec> vectors;
  std::vector<arma::mat> matrices;
  std::vector<arma::cube> cubes;
  std::vector<std::vector<arma::vec>> list_vectors;
  std::vector<std::vector<arma::mat>> list_matrices;

  std::vector<std::string> names_doubles, names_vectors, names_matrices,
  names_cubes, names_list_vectors, names_list_matrices;

  arma::mat jacob;

  virtual void transform(arguments_optim& x) = 0;

  virtual void update_grad(arguments_optim& x) = 0;

  virtual void dtransform(arguments_optim& x) = 0;

  virtual void update_dgrad(arguments_optim& x) = 0;

  virtual void jacobian(arguments_optim& x) = 0;

  virtual void update_vcov(arguments_optim& x) {}

  virtual void dconstraints(arguments_optim& x) {}

  virtual void outcomes(arguments_optim& x) = 0;

};

#include "transformations/identity.h"
#include "transformations/softmax.h"
#include "transformations/exponential.h"
#include "transformations/logarithm.h"
#include "transformations/crossprod.h"
#include "transformations/column_space.h"
#include "transformations/factor_cor.h"
#include "transformations/matrix_inverse.h"
#include "transformations/XY.h"
#include "transformations/XYt.h"
#include "transformations/deltaparam.h"
#include "transformations/mvnormal.h"
#include "transformations/normal.h"
#include "transformations/multinomial.h"
#include "transformations/sum_vectors.h"
#include "transformations/sqrt_vector.h"
#include "transformations/pos_incrsng.h"

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
  { "XYt", choose_XYt },
  { "deltaparam", choose_deltaparam },
  { "mvnormal",   choose_mvnormal   },
  { "sum_vectors",   choose_sum_vectors   },
  { "sqrt_vector",   choose_sqrt_vector   },
  { "pos_incrsng",   choose_pos_incrsng   }
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

    // Reset the initial transparameters vector.
    // This will include fixed values and allow overriding parameters
    x.transparameters = x.transparameters_init;
    x.transparameters(x.transparam2param) = x.parameters;

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->transform(x);

    }

  }

  void update_grad(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Update the gradient after each parameter transformation:
    // Use the gradient of the transformed parameters (grad) to get the final gradient (g):
    // x.g = x.jacob.t() * x.grad;

    // The following sequence avoids this multiplication to reduce computing cost
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
    x.jacob.set_size(x.transparameters.n_elem, x.transparameters.n_elem);

    for(arma::uword i : x.idx_transforms) {
      // for(int i=0; i < x.ntransforms; ++i) {

        arma::uvec indices_in = xtransformations[i]->indices_in;
        arma::uvec indices_out = xtransformations[i]->indices_out;
        xtransformations[i]->jacobian(x);

        for(arma::uword k = 0L; k < indices_in.n_elem; ++k) {
          for(arma::uword j = 0L; j < indices_out.n_elem; ++j) {
            x.jacob(indices_out[j], indices_in[k]) = xtransformations[i]->jacob(j, k);
          }
        }

      }

  }

  void update_vcov(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    const arma::uword ntrans = x.transparameters.n_elem;

    // Initialize the covariance and Jacobian matrices

    x.vcov.zeros(ntrans, ntrans);
    x.jacob.zeros(ntrans, ntrans);

    // Indicates which transformed-parameter positions currently have a valid
    // variance-covariance row and column.
    arma::uvec available(ntrans, arma::fill::zeros);

    // Variance-covariance matrix of the untransformed parameters

    // x.transparam2param may not contain contiguous indices, so insert the
    // covariance matrix element by element.
    for(arma::uword j = 0L; j < x.transparam2param.n_elem; ++j) {

      arma::uword jj = x.transparam2param[j];
      available[jj] = 1L;

      for(arma::uword i = 0L; i < x.transparam2param.n_elem; ++i) {

        arma::uword ii = x.transparam2param[i];

        x.vcov(ii, jj) = x.v(i, j);

      }

    }

    // Update the variance-covariance matrix through the transformations

    // x.idx_transforms must be ordered so that the transformations producing an
    // input parameter are evaluated before transformations using that parameter.
    for(arma::uword t : x.idx_transforms) {

      arma::uvec indices_in = xtransformations[t]->indices_in;
      arma::uvec indices_out = xtransformations[t]->indices_out;

      // Local Jacobian

      xtransformations[t]->jacobian(x);

      const arma::mat& J = xtransformations[t]->jacob;

      if(J.n_rows != indices_out.n_elem ||
         J.n_cols != indices_in.n_elem) {
        Rcpp::stop(
          "The Jacobian dimensions of transformation " +
            std::to_string(t + 1L) +
            " do not match its input and output parameter indices."
        );
      }

      // Store the local Jacobian in the complete Jacobian matrix.
      for(arma::uword j = 0L; j < indices_in.n_elem; ++j) {

        for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {

          x.jacob(indices_out[i], indices_in[j]) = J(i, j);

        }

      }

      // Parameters whose covariance is already available

      // Remove the current output coordinates from the set of already available
      // parameters. Their old covariance, if any, must not be used when the
      // transformation replaces them.
      arma::uvec available_previous = available;

      for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {
        available_previous[indices_out[i]] = 0L;
      }

      arma::uvec indices_previous =
        arma::find(available_previous > 0L);

      // Covariance among the transformation inputs

      arma::mat vcov_in(indices_in.n_elem,
                        indices_in.n_elem,
                        arma::fill::zeros);

      for(arma::uword j = 0L; j < indices_in.n_elem; ++j) {

        for(arma::uword i = 0L; i < indices_in.n_elem; ++i) {

          vcov_in(i, j) =
            x.vcov(indices_in[i], indices_in[j]);

        }

      }

      // Covariance between inputs and all previous parameters

      arma::mat vcov_in_previous(indices_in.n_elem,
                                 indices_previous.n_elem,
                                 arma::fill::zeros);

      for(arma::uword j = 0L; j < indices_previous.n_elem; ++j) {

        for(arma::uword i = 0L; i < indices_in.n_elem; ++i) {

          vcov_in_previous(i, j) =
            x.vcov(indices_in[i], indices_previous[j]);

        }

      }

      // Propagate the covariance

      // Cov(output, previous) =
      //   J * Cov(input, previous)
      arma::mat vcov_out_previous =
        J * vcov_in_previous;

      // Cov(output, output) =
      //   J * Cov(input, input) * J'
      arma::mat vcov_out =
        J * vcov_in * J.t();

      // Remove small numerical asymmetries.
      vcov_out = 0.5*(vcov_out + vcov_out.t());

      // Insert output/previous cross-covariances

      for(arma::uword j = 0L; j < indices_previous.n_elem; ++j) {

        arma::uword previous_index = indices_previous[j];

        for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {

          arma::uword output_index = indices_out[i];
          double value = vcov_out_previous(i, j);

          x.vcov(output_index, previous_index) = value;
          x.vcov(previous_index, output_index) = value;

        }

      }

      //  Insert output/output covariance

      for(arma::uword j = 0L; j < indices_out.n_elem; ++j) {

        arma::uword jj = indices_out[j];

        for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {

          arma::uword ii = indices_out[i];

          x.vcov(ii, jj) = vcov_out(i, j);

        }

      }

      // Make the outputs available to subsequent transformations

      for(arma::uword i = 0L; i < indices_out.n_elem; ++i) {
        available[indices_out[i]] = 1L;
      }

    }

    // Standard errors

    arma::vec variances = arma::vec(x.vcov.diag());

    // Protect only against negligible negative values caused by floating-point
    // error. Substantively negative variances remain negative and therefore
    // produce NaN, making the problem visible.
    for(arma::uword i = 0L; i < variances.n_elem; ++i) {

      if(variances[i] < 0.00 && variances[i] > -1e-12) {
        variances[i] = 0.00;
      }

    }

    x.se = arma::sqrt(variances);

  }

  void dconstraints(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    // Fill-in the constraints:
    for(arma::uword i : x.idx_transforms) {
      // for(int i=0L; i < x.ntransforms ; ++i) {
      xtransformations[i]->dconstraints(x);
    }

  }

  void outcomes(arguments_optim& x, std::vector<transformations*>& xtransformations) {

    std::get<0>(x.outputs_transform).resize(x.ntransforms);
    std::get<1>(x.outputs_transform).resize(x.ntransforms);
    std::get<2>(x.outputs_transform).resize(x.ntransforms);
    std::get<3>(x.outputs_transform).resize(x.ntransforms);
    std::get<4>(x.outputs_transform).resize(x.ntransforms);
    std::get<5>(x.outputs_transform).resize(x.ntransforms);

    std::get<6>(x.outputs_transform).resize(x.ntransforms);
    std::get<7>(x.outputs_transform).resize(x.ntransforms);
    std::get<8>(x.outputs_transform).resize(x.ntransforms);
    std::get<9>(x.outputs_transform).resize(x.ntransforms);
    std::get<10>(x.outputs_transform).resize(x.ntransforms);
    std::get<11>(x.outputs_transform).resize(x.ntransforms);

    for(int i=0; i < x.ntransforms; ++i) {

      xtransformations[i]->outcomes(x);

      std::get<0>(x.outputs_transform)[i] = xtransformations[i]->doubles;
      std::get<1>(x.outputs_transform)[i] = xtransformations[i]->vectors;
      std::get<2>(x.outputs_transform)[i] = xtransformations[i]->matrices;
      std::get<3>(x.outputs_transform)[i] = xtransformations[i]->cubes;
      std::get<4>(x.outputs_transform)[i] = xtransformations[i]->list_vectors;
      std::get<5>(x.outputs_transform)[i] = xtransformations[i]->list_matrices;

      std::get<6>(x.outputs_transform)[i] = xtransformations[i]->names_doubles;
      std::get<7>(x.outputs_transform)[i] = xtransformations[i]->names_vectors;
      std::get<8>(x.outputs_transform)[i] = xtransformations[i]->names_matrices;
      std::get<9>(x.outputs_transform)[i] = xtransformations[i]->names_cubes;
      std::get<10>(x.outputs_transform)[i] = xtransformations[i]->names_list_vectors;
      std::get<11>(x.outputs_transform)[i] = xtransformations[i]->names_list_matrices;

    }

  }

};

