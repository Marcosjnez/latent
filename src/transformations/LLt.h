#ifndef LATENT_TRANSFORMATIONS_LLT_H
#define LATENT_TRANSFORMATIONS_LLT_H

#include <cmath>
#include <stdexcept>

// Include from transformations.h after the definition of transformations.
// choose_LLt setup: target (symmetric binary matrix, diagonal ones),
// indices_in / indices_out (one index vector each, zero-based),
// standardize (default true), output ("gram", default, or "factor").
//
// For i>j, target(i,j)=0 means (L*L.t())(i,j)=0. Its factor entry is
// computed recursively; free lower entries are unconstrained real numbers.
// Inputs follow R column-major order, in one of the following layouts:
//   * all p*p positions of temp (unused positions have zero derivatives);
//   * the free lower triangle, including a log pivot on every diagonal;
//   * when standardizing, just the free strict lower triangle. This fixes
//     each raw pivot at one and removes the redundant row-scale parameters.
// The first two layouts implement exp(diag(temp)) even when standardizing.
// Outputs are all p*p entries, in column-major order. For standardized rows
// L is the row-normalized factor, so the Gram output is a correlation matrix.
// Positive pivots give positive definiteness in exact arithmetic. Forming
// the Gram matrix in double precision can lose rank near the boundary.
// Use distinct input and output blocks for latent's package-wide Jacobian.
// Mirrored output labels may be shared for the Gram matrix. The two halves
// of the triangular factor generally have different values. The outcome L
// is available after fitting; jacobian() builds the dense Jacobian on demand.

class LLt final : public transformations {
  using index = arma::uword;

  arma::mat target;
  bool standardize, factor_output;
  index p;
  arma::uvec positions, input_slots, overlap, constant_gram;
  arma::mat B, L, dB, dL, output_grad;
  arma::vec row_scale, scaled_norm, row_speed;

  // The denominator is ||B_i||. Splitting it keeps the norm calculation
  // finite even if squaring the raw entries would overflow.
  arma::rowvec divide_by_norm(const arma::rowvec& row, index i) const {
    if (row_scale[i] < 1.0) {
      const long double norm =
        static_cast<long double>(row_scale[i]) * scaled_norm[i];
      arma::rowvec result(row.n_elem);
      for (index k = 0; k < row.n_elem; ++k) {
        result[k] = static_cast<double>(row[k] / norm);
      }
      return result;
    }
    return (row / row_scale[i]) / scaled_norm[i];
  }

  arma::mat output_seed(const arma::vec& grad) const {
    arma::mat seed = matrix_output_gradient(grad, p);
    if (!factor_output) seed.elem(constant_gram).zeros();
    return seed;
  }

  // The tangent of the triangular completion is itself a triangular solve.
  // For a fixed position differentiate B_ij*B_jj + sum_{k<j} B_ik*B_jk = 0.
  void factor_differential(const arma::vec& seed, arma::mat& raw_dir,
                           arma::mat& factor_dir, arma::vec& speed) const {
    raw_dir.zeros(p, p);
    raw_dir.elem(positions) = seed;
    raw_dir.diag() %= B.diag();
    for (index i = 0; i < p; ++i) {
      for (index j = 0; j < i; ++j) {
        if (target(i,j) != 0.0) continue;
        long double sum = static_cast<long double>(B(i,j)) * raw_dir(j,j);
        for (index k = 0; k < j; ++k) {
          sum += static_cast<long double>(raw_dir(i,k)) * B(j,k)
               + static_cast<long double>(B(i,k)) * raw_dir(j,k);
        }
        raw_dir(i,j) = static_cast<double>(-sum / B(j,j));
      }
    }

    factor_dir = raw_dir;
    if (!standardize) return;
    speed.set_size(p);
    for (index i = 0; i < p; ++i) {
      const arma::rowvec scaled = divide_by_norm(raw_dir.row(i), i);
      speed[i] = arma::dot(L.row(i), scaled);  // d||B_i|| / ||B_i||
      factor_dir.row(i) = scaled - speed[i] * L.row(i);
    }
  }

  arma::mat output_differential(const arma::mat& factor_dir) const {
    if (factor_output) return factor_dir;
    arma::mat half = factor_dir * L.t();
    arma::mat result = half + half.t();
    result.elem(constant_gram).zeros();
    return result;
  }

  // Reverse the triangular solve, optionally differentiating the same sweep.
  // An implied entry feeds both earlier row prefixes and the earlier pivot.
  void reverse_completion(arma::mat& adj, arma::mat* dadj) const {
    for (index row = p; row > 0; --row) {
      const index i = row - 1;
      for (index col = i; col > 0; --col) {
        const index j = col - 1;
        if (target(i,j) != 0.0) continue;
        const long double pivot = B(j,j);
        const long double weight = -static_cast<long double>(adj(i,j)) / pivot;
        if (dadj) {
          const long double dweight =
            -((*dadj)(i,j) + weight * dB(j,j)) / pivot;
          for (index k = 0; k < j; ++k) {
            (*dadj)(i,k) += dweight * B(j,k) + weight * dB(j,k);
            (*dadj)(j,k) += dweight * B(i,k) + weight * dB(i,k);
          }
          (*dadj)(j,j) += dweight * B(i,j) + weight * dB(i,j);
        }
        for (index k = 0; k < j; ++k) {
          adj(i,k) += weight * B(j,k);
          adj(j,k) += weight * B(i,k);
        }
        adj(j,j) += weight * B(i,j);
      }
    }
    if (dadj) dadj->diag() = dadj->diag() % B.diag() + adj.diag() % dB.diag();
    adj.diag() %= B.diag();
  }

  // Pull back the supplied output gradient. With a directional seed, return
  // d(J.t()*g) = (dJ).t()*g + J.t()*dg, without constructing either Jacobian.
  arma::vec pullback(const arma::mat* dseed) const {
    arma::mat adj, dadj;
    if (factor_output) {
      adj = output_grad;
      if (dseed) dadj = *dseed;
    } else {
      const arma::mat sym = output_grad + output_grad.t();
      adj = sym * L;
      if (dseed) dadj = ((*dseed) + dseed->t()) * L + sym * dL;
    }
    if (standardize) {
      for (index i = 0; i < p; ++i) {
        const arma::rowvec g = adj.row(i);
        const double radial = arma::dot(L.row(i), g);
        adj.row(i) = divide_by_norm(g - radial * L.row(i), i);
        if (dseed) {
          const arma::rowvec dg = dadj.row(i);
          const double dradial = arma::dot(dL.row(i), g) + arma::dot(L.row(i), dg);
          dadj.row(i) = divide_by_norm(dg - dradial * L.row(i) - radial * dL.row(i), i)
                      - row_speed[i] * adj.row(i);
        }
      }
    }
    reverse_completion(adj, dseed ? &dadj : nullptr);
    return dseed ? arma::vec(dadj.elem(positions)) : arma::vec(adj.elem(positions));
  }

  void add_input_gradient(arma::vec& grad, const arma::vec& values) const {
    // Output assignment overwrites any overlapping input. Other terms that
    // use an input still accumulate, including repeated input labels.
    grad.elem(overlap).zeros();
    for (index k = 0; k < input_slots.n_elem; ++k) {
      grad[indices_in[input_slots[k]]] += values[k];
    }
  }

public:
  LLt(const arma::mat& mask, bool unit, bool factor,
      const arma::uvec& inputs, const arma::uvec& outputs)
    : target(mask), standardize(unit), factor_output(factor), p(mask.n_rows) {
    if (p == 0 || mask.n_cols != p || !mask.is_finite() ||
        !mask.is_symmetric(0.0) || arma::any(mask.diag() != 1.0) ||
        arma::accu((mask != 0.0) % (mask != 1.0)) != 0) {
      throw std::invalid_argument("LLt: target must be square, symmetric and binary with diagonal ones.");
    }
    if (outputs.n_elem != mask.n_elem) {
      throw std::invalid_argument("LLt: the output block must contain p*p indices.");
    }
    indices_in = inputs;
    set_matrix_output(outputs);
    overlap = arma::intersect(inputs, outputs);
    arma::mat active = arma::trimatl(mask);
    positions = arma::find(active);
    if (inputs.n_elem == mask.n_elem) {
      input_slots = positions;
    } else {
      if (standardize && inputs.n_elem == positions.n_elem - p) {
        active.diag().zeros();
        positions = arma::find(active);
      }
      if (inputs.n_elem != positions.n_elem) {
        throw std::invalid_argument("LLt: input block length does not match the target and input layout.");
      }
      input_slots.set_size(positions.n_elem);
      for (index k = 0; k < positions.n_elem; ++k) input_slots[k] = k;
    }
    arma::mat variable = mask;
    if (standardize) variable.diag().zeros();
    constant_gram = arma::find(variable == 0.0);
  }

  void transform(arguments_optim& x) override {
    B.zeros(p, p);
    B.elem(positions) = x.transparameters.elem(indices_in.elem(input_slots));
    B.diag() = arma::exp(B.diag());
    if (!B.is_finite() || arma::any(B.diag() <= 0.0)) {
      throw std::range_error("LLt: a free input or exponentiated pivot is outside the finite range.");
    }
    for (index i = 0; i < p; ++i) {
      for (index j = 0; j < i; ++j) {
        if (target(i,j) != 0.0) continue;
        long double sum = 0.0L;
        for (index k = 0; k < j; ++k) {
          sum += static_cast<long double>(B(i,k)) * B(j,k);
        }
        B(i,j) = static_cast<double>(-sum / B(j,j));
      }
    }
    if (!B.is_finite()) throw std::range_error("LLt: triangular completion overflowed.");

    L = B;
    if (standardize) {
      row_scale.set_size(p);
      scaled_norm.set_size(p);
      for (index i = 0; i < p; ++i) {
        row_scale[i] = arma::abs(B.row(i)).max();
        L.row(i) /= row_scale[i];
        scaled_norm[i] = arma::norm(L.row(i), 2);
        L.row(i) /= scaled_norm[i];
      }
      if (arma::any(L.diag() <= 0.0)) {
        throw std::range_error("LLt: a normalized pivot underflowed.");
      }
    }
    if (factor_output) {
      x.transparameters.elem(indices_out) = arma::vectorise(L);
    } else {
      arma::mat gram = L * L.t();
      gram.elem(constant_gram).zeros();
      if (standardize) gram.diag().ones();
      if (!gram.is_finite() || arma::any(gram.diag() <= 0.0)) {
        throw std::range_error("LLt: a covariance is outside the finite range.");
      }
      x.transparameters.elem(indices_out) = arma::vectorise(gram);
    }
  }

  void dtransform(arguments_optim& x) override {
    const arma::vec seed = x.dtransparameters.elem(indices_in.elem(input_slots));
    factor_differential(seed, dB, dL, row_speed);
    x.dtransparameters.elem(indices_out) = arma::vectorise(output_differential(dL));
  }

  void update_grad(arguments_optim& x) override {
    output_grad = output_seed(x.grad);
    add_input_gradient(x.grad, pullback(nullptr));
  }

  void update_dgrad(arguments_optim& x) override {
    const arma::mat dseed = output_seed(x.dgrad);
    add_input_gradient(x.dgrad, pullback(&dseed));
  }

  void jacobian(arguments_optim&) override {
    jacob.zeros(target.n_elem, indices_in.n_elem);
    arma::vec seed(positions.n_elem, arma::fill::zeros);
    arma::mat raw_dir, factor_dir;
    arma::vec speed;
    for (index k = 0; k < positions.n_elem; ++k) {
      seed[k] = 1.0;
      factor_differential(seed, raw_dir, factor_dir, speed);
      jacob.col(input_slots[k]) = arma::vectorise(output_differential(factor_dir));
      seed[k] = 0.0;
    }
  }

  void outcomes(arguments_optim&) override {
    matrices = {L};
    names_matrices = {"L"};
  }
};

inline LLt* choose_LLt(const Rcpp::List& setup) {
  const std::vector<arma::uvec> in = setup["indices_in"];
  const std::vector<arma::uvec> out = setup["indices_out"];
  if (in.size() != 1 || out.size() != 1) {
    throw std::invalid_argument("LLt: supply one input block and one output block.");
  }
  const arma::mat target = Rcpp::as<arma::mat>(setup["target"]);
  bool unit = true;
  if (setup.containsElementNamed("standardize")) {
    const double flag = Rcpp::as<double>(setup["standardize"]);
    if (flag != 0.0 && flag != 1.0) throw std::invalid_argument("LLt: standardize must be TRUE or FALSE.");
    unit = flag == 1.0;
  }
  const std::string output = setup.containsElementNamed("output")
    ? Rcpp::as<std::string>(setup["output"]) : "gram";
  if (output != "gram" && output != "factor") {
    throw std::invalid_argument("LLt: output must be 'gram' or 'factor'.");
  }
  return new LLt(target, unit, output == "factor", in[0], out[0]);
}

#endif  // LATENT_TRANSFORMATIONS_LLT_H
