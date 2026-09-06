#ifndef LATENT_POLY_ACOV_PARALLEL_H
#define LATENT_POLY_ACOV_PARALLEL_H

#include <algorithm>
#include <exception>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace latent_asymptotic_poly {

// Workers must only access C++ storage. R/Rcpp calls belong on the main thread.
// Preserve exceptions until every worker has finished; never throw across an
// OpenMP boundary. Each task owns a separate output region and exception slot.
template<typename Function>
inline int acov_parallel_for(const arma::uword size, const int cores,
                              const Function& function) {
#ifdef _OPENMP
  if(cores > 1 && size > 1L && !omp_in_parallel()) {
    const int threads = static_cast<int>(std::min<arma::uword>(size, cores));
    std::vector<std::exception_ptr> errors(size);
    int used = 1;
#pragma omp parallel num_threads(threads)
    {
#pragma omp single
      used = omp_get_num_threads();
#pragma omp for schedule(dynamic, 1)
      for(arma::uword i = 0L; i < size; ++i) {
        try {
          function(i);
        } catch(...) {
          errors[i] = std::current_exception();
        }
      }
    }
    for(const auto& error : errors) {
      if(error) std::rethrow_exception(error);
    }
    return used;
  }
#endif
  (void)cores;
  for(arma::uword i = 0L; i < size; ++i) function(i);
  return 1;
}

// A symmetric crossproduct without nested BLAS calls in OpenMP workers. Keep
// the existing Armadillo/BLAS path for one core and small matrices. Tile pairs
// own disjoint cells of the output, so no floating-point reduction is shared.
inline arma::mat acov_crossprod(const arma::mat& scores, const int cores,
                                int* cores_used = nullptr) {
#ifdef _OPENMP
  const arma::uword p = scores.n_cols;
  const arma::uword n = scores.n_rows;
  if(cores > 1 && p >= 64L && !omp_in_parallel()) {
    const arma::uword width = 16L;
    const arma::uword blocks = (p+width-1L)/width;
    std::vector<std::pair<arma::uword, arma::uword>> tiles;
    tiles.reserve(blocks*(blocks+1L)/2L);
    for(arma::uword j = 0L; j < blocks; ++j) {
      for(arma::uword k = j; k < blocks; ++k) {
        tiles.emplace_back(j*width, k*width);
      }
    }
    arma::mat result(p, p, arma::fill::none);
    const int used = acov_parallel_for(tiles.size(), cores,
      [&](const arma::uword tile) {
      const arma::uword first = tiles[tile].first;
      const arma::uword second = tiles[tile].second;
      const arma::uword last_first = std::min(p, first+width);
      const arma::uword last_second = std::min(p, second+width);
      for(arma::uword j = first; j < last_first; ++j) {
        const double* x = scores.colptr(j);
        for(arma::uword k = std::max(second, j); k < last_second; ++k) {
          const double* y = scores.colptr(k);
          double value = 0.0;
#pragma omp simd reduction(+:value)
          for(arma::uword r = 0L; r < n; ++r) value += x[r]*y[r];
          result(j, k) = value;
          result(k, j) = value;
        }
      }
    });
    if(cores_used) *cores_used = std::max(*cores_used, used);
    return result;
  }
#endif
  (void)cores;
  if(cores_used) *cores_used = std::max(*cores_used, 1);
  return scores.t()*scores;
}

} // namespace latent_asymptotic_poly
#endif
