#ifndef LATENT_POLY_ACOV_PARALLEL_H
#define LATENT_POLY_ACOV_PARALLEL_H

#include <algorithm>
#include <exception>
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
#pragma omp for schedule(static)
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

} // namespace latent_asymptotic_poly
#endif
