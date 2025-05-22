/*
 * Author: Marcos Jiménez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

// #define ARMA_NO_DEBUG
// #define ARMA_DONT_USE_OPENMP
// #define ARMA_DONT_OPTIMISE_BAND
// #define ARMA_OPENMP_THREADS 10
// #define ARMA_DONT_OPTIMISE_SYMPD
// #define ARMA_DONT_USE_SUPERLU
// #define ARMA_DONT_USE_BLAS
// #define ARMA_OPENMP_THRESHOLD 240
// #define ARMA_64BIT_WORD
// #define ARMA_MAT_PREALLOC 4

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_num_threads() 1
  #define omp_set_num_threads() 1
#endif

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "structures.h"
#include "auxiliary.h"
#include "manifolds.h"
#include "transformations.h"
#include "estimators.h"
#include "optim.h"
#include "optimizer.h"
#include "polychorics.h"
#include "polyfast.h"
#include "asymptotic_cov.h"
#include "grad_comp.h"
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <R_ext/Print.h>   // For Rprintf()

// [[Rcpp::export]]
Rcpp::List optimizer(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer);

// [[Rcpp::export]]
Rcpp::List polyfast(arma::mat data, std::string missing = "pairwise.complete.cases",
                    const std::string acov = "none",
                    const std::string smooth = "none", double min_eigval = 0.001,
                    const int nboot = 1000L, const bool fit = false,
                    const int cores = 1L);

// [[Rcpp::export]]
Rcpp::List grad_comp(arma::vec parameters,
                     Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer,
                     double eps = 1e-04);
