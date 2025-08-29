/*
 * Author: Marcos Jiménez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 17/07/2025
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

#include <unordered_map>
#include <functional>
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "structures.h"
#include "auxiliary.h"
#include "polychorics.h"
#include "polyfast.h"
#include "asymptotic_cov.h"
#include "manifolds.h"
#include "transformations.h"
#include "estimators.h"
#include "optim.h"
#include "optimizer.h"
#include "grad_comp.h"
#include "vcov.h"
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
arma::mat asymptotic_normal(arma::mat P);

// [[Rcpp::export]]
arma::mat asymptotic_general(arma::mat X);

// [[Rcpp::export]]
arma::mat DACOV2(int n, arma::mat poly,
                 std::vector<std::vector<std::vector<int>>> tabs,
                 std::vector<std::vector<double>> taus,
                 std::vector<std::vector<double>> mvphis);

// [[Rcpp::export]]
arma::mat orth(arma::mat X);

// [[Rcpp::export]]
arma::mat oblq(arma::mat X);

// [[Rcpp::export]]
arma::mat poblq(arma::mat X, arma::mat target);

// [[Rcpp::export]]
arma::mat rorth(int p, int q);

// [[Rcpp::export]]
arma::mat roblq(int p, int q);

// [[Rcpp::export]]
arma::mat rpoblq(int p, int q, arma::mat target);

// [[Rcpp::export]]
Rcpp::List grad_comp(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer,
                     std::string compute = "all",
                     double eps = 1e-04);

// [[Rcpp::export]]
Rcpp::List vcov_all(Rcpp::List control_manifold,
                    Rcpp::List control_transform,
                    Rcpp::List control_estimator,
                    Rcpp::List control_optimizer,
                    arma::mat H);

