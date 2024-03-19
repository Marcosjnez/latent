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
#include "efa_manifolds.h"
#include "efa_criteria.h"
#include "auxiliary_efa_optim.h"
#include "efa_optim.h"
#include "rotation_manifolds.h"
#include "rotation_criteria.h"
#include "auxiliary_rotation_optim.h"
#include "rotation_checks.h"
#include "rotation_optim.h"
#include "rotate.h"
#include "polychorics.h"
#include "cor_manifolds.h"
#include "cor_criteria.h"
#include "auxiliary_cor_optim.h"
#include "cor_optim.h"
#include "asymptotic_cov.h"
#include "efa_checks.h"
#include "polyfast.h"
#include "checks_cor.h"
#include "efa_fit.h"
#include "efa.h"
#include "efast.h"
#include "sl.h"
#include "bifad.h"
#include "GSLiD.h"
#include "botmin.h"
#include "bifactor.h"
#include "method_derivatives.h"
#include "se.h"
#include "dimensionality.h"
#include "check_deriv.h"
#include "cfa_checks.h"
#include "cfa_criteria.h"
#include "cfa_manifolds.h"
#include "auxiliary_cfa_optim.h"
#include "cfa_optim.h"
#include "cfa.h"

// [[Rcpp::export]]
arma::mat random_orth(int p, int q);

// [[Rcpp::export]]
arma::mat random_oblq(int p, int q);

// [[Rcpp::export]]
arma::mat random_poblq(int p, int q, arma::uvec oblq_factors);

// [[Rcpp::export]]
arma::mat retr_orth(arma::mat X);

// [[Rcpp::export]]
arma::mat retr_oblq(arma::mat X);

// [[Rcpp::export]]
arma::mat retr_poblq(arma::mat X,
                     Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
                     Rcpp::Nullable<arma::mat> PhiTarget = R_NilValue);

// [[Rcpp::export]]
Rcpp::List sl(arma::mat X, int n_generals, int n_groups,
              std::string cor = "pearson",
              std::string estimator = "uls",
              std::string missing = "pairwise.complete.cases",
              Rcpp::Nullable<int> nobs = R_NilValue,
              Rcpp::Nullable<Rcpp::List> first_efa = R_NilValue,
              Rcpp::Nullable<Rcpp::List> second_efa = R_NilValue,
              int cores = 1L);

// [[Rcpp::export]]
Rcpp::List rotate(arma::mat lambda,
                  Rcpp::CharacterVector rotation = Rcpp::CharacterVector::create("oblimin"),
                  std::string projection = "oblq",
                  arma::vec gamma = 0, arma::vec epsilon = Rcpp::NumericVector::create(0.01),
                  arma::vec k = 0, double w = 1,
                  Rcpp::Nullable<arma::mat> Target = R_NilValue,
                  Rcpp::Nullable<arma::mat> Weight = R_NilValue,
                  Rcpp::Nullable<arma::mat> PhiTarget = R_NilValue,
                  Rcpp::Nullable<arma::mat> PhiWeight = R_NilValue,
                  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks = R_NilValue,
                  Rcpp::Nullable<arma::vec> block_weights = R_NilValue,
                  Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
                  std::string normalization = "none",
                  Rcpp::Nullable<Rcpp::List> rot_control = R_NilValue,
                  int random_starts = 1, int cores = 1);

// [[Rcpp::export]]
Rcpp::List efast(arma::mat X, int nfactors, std::string cor = "pearson",
                 std::string estimator = "uls",
                 Rcpp::CharacterVector rotation = Rcpp::CharacterVector::create("oblimin"),
                 std::string projection = "oblq",
                 std::string missing = "pairwise.complete.cases",
                 Rcpp::Nullable<int> nobs = R_NilValue,
                 Rcpp::Nullable<arma::mat> Target = R_NilValue,
                 Rcpp::Nullable<arma::mat> Weight = R_NilValue,
                 Rcpp::Nullable<arma::mat> PhiTarget = R_NilValue,
                 Rcpp::Nullable<arma::mat> PhiWeight = R_NilValue,
                 Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks = R_NilValue,
                 Rcpp::Nullable<arma::vec> block_weights = R_NilValue,
                 Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
                 arma::vec gamma = 0, arma::vec epsilon = Rcpp::NumericVector::create(0.01),
                 arma::vec k = 0, double w = 1,
                 int random_starts = 1, int cores = 1,
                 Rcpp::Nullable<arma::vec> init = R_NilValue,
                 Rcpp::Nullable<Rcpp::List> efa_control = R_NilValue,
                 Rcpp::Nullable<Rcpp::List> rot_control = R_NilValue);

// [[Rcpp::export]]
arma::mat get_target(arma::mat loadings, Rcpp::Nullable<arma::mat> Phi, double cutoff = 0);

// Rcpp::List bifad(arma::mat R, int n_generals, int n_groups,
//                  std::string projection = "orth",
//                  Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
//                  double cutoff = 0,
//                  std::string normalization = "none",
//                  Rcpp::Nullable<int> nobs = R_NilValue,
//                  Rcpp::Nullable<Rcpp::List> first_efa = R_NilValue,
//                  Rcpp::Nullable<Rcpp::List> second_efa = R_NilValue,
//                  Rcpp::Nullable<Rcpp::List> rot_control = R_NilValue,
//                  int random_starts = 1, int cores = 1);

// [[Rcpp::export]]
Rcpp::List bifactor(arma::mat X, int n_generals, int n_groups,
                   std::string method = "GSLiD",
                   std::string cor = "pearson",
                   std::string estimator = "uls",
                   std::string projection = "oblq",
                   std::string missing = "pairwise.complete.cases",
                   Rcpp::Nullable<int> nobs = R_NilValue,
                   Rcpp::Nullable<arma::mat> PhiTarget = R_NilValue,
                   Rcpp::Nullable<arma::mat> PhiWeight = R_NilValue,
                   Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks = R_NilValue,
                   Rcpp::Nullable<arma::vec> block_weights = R_NilValue,
                   Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
                   Rcpp::Nullable<arma::mat> init_Target = R_NilValue,
                   int maxit = 20, double cutoff = 0, std::string normalization = "none",
                   double w = 1, int random_starts = 1, int cores = 1,
                   Rcpp::Nullable<arma::vec> init = R_NilValue,
                   Rcpp::Nullable<Rcpp::List> efa_control = R_NilValue,
                   Rcpp::Nullable<Rcpp::List> rot_control = R_NilValue,
                   Rcpp::Nullable<Rcpp::List> first_efa = R_NilValue,
                   Rcpp::Nullable<Rcpp::List> second_efa = R_NilValue,
                   bool verbose = true);

// [[Rcpp::export]]
arma::mat asymp_cov(arma::mat R,
                    Rcpp::Nullable<arma::mat> X = R_NilValue,
                    double eta = 1, std::string type = "normal");

// [[Rcpp::export]]
Rcpp::List se(Rcpp::List fit = R_NilValue,
              Rcpp::Nullable<int> nobs = R_NilValue,
              Rcpp::Nullable<arma::mat> X = R_NilValue,
              std::string type = "normal", double eta = 1);

// [[Rcpp::export]]
Rcpp::List parallel(arma::mat X, int nboot = 100, std::string cor = "pearson",
                    std::string missing = "pairwise.complete.cases",
                    Rcpp::Nullable<arma::vec> quant = R_NilValue,
                    bool mean = false, bool replace = false,
                    Rcpp::Nullable<std::vector<std::string>> PA = R_NilValue,
                    bool hierarchical = false, Rcpp::Nullable<Rcpp::List> efa = R_NilValue,
                    int cores = 1);

// Rcpp::List cv_eigen(arma::mat X, int N = 100, bool hierarchical = false,
                    // Rcpp::Nullable<Rcpp::List> efa = R_NilValue, int cores = 1);

// Rcpp::List check_deriv(arma::mat L, arma::mat Phi,
//                        arma::mat dL, arma::mat dP,
//                        Rcpp::CharacterVector rotation = Rcpp::CharacterVector::create("oblimin"),
//                        std::string projection = "oblq",
//                        Rcpp::Nullable<arma::mat> Target = R_NilValue,
//                        Rcpp::Nullable<arma::mat> Weight = R_NilValue,
//                        Rcpp::Nullable<arma::mat> PhiTarget = R_NilValue,
//                        Rcpp::Nullable<arma::mat> PhiWeight = R_NilValue,
//                        Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks = R_NilValue,
//                        Rcpp::Nullable<arma::vec> block_weights = R_NilValue,
//                        Rcpp::Nullable<arma::uvec> oblq_factors = R_NilValue,
//                        arma::vec gamma = 0, arma::vec epsilon = Rcpp::NumericVector::create(0.01),
//                        arma::vec k = 0, double w = 1);

// [[Rcpp::export]]
Rcpp::List polyfast(arma::mat data, std::string missing = "pairwise.complete.cases",
                    const std::string acov = "none",
                    const std::string smooth = "none", double min_eigval = 0.001,
                    const int nboot = 1000L, const bool fit = false,
                    const int cores = 1L);

// [[Rcpp::export]]
Rcpp::List cfa(arma::vec parameters,
               std::vector<arma::mat> X,
               arma::ivec nfactors,
               arma::ivec nobs,
               std::vector<arma::mat> lambda,
               std::vector<arma::mat> phi,
               std::vector<arma::mat> psi,
               std::vector<arma::uvec> lambda_indexes,
               std::vector<arma::uvec> phi_indexes,
               std::vector<arma::uvec> psi_indexes,
               std::vector<arma::uvec> target_indexes,
               std::vector<arma::uvec> targetphi_indexes,
               std::vector<arma::uvec> targetpsi_indexes,
               std::vector<arma::uvec> free_indices_phi,
               std::vector<arma::uvec> free_indices_psi,
               Rcpp::CharacterVector cor = Rcpp::CharacterVector::create("pearson"),
               Rcpp::CharacterVector estimator = Rcpp::CharacterVector::create("uls"),
               Rcpp::CharacterVector projection = Rcpp::CharacterVector::create("id"),
               Rcpp::CharacterVector missing = Rcpp::CharacterVector::create("pairwise.complete.cases"),
               std::string se = "robust",
               int random_starts = 1L, int cores = 1L,
               Rcpp::Nullable<Rcpp::List> control = R_NilValue);
