/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include <Rcpp/Benchmark/Timer.h>
// #include "structures.h"
// #include "manifolds.h"
// #include "criteria.h"
// #include "checks.h"
// #include "multiple_rotations.h"
// #include "EFA.h"

Rcpp::List bifactor(arma::mat X, int n_generals, int n_groups,
                    std::string method, std::string cor,
                    std::string estimator, std::string projection,
                    std::string missing,
                    Rcpp::Nullable<int> nullable_nobs,
                    Rcpp::Nullable<arma::mat> nullable_PhiTarget,
                    Rcpp::Nullable<arma::mat> nullable_PhiWeight,
                    Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks,
                    Rcpp::Nullable<arma::vec> nullable_block_weights,
                    Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                    Rcpp::Nullable<arma::mat> nullable_Target,
                    int maxit, double cutoff, std::string normalization,
                    double w, int random_starts, int cores,
                    Rcpp::Nullable<arma::vec> nullable_init,
                    Rcpp::Nullable<Rcpp::List> nullable_efa_control,
                    Rcpp::Nullable<Rcpp::List> nullable_rot_control,
                    Rcpp::Nullable<Rcpp::List> nullable_first_efa,
                    Rcpp::Nullable<Rcpp::List> nullable_second_efa,
                    bool verbose) {

  Rcpp::Timer timer;
  Rcpp::List result, SL_result;

  arguments_cor xcor;
  xcor.X = X;
  xcor.cor = cor;
  xcor.estimator = estimator;
  xcor.p = X.n_cols;
  xcor.q = n_generals + n_groups;
  xcor.missing = missing;
  xcor.cores = cores;
  if(nullable_nobs.isNotNull()) {
    xcor.nobs = Rcpp::as<int>(nullable_nobs);
  }
  check_cor(xcor);
  Rcpp::List correlation_result = xcor.correlation_result;

  result["correlation"] = correlation_result;

  arguments_efa xefa;
  xefa.X = xcor.X;
  xefa.R = xcor.R;
  xefa.W = xcor.W;
  xefa.cor = xcor.cor;
  xcor.estimator = xcor.estimator;
  xefa.p = xcor.p;
  xefa.q = xcor.q;
  xefa.missing = xcor.missing;
  xefa.cores = xcor.cores;
  xefa.nobs = xcor.nobs;

  xefa.upper = arma::diagvec(xefa.R);
  xefa.nullable_efa_control = nullable_efa_control;
  xefa.nullable_first_efa = nullable_first_efa;
  xefa.nullable_second_efa = nullable_second_efa;
  xefa.nullable_init = nullable_init;

  check_efa(xefa);

  if(maxit < 1) Rcpp::stop("maxit must be an integer greater than 0");
  if(cutoff < 0) Rcpp::stop("cutoff must be nonnegative");

  int n = xefa.R.n_rows;
  int nfactors = n_generals + n_groups;

  if(method == "botmin") {

    Rcpp::List botmin_result = botmin(xefa.R, n_generals, n_groups,
                                      xefa.estimator, projection,
                                      nullable_nobs,
                                      nullable_oblq_factors,
                                      cutoff, random_starts, cores,
                                      xefa.nullable_efa_control,
                                      nullable_rot_control);

    result = botmin_result;

  } else if(method == "bifad") {

    Rcpp::List bifad_result = bifad(xefa.R, n_generals, n_groups,
                                    projection,
                                    nullable_oblq_factors,
                                    cutoff, normalization,
                                    nullable_nobs,
                                    xefa.nullable_first_efa,
                                    nullable_second_efa,
                                    nullable_rot_control,
                                    random_starts, cores);

    result = bifad_result;

  } else if(method == "SL") {

    xefa.q = n_groups;
    SL_result = sl(xefa.R, n_generals, n_groups, cor, estimator, "none", nullable_nobs,
                   xefa.nullable_first_efa, xefa.nullable_second_efa, cores);
    result["SL"] = SL_result;

  } else if(method == "GSLiD") {

    // Create initial target with Schmid-Leiman (SL) if there is no custom initial target:

    if(nullable_Target.isNull()) {

      xefa.q = n_groups;
      SL_result = sl(xefa.R, n_generals, n_groups, cor, estimator, "none", nullable_nobs,
                     xefa.nullable_first_efa, xefa.nullable_second_efa, cores);

      // Create the factor correlation matrix for the SL solution:

      arma::mat new_Phi(nfactors, nfactors, arma::fill::eye);

      if(n_generals > 1) {

        Rcpp::List second_order_solution = SL_result["second_order_solution"];
        Rcpp::List second_order_solution_rotation = second_order_solution["rotation"];
        arma::mat Phi_generals = second_order_solution_rotation["phi"];
        new_Phi(arma::span(0, n_generals-1), arma::span(0, n_generals-1)) = Phi_generals;

      }

      // SL loadings:

      arma::mat SL_loadings = SL_result["lambda"];
      arma::mat loadings = SL_loadings;

      // Create initial target:

      arma::mat Target;
      update_target(n_generals, n, nfactors, loadings, new_Phi, cutoff, Target);
      SEXP Target_ = Rcpp::wrap(Target);
      nullable_Target = Target_;

    }

    Rcpp::List GSLiD_result = GSLiD(xefa.R, n_generals, n_groups,
                                    estimator, projection,
                                    nullable_nobs,
                                    nullable_Target,
                                    nullable_PhiTarget, nullable_PhiWeight,
                                    w, maxit, random_starts, cores,
                                    nullable_init,
                                    xefa.nullable_efa_control,
                                    nullable_rot_control,
                                    nullable_blocks,
                                    nullable_block_weights,
                                    nullable_oblq_factors,
                                    cutoff, verbose);

    result = GSLiD_result; // bifactor and modelInfo outputs
    result["SL"] = SL_result;

  } else {

    Rcpp::stop("Unkown bifactor method");

  }

  timer.step("elapsed");

  result["elapsed"] = timer;

  result.attr("class") = "bifactor";
  return result;

}
