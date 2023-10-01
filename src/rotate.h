/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include <Rcpp/Benchmark/Timer.h>
// #include "structures.h"
// #include "auxiliary_manifolds.h"
// #include "manifolds.h"
// #include "criteria.h"
// #include "NPF.h"
// #include "auxiliary_checks.h"
// #include "checks.h"

Rcpp::List rotate(arma::mat lambda, Rcpp::CharacterVector char_rotation,
                  std::string projection,
                  arma::vec gamma, arma::vec epsilon, arma::vec k, double w,
                  Rcpp::Nullable<arma::mat> nullable_Target,
                  Rcpp::Nullable<arma::mat> nullable_Weight,
                  Rcpp::Nullable<arma::mat> nullable_PhiTarget,
                  Rcpp::Nullable<arma::mat> nullable_PhiWeight,
                  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks,
                  Rcpp::Nullable<arma::vec> nullable_block_weights,
                  Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                  std::string normalization,
                  Rcpp::Nullable<Rcpp::List> nullable_rot_control,
                  int random_starts, int cores) {

  Rcpp::Timer timer;

  std::vector<std::string> rotation = Rcpp::as<std::vector<std::string>>(char_rotation);

  // Structure of rotation arguments:

  arguments_rotate x;
  x.p = lambda.n_rows, x.q = lambda.n_cols;
  x.lambda = lambda;
  x.Phi.set_size(x.q, x.q); x.Phi.eye();
  x.gamma = gamma, x.epsilon = epsilon, x.k = k, x.w = w;
  x.rotations = rotation;
  x.projection = projection;
  x.nullable_Target = nullable_Target;
  x.nullable_Weight = nullable_Weight;
  x.nullable_PhiTarget = nullable_PhiTarget;
  x.nullable_PhiWeight = nullable_PhiWeight;
  x.nullable_blocks = nullable_blocks;
  x.nullable_oblq_factors = nullable_oblq_factors;
  x.nullable_block_weights = nullable_block_weights;
  x.nullable_rot_control = nullable_rot_control;

  // Check inputs and compute constants for rotation criteria:
  check_rotate(x, random_starts, cores);

  // Select one manifold:
  rotation_manifold* manifold = choose_manifold(x.projection);
  // Select one specific criteria or mixed criteria:
  rotation_criterion* criterion = choose_criterion(x.rotations, x.projection, x.cols_list);
  // Select the optimization rutine:
  rotation_optim* algorithm = choose_optim(x.optim);

  arma::vec weigths;
  if (normalization == "kaiser") {

    weigths = sqrt(sum(x.lambda % x.lambda, 1));
    x.lambda.each_col() /= weigths;

  }

  arma::vec xf(random_starts);
  NTR x1;
  std::vector<NTR> x2(random_starts);

  // Perform multiple rotations with random starting values:

#ifdef _OPENMP
  omp_set_num_threads(cores);
#pragma omp parallel for
#endif
  for (int i=0; i < random_starts; ++i) {

    arguments_rotate args = x;
    args.T = random_orth(args.q, args.q);

    x2[i] = algorithm->optim(args, manifold, criterion);
    // x2[i] = ntr(args, manifold, criterion);
    // x2[i] = gd(args, manifold, criterion);

    xf[i] = std::get<3>(x2[i]);
    // xf[i] = x2[i].f;

  }

  // Choose the rotation with the smallest objective value:

  arma::uword index_minimum = index_min(xf);
  x1 = x2[index_minimum];

  arma::mat L = std::get<0>(x1);
  arma::mat Phi = std::get<1>(x1);
  if(Phi.is_empty()) {Phi.set_size(x.q, x.q); Phi.eye();}
  arma::mat T = std::get<2>(x1);
  double f = std::get<3>(x1);
  int iterations = std::get<4>(x1);
  bool convergence = std::get<5>(x1);

  // Force average positive lambda in all factors:

  // arma::vec v = arma::sign(arma::sum(L, 0));
  // L.each_row() /= v;

  for (int j=0; j < x.q; ++j) {
    if (sum(L.col(j)) < 0) {
      L.col(j)   *= -1;
      Phi.col(j) *= -1;
      Phi.row(j) *= -1;
    }
  }

  if(!convergence) {

    Rcpp::Rcout << "\n" << std::endl;
    Rcpp::warning("Failed rotation convergence");

  }

  if (normalization == "kaiser") {

    L.each_col() %= weigths;

  }

  arma::vec propVar = arma::diagvec(Phi * L.t() * L)/x.p;

  Rcpp::List modelInfo;
  modelInfo["lambda"] = lambda;
  modelInfo["rotation"] = rotation;
  modelInfo["projection"] = projection;
  modelInfo["n_vars"] = lambda.n_rows;
  modelInfo["nfactors"] = lambda.n_cols;
  // modelInfo["df"] = df;
  // modelInfo["df_null"] = df_null;
  // modelInfo["f_null"] = f_null;
  modelInfo["k"] = k;
  modelInfo["gamma"] = gamma;
  modelInfo["epsilon"] = epsilon;
  modelInfo["w"] = w;
  modelInfo["normalization"] = normalization;
  modelInfo["Target"] = nullable_Target;
  modelInfo["Weight"] = nullable_Weight;
  modelInfo["PhiTarget"] = nullable_PhiTarget;
  modelInfo["PhiWeight"] = nullable_PhiWeight;
  modelInfo["blocks"] = nullable_blocks;
  modelInfo["block_weights"] = nullable_block_weights;
  modelInfo["oblq_factors"] = nullable_oblq_factors;

  Rcpp::List result;
  result["lambda"] = L;
  result["phi"] = Phi;
  result["propVar"] = propVar;
  result["T"] = T;
  result["f"] = f;
  result["iterations"] = iterations;
  result["convergence"] = convergence;

  timer.step("elapsed");

  result["elapsed"] = timer;

  result.attr("class") = "rotation";
  return result;

}
