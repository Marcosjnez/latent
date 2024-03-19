/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 28/05/2022
 *
 */

// #include <Rcpp/Benchmark/Timer.h>
// #include "structures.h"
// #include "manifolds.h"
// #include "multiple_rotations.h"
// #include "NPF.h"
// #include "criteria.h"
// #include "EFA_fit.h"
// #include "checks.h"

Rcpp::List rotate_efa(arguments_rotate x, rotation_manifold *manifold, rotation_criterion *criterion,
                      int random_starts, int cores) {

  arma::vec xf(random_starts);
  NTR x1;
  std::vector<NTR> x2(random_starts);
  // Select the optimization rutine:
  rotation_optim* algorithm = choose_optim(x.optim);

  // Perform multiple rotations with random starting values:

#ifdef _OPENMP
  omp_set_num_threads(cores);
#pragma omp parallel for
#endif
  for (int i=0; i < random_starts; ++i) {

    arguments_rotate args = x;

    args.T = random_orth(args.q, args.q);
    // args.T = arma::mat(args.q, args.q, arma::fill::eye);

    x2[i] = algorithm->optim(args, manifold, criterion);
    // x2[i] = ntr(args, manifold, criterion);
    // x2[i] = gd(args, manifold, criterion);

    xf[i] = std::get<3>(x2[i]);

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

  if(!convergence) {

    Rcpp::Rcout << "\n" << std::endl;
    Rcpp::warning("Failed rotation convergence");

  }

  Rcpp::List result;
  result["lambda"] = L;
  result["phi"] = Phi;
  result["T"] = T;
  result["f"] = f;
  result["iterations"] = iterations;
  result["convergence"] = convergence;

  result.attr("class") = "rotation";
  return result;

}

Rcpp::List efast(arma::mat X, int nfactors, std::string cor, std::string estimator,
                 Rcpp::CharacterVector char_rotation, std::string projection,
                 std::string missing,
                 Rcpp::Nullable<int> nullable_nobs,
                 Rcpp::Nullable<arma::mat> nullable_Target,
                 Rcpp::Nullable<arma::mat> nullable_Weight,
                 Rcpp::Nullable<arma::mat> nullable_PhiTarget,
                 Rcpp::Nullable<arma::mat> nullable_PhiWeight,
                 Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks,
                 Rcpp::Nullable<arma::vec> nullable_block_weights,
                 Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                 arma::vec gamma, arma::vec epsilon, arma::vec k, double w,
                 int random_starts, int cores,
                 Rcpp::Nullable<arma::vec> nullable_init,
                 Rcpp::Nullable<Rcpp::List> nullable_efa_control,
                 Rcpp::Nullable<Rcpp::List> nullable_rot_control) {

  Rcpp::Timer timer;
  Rcpp::List result;

  // cor structure:
  arguments_cor xcor;
  xcor.X = X;
  xcor.cor = cor;
  xcor.estimator = estimator;
  xcor.p = X.n_cols;
  xcor.q = nfactors;
  xcor.missing = missing;
  xcor.cores = cores;
  if(nullable_nobs.isNotNull()) {
    xcor.nobs = Rcpp::as<int>(nullable_nobs);
  }
  check_cor(xcor);
  Rcpp::List correlation_result = xcor.correlation_result;

  std::vector<std::string> rotation = Rcpp::as<std::vector<std::string>>(char_rotation);

  // efa structure:
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
  xefa.random_starts = random_starts;

  xefa.upper = arma::diagvec(xefa.R);
  xefa.nullable_efa_control = nullable_efa_control;
  xefa.nullable_init = nullable_init;

  check_efa(xefa);

  // Select one manifold:
  efa_manifold* efa_manifold = choose_efa_manifold(xefa.manifold);
  // Select the estimator:
  efa_criterion* efa_criterion = choose_efa_criterion(xefa.estimator);

  Rcpp::List efa_result = efa(xefa, efa_manifold, efa_criterion,
                              xefa.random_starts, xefa.cores);

  xefa.heywood = efa_result["heywood"];

  double df_null = xefa.p*(xefa.p-1)/2;
  double df = xefa.p*(xefa.p+1)/2 - (xefa.p*xefa.q + xefa.p - xefa.q*(xefa.q-1)/2);

  double f_null;
  if(estimator == "uls" || estimator == "pa") {
    f_null = 0.5*(arma::accu(xefa.R % xefa.R) - xefa.p);
  } else if(estimator == "dwls") {
    f_null = 0.5*arma::accu(xefa.R % xefa.R % xefa.W);
  } else if(estimator == "ml") {
    f_null = -arma::log_det_sympd(xefa.R);
  } else if(estimator == "minrank") {
    f_null = 0;
  }

  Rcpp::List modelInfo;
  modelInfo["correlation"] = xefa.R;
  modelInfo["smoothed"] = xefa.smoothed;
  modelInfo["estimator"] = xefa.estimator;
  modelInfo["projection"] = projection;
  modelInfo["rotation"] = rotation;
  modelInfo["nvars"] = xefa.p;
  modelInfo["nfactors"] = xefa.q;
  modelInfo["nobs"] = xefa.nobs;
  modelInfo["df"] = df;
  modelInfo["df_null"] = df_null;
  modelInfo["f_null"] = f_null;
  modelInfo["k"] = k;
  modelInfo["gamma"] = gamma;
  modelInfo["epsilon"] = epsilon;
  arma::vec clf_epsilon = {0.01};
  modelInfo["clf_epsilon"] = clf_epsilon;
  modelInfo["w"] = w;
  modelInfo["normalization"] = xefa.normalization;
  modelInfo["Target"] = nullable_Target;
  modelInfo["Weight"] = nullable_Weight;
  modelInfo["PhiTarget"] = nullable_PhiTarget;
  modelInfo["PhiWeight"] = nullable_PhiWeight;
  modelInfo["blocks"] = nullable_blocks;
  modelInfo["block_weights"] = nullable_block_weights;
  modelInfo["oblq_factors"] = nullable_oblq_factors;
  modelInfo["lower"] = xefa.lower;
  modelInfo["upper"] = xefa.upper;

  arma::mat lambda = efa_result["lambda"];

  // Structure of rotation arguments:

  arguments_rotate x;
  x.lambda = lambda;
  x.p = xefa.p, x.q = xefa.q;
  // x.lambda.set_size(x.p, x.q);
  x.Phi.set_size(x.q, x.q); x.Phi.eye();
  x.gamma = gamma, x.epsilon = epsilon, x.k = k, x.w = w, x.clf_epsilon = clf_epsilon;
  x.rotations = rotation;
  x.projection = projection;
  x.nullable_Target = nullable_Target;
  x.nullable_Weight = nullable_Weight;
  x.nullable_PhiTarget = nullable_PhiTarget;
  x.nullable_PhiWeight = nullable_PhiWeight;
  x.nullable_oblq_factors = nullable_oblq_factors;
  x.nullable_block_weights = nullable_block_weights;
  x.nullable_blocks = nullable_blocks;
  x.nullable_rot_control = nullable_rot_control;

  // Check rotation inputs and compute constants for rotation criteria:

  check_rotate(x, random_starts, cores);

  // Select one manifold:
  rotation_manifold* manifold = choose_manifold(x.projection);
  // Select one specific criteria or mixed criteria:
  rotation_criterion* criterion = choose_criterion(x.rotations, x.projection, x.cols_list);

  Rcpp::List rotation_result;

  arma::vec weigths;
  if (xefa.normalization == "kaiser") {

    weigths = sqrt(arma::sum(x.lambda % x.lambda, 1));
    x.lambda.each_col() /= weigths;

  } else if(xefa.normalization == "cureton") {

  } else if(xefa.normalization == "none") {

  } else {

    Rcpp::stop("Unkown normalization");

  }

  if(rotation.size() == 1 && rotation[0] == "none" || random_starts < 1) {

    arma::vec propVar = arma::diagvec(x.lambda.t() * x.lambda)/x.p;
    efa_result["lambda"] = x.lambda;
    efa_result["propVar"] = propVar;

    result["efa"] = efa_result;
    result["modelInfo"] = modelInfo;
    timer.step("elapsed");
    result["elapsed"] = timer;

    result.attr("class") = "efa";
    return result;

  } else {

    rotation_result = rotate_efa(x, manifold, criterion,
                                 random_starts, cores);

  }

  arma::mat L = rotation_result["lambda"];
  arma::mat Phi = rotation_result["phi"];

  if (xefa.normalization == "kaiser") {

    L.each_col() %= weigths;

  } else if(xefa.normalization == "cureton") {

  }

  // Force average positive lambda in all factors:

  // arma::vec v = arma::sign(arma::sum(L, 0));
  // L.each_row() /= v;

  for (int j=0; j < x.q; ++j) {
    if (arma::sum(L.col(j)) < 0) {
      L.col(j)   *= -1;
      Phi.col(j) *= -1;
      Phi.row(j) *= -1;
    }
  }

  // arma::vec uniquenesses = efa_result["uniquenesses"];
  // xefa.lambda = L;
  // xefa.phi = Phi;
  // xefa.psi = arma::diagmat(uniquenesses);
  // efa_criterion->H(xefa);
  // efa_criterion->H2(xefa);
  // xefa.lambda_indexes
  // Standard errors:
  // arma::mat inv_hessian;
  // if(!opt.hessian.is_sympd()) {
  //   arma::vec eigval;
  //   arma::mat eigvec;
  //   eig_sym(eigval, eigvec, opt.hessian);
  //   arma::vec d = arma::clamp(eigval, 0.01, eigval.max());
  //   inv_hessian = eigvec * arma::diagmat(1/d) * eigvec.t();
  // } else {
  //   inv_hessian = arma::inv_sympd(opt.hessian);
  // }
  // // arma::mat inv_hessian = arma::inv(opt.hessian);
  // arma::uvec indexes = trimatl_ind(arma::size(xcfa[0].R), 0);
  // arma::mat full_Sigma = asymptotic_normal(xcfa[0].R);
  // arma::mat Sigma = full_Sigma(indexes, indexes);
  // arma::mat dLPU_dS = opt.dLPU_dS[0].cols(indexes);
  // for(int i=1; i < opt.nblocks; ++i) {
  //   arma::uvec indexes = trimatl_ind(arma::size(xcfa[i].R), 0);
  //   arma::mat full_Sigma = asymptotic_normal(xcfa[i].R);
  //   Sigma = diagConc(Sigma, full_Sigma(indexes, indexes));
  //   dLPU_dS = arma::join_rows(dLPU_dS, opt.dLPU_dS[i].cols(indexes));
  // }
  // arma::mat B = dLPU_dS * Sigma * dLPU_dS.t();
  // arma::mat COV = inv_hessian * B * inv_hessian;
  // double denominator = (opt.total_nobs-1L)/opt.nblocks;
  // opt.se = sqrt(arma::diagvec(COV)/denominator);

  arma::vec propVar = arma::diagvec(Phi * L.t() * L)/x.p;

  rotation_result["lambda"] = L;
  rotation_result["phi"] = Phi;
  rotation_result["Rhat"] = efa_result["Rhat"];
  rotation_result["uniquenesses"] = efa_result["uniquenesses"];
  rotation_result["Rhat"] = efa_result["Rhat"];
  rotation_result["residuals"] = efa_result["residuals"];
  rotation_result["propVar"] = propVar;

  result["correlation"] = correlation_result;
  result["efa"] = efa_result;
  result["rotation"] = rotation_result;
  result["modelInfo"] = modelInfo;
  // result["oblq_indexes"] = x.oblq_indexes;
  // result["orth_indexes"] = x.orth_indexes;
  // result["loblq_indexes"] = x.loblq_indexes;

  timer.step("elapsed");
  result["elapsed"] = timer;

  result.attr("class") = "efa";
  return result;
}

// Do not export this (overloaded to support std::vector<std::string> rotation):
Rcpp::List efast(arma::mat X, int nfactors, std::string cor, std::string estimator,
                 std::vector<std::string> rotation, std::string projection,
                 std::string missing,
                 Rcpp::Nullable<int> nullable_nobs,
                 Rcpp::Nullable<arma::mat> nullable_Target,
                 Rcpp::Nullable<arma::mat> nullable_Weight,
                 Rcpp::Nullable<arma::mat> nullable_PhiTarget,
                 Rcpp::Nullable<arma::mat> nullable_PhiWeight,
                 Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> nullable_blocks,
                 Rcpp::Nullable<arma::vec> nullable_block_weights,
                 Rcpp::Nullable<arma::uvec> nullable_oblq_factors,
                 arma::vec gamma, arma::vec epsilon, arma::vec k, double w,
                 int random_starts, int cores,
                 Rcpp::Nullable<arma::vec> nullable_init,
                 Rcpp::Nullable<Rcpp::List> nullable_efa_control,
                 Rcpp::Nullable<Rcpp::List> nullable_rot_control) {

  Rcpp::Timer timer;
  Rcpp::List result;

  // cor structure:
  arguments_cor xcor;
  xcor.X = X;
  xcor.cor = cor;
  xcor.estimator = estimator;
  xcor.p = X.n_cols;
  xcor.q = nfactors;
  xcor.missing = missing;
  xcor.cores = cores;
  if(nullable_nobs.isNotNull()) {
    xcor.nobs = Rcpp::as<int>(nullable_nobs);
  }
  check_cor(xcor);
  Rcpp::List correlation_result = xcor.correlation_result;

  // efa structure:
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
  xefa.nullable_init = nullable_init;
  xefa.random_starts = random_starts;

  check_efa(xefa);

  // Select one manifold:
  efa_manifold* efa_manifold = choose_efa_manifold(xefa.manifold);
  // Select the estimator:
  efa_criterion* efa_criterion = choose_efa_criterion(xefa.estimator);

  Rcpp::List efa_result = efa(xefa, efa_manifold, efa_criterion,
                              xefa.random_starts, xefa.cores);

  xefa.heywood = efa_result["heywood"];

  arma::mat lambda = efa_result["lambda"];

  // Structure of rotation arguments:

  arguments_rotate x;
  x.lambda = lambda;
  x.p = xefa.p, x.q = xefa.q;
  // x.lambda.set_size(x.p, x.q);
  x.Phi.set_size(x.q, x.q); x.Phi.eye();
  x.gamma = gamma, x.epsilon = epsilon, x.k = k, x.w = w, x.clf_epsilon = {0.01};
  x.rotations = rotation;
  x.projection = projection;
  x.nullable_Target = nullable_Target;
  x.nullable_Weight = nullable_Weight;
  x.nullable_PhiTarget = nullable_PhiTarget;
  x.nullable_PhiWeight = nullable_PhiWeight;
  x.nullable_oblq_factors = nullable_oblq_factors;
  x.nullable_block_weights = nullable_block_weights;
  x.nullable_blocks = nullable_blocks;
  x.nullable_rot_control = nullable_rot_control;

  // Check inputs and compute constants for rotation criteria:

  check_rotate(x, random_starts, cores);

  // Model Info:

  double df_null = x.p*(x.p-1)/2;
  double df = x.p*(x.p+1)/2 - (x.p*x.q + x.p - x.q*(x.q-1)/2);

  double f_null;
  if(estimator == "uls" || estimator == "pa") {
    f_null = 0.5*(arma::accu(xefa.R % xefa.R) - x.p);
  } else if(estimator == "dwls") {
    f_null = 0.5*arma::accu(xefa.R % xefa.R % xefa.W);
  } else if(estimator == "ml") {
    f_null = -arma::log_det_sympd(xefa.R);
  } else if(estimator == "minrank") {
    f_null = 0;
  }

  Rcpp::List modelInfo;
  modelInfo["correlation"] = xefa.R;
  modelInfo["smoothed"] = xefa.smoothed;
  modelInfo["estimator"] = xefa.estimator;
  modelInfo["projection"] = x.projection;
  modelInfo["rotation"] = x.rotations;
  modelInfo["nvars"] = xefa.p;
  modelInfo["nfactors"] = xefa.q;
  modelInfo["nobs"] = xefa.nobs;
  modelInfo["df"] = df;
  modelInfo["df_null"] = df_null;
  modelInfo["f_null"] = f_null;
  modelInfo["k"] = x.k;
  modelInfo["gamma"] = x.gamma;
  modelInfo["epsilon"] = x.epsilon;
  modelInfo["clf_epsilon"] = x.clf_epsilon;
  modelInfo["w"] = x.w;
  modelInfo["normalization"] = x.normalization;
  modelInfo["Target"] = x.nullable_Target;
  modelInfo["Weight"] = x.nullable_Weight;
  modelInfo["PhiTarget"] = x.nullable_PhiTarget;
  modelInfo["PhiWeight"] = x.nullable_PhiWeight;
  modelInfo["blocks"] = x.nullable_blocks;
  modelInfo["block_weights"] = x.nullable_block_weights;
  modelInfo["oblq_factors"] = x.nullable_oblq_factors;
  modelInfo["lower"] = xefa.lower;
  modelInfo["upper"] = xefa.upper;


  // Select one manifold:
  rotation_manifold* manifold = choose_manifold(x.projection);
  // Select one specific criteria or mixed criteria:
  rotation_criterion* criterion = choose_criterion(x.rotations, x.projection, x.cols_list);

  Rcpp::List rotation_result;

  arma::vec weigths;
  if (xefa.normalization == "kaiser") {

    weigths = sqrt(arma::sum(x.lambda % x.lambda, 1));
    x.lambda.each_col() /= weigths;

  } else if(xefa.normalization == "cureton") {

  } else if(xefa.normalization == "none") {

  } else {

    Rcpp::stop("Unkown normalization");

  }

  if(rotation.size() == 1 && rotation[0] == "none" || random_starts < 1) {

    arma::vec propVar = arma::diagvec(x.lambda.t() * x.lambda)/x.p;
    efa_result["lambda"] = x.lambda;
    efa_result["propVar"] = propVar;

    result["efa"] = efa_result;
    result["modelInfo"] = modelInfo;
    timer.step("elapsed");
    result["elapsed"] = timer;

    result.attr("class") = "efa";
    return result;

  } else {

    rotation_result = rotate_efa(x, manifold, criterion,
                                 random_starts, cores);

  }

  arma::mat L = rotation_result["lambda"];
  arma::mat Phi = rotation_result["phi"];

  if (xefa.normalization == "kaiser") {

    L.each_col() %= weigths;

  } else if(xefa.normalization == "cureton") {

  }

  // Force average positive lambda in all factors:

  // arma::vec v = arma::sign(arma::sum(L, 0));
  // L.each_row() /= v;

  for (int j=0; j < x.q; ++j) {
    if (arma::sum(L.col(j)) < 0) {
      L.col(j)   *= -1;
      Phi.col(j) *= -1;
      Phi.row(j) *= -1;
    }
  }

  arma::vec propVar = arma::diagvec(Phi * L.t() * L)/x.p;

  rotation_result["lambda"] = L;
  rotation_result["phi"] = Phi;
  rotation_result["Rhat"] = efa_result["Rhat"];
  rotation_result["uniquenesses"] = efa_result["uniquenesses"];
  rotation_result["Rhat"] = efa_result["Rhat"];
  rotation_result["residuals"] = efa_result["residuals"];
  rotation_result["propVar"] = propVar;

  result["correlation"] = correlation_result;
  result["efa"] = efa_result;
  result["rotation"] = rotation_result;
  result["modelInfo"] = modelInfo;

  timer.step("elapsed");
  result["elapsed"] = timer;

  result.attr("class") = "efa";
  return result;
}
