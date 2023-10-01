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

void extract_efa(arguments_efa& x, efa_manifold *manifold, efa_criterion *criterion) {

  efa_NTR x1;
  // Select the optimization routine:
  efa_optim* algorithm = choose_efa_optim(x.optim);

  arguments_efa args = x;
  // args.parameters = x.init;

  x1 = algorithm->optim(args, manifold, criterion);

  x.lambda = std::get<0>(x1);
  x.uniquenesses = std::get<1>(x1);
  x.Rhat = std::get<2>(x1);
  x.f = std::get<3>(x1);
  x.iterations = std::get<4>(x1);
  x.convergence = std::get<5>(x1);

  // if(!x.convergence) {
  //
  //   Rcpp::Rcout << "\n" << std::endl;
  //   Rcpp::warning("Failed rotation convergence");
  //
  // }

}

Rcpp::List efa(arma::vec psi, arma::mat R, int nfactors, std::string estimator,
               int efa_max_iter, double efa_factr, int lmm) {

  Rcpp::List result;

  arma::mat w, Rhat;
  arma::vec uniquenesses;

  int iteration = 0;

  if (estimator == "uls") {

    Rcpp::List optim_result = optim_rcpp(psi, R, nfactors, estimator, efa_max_iter, efa_factr, lmm);

    arma::vec psi_temp = optim_result["par"];
    psi = psi_temp;
    arma::mat reduced_R = R - arma::diagmat(psi);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
    arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
    for(int i=0; i < nfactors; ++i) {
      if(eigenvalues(i) < 0) eigenvalues(i) = 0;
    }
    arma::mat D = arma::diagmat(arma::sqrt(eigenvalues));

    w = A * D;
    arma::mat ww = w * w.t();

    uniquenesses = 1 - arma::diagvec(ww);

    Rhat = ww;
    Rhat.diag() = R.diag();

    bool convergence = false;
    int convergence_result = optim_result["convergence"];

    if(convergence_result == 0) convergence = true;

    result["f"] = optim_result["value"];
    result["convergence"] = convergence;

  } else if (estimator == "ml") {

    Rcpp::List optim_result = optim_rcpp(psi, R, nfactors, estimator, efa_max_iter, efa_factr, lmm);
    arma::vec psi_temp = optim_result["par"];
    psi = psi_temp;

    arma::vec sqrt_psi = arma::sqrt(psi);
    arma::mat sc = arma::diagmat(1/sqrt_psi);
    arma::mat Sstar = sc * R * sc;

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Sstar);

    arma::vec eigval2 = arma::reverse(eigval);
    arma::mat eigvec2 = arma::reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
    arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1)) - 1;
    for(int i=0; i < nfactors; ++i) {
      if(eigenvalues[i] < 0) eigenvalues[i] = 0;
    }
    arma::mat D = arma::diagmat(arma::sqrt(eigenvalues));

    w = A * D;
    w = arma::diagmat(sqrt_psi) * w;
    arma::mat ww = w * w.t();
    uniquenesses = 1 - arma::diagvec(ww);

    Rhat = ww;
    Rhat.diag() = R.diag();

    bool convergence = false;
    int convergence_result = optim_result["convergence"];

    if(convergence_result == 0) convergence = true;

    int p = R.n_cols;
    double value = optim_result["value"];
    double f = arma::log_det_sympd(Rhat) - arma::log_det_sympd(R) + arma::trace(R * arma::inv_sympd(Rhat)) - p;
    result["f"] = f;
    result["convergence"] = convergence;

  } else if (estimator == "pa") {

    Rcpp::List pa_result = principal_axis(psi, R, nfactors, 1e-03, efa_max_iter);

    arma::mat w_temp = pa_result["lambda"];
    arma::vec uniquenesses_temp = pa_result["uniquenesses"];
    arma::mat Rhat_temp = pa_result["Rhat"];

    result["f"] = pa_result["f"];
    result["convergence"] = pa_result["convergence"];
    w = w_temp;
    uniquenesses = uniquenesses_temp;
    Rhat = Rhat_temp;

    result["iterations"] = pa_result["iterations"];

  } else if (estimator == "minrank") {

    arma::vec communalities = sdp_cpp(R);

    psi = 1 - communalities;

    arma::mat reduced_R = R - arma::diagmat(psi);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, reduced_R);

    arma::vec eigval2 = arma::reverse(eigval);
    arma::mat eigvec2 = arma::reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
    arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
    for(int i=0; i < nfactors; ++i) {
      if(eigenvalues(i) < 0) eigenvalues(i) = 0;
    }
    arma::mat D = arma::diagmat(arma::sqrt(eigenvalues));

    w = A * D;
    arma::mat ww = w * w.t();
    uniquenesses = 1 - arma::diagvec(ww);

    Rhat = ww;
    Rhat.diag() = R.diag();

  } else {

    Rcpp::stop("Unkown estimator");

  }

  bool heywood = arma::any( uniquenesses < 0 );

  // Force average positive lambda in all factors:

  for (int j=0; j < nfactors; ++j) {
    if (sum(w.col(j)) < 0) {
      w.col(j)   *= -1;
    }
  }

  result["lambda"] = w;
  result["uniquenesses"] = uniquenesses;
  result["Rhat"] = Rhat;
  result["residuals"] = R - Rhat;
  result["Heywood"] = heywood;
  result["estimator"] = estimator;

  return result;
}

Rcpp::List efa(arguments_efa x, efa_manifold* manifold, efa_criterion* criterion,
               int random_starts, int cores) {

  Rcpp::List result;

  if(x.estimator == "ml" || x.estimator == "uls" || x.estimator == "dwls") {

    extract_efa(x, manifold, criterion);

  } else if(x.estimator == "pa") {

    Rcpp::List pa_result = principal_axis(x.parameters, x.R, x.q, x.eps, x.maxit);

    arma::mat w_temp = pa_result["lambda"];
    arma::vec uniquenesses_temp = pa_result["uniquenesses"];
    arma::mat Rhat_temp = pa_result["Rhat"];

    result["f"] = pa_result["f"];
    result["convergence"] = pa_result["convergence"];
    x.lambda = w_temp;
    x.uniquenesses = uniquenesses_temp;
    x.Rhat = Rhat_temp;
    x.iterations = pa_result["iterations"];

  } else if (x.estimator == "minrank") {

    arma::vec communalities = sdp_cpp(x.R);
    x.parameters = 1 - communalities;
    arma::mat reduced_R = x.R - diagmat(x.parameters);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, x.q-1));
    arma::vec eigenvalues = eigval2(arma::span(0, x.q-1));
    for(int i=0; i < x.q; ++i) {
      if(eigenvalues(i) < 0) eigenvalues(i) = 0;
    }
    arma::mat D = arma::diagmat(sqrt(eigenvalues));

    x.lambda = A * D;
    x.Rhat = x.lambda * x.lambda.t();
    x.uniquenesses = 1 - arma::diagvec(x.Rhat);
    x.Rhat.diag() = x.R.diag();
    x.iterations = 0;

  } else {

    Rcpp::stop("Unkown estimator");

  }

  x.heywood = arma::any( arma::vectorise(x.uniquenesses) <= 0 );

  result["lambda"] = x.lambda;
  result["uniquenesses"] = x.uniquenesses;
  result["Rhat"] = x.Rhat;
  result["residuals"] = x.R - x.Rhat;
  result["heywood"] = x.heywood;
  result["f"] = x.f;
  result["estimator"] = x.estimator;
  result["iterations"] = x.iterations;
  result["convergence"] = x.convergence;

  return result;

}
