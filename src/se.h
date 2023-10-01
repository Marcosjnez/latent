/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

// #include "structures.h"
// #include "checks.h"
// #include "manifolds.h"
// #include "criteria.h"
// #include "multiple_rotations.h"
// #include "method_derivatives.h"

Rcpp::List se(Rcpp::List fit,
              Rcpp::Nullable<int> nullable_nobs,
              Rcpp::Nullable<arma::mat> nullable_X,
              std::string type, double eta) {

  Rcpp::List modelInfo = fit["modelInfo"];
  Rcpp::Nullable<int> modelInfo_nobs = modelInfo["nobs"];
  int nobs;

  if(modelInfo_nobs.isNotNull() & nullable_nobs.isNull()) {
    nobs = Rcpp::as<int>(modelInfo_nobs);
  } else if(nullable_nobs.isNotNull()) {
    nobs = Rcpp::as<int>(nullable_nobs);
  } else {
    Rcpp::stop("Please, specify the sample size in the nobs argument");
  }

  arguments_rotate x;
  int random_starts = 1L, cores = 1L;

  // Overwrite x according to modelInfo:

  Rcpp::List rot = fit["rotation"];
  arma::mat L_ = rot["lambda"]; x.L = L_;
  arma::mat Phi_ = rot["phi"]; x.Phi = Phi_;
  arma::mat T_ = rot["T"]; x.T = T_;

  Rcpp::List efa = fit["efa"];
  arma::mat lambda_ = efa["lambda"]; x.lambda = lambda_;
  std::string projection_ = modelInfo["projection"]; x.projection = projection_;
  std::vector<std::string> rotation_ = modelInfo["rotation"]; x.rotations = rotation_;

  x.p = x.lambda.n_rows;
  x.q = x.lambda.n_cols;
  // x.lambda.set_size(x.p, x.q);
  // x.Phi.set_size(x.q, x.q); x.Phi.eye();
  arma::vec gamma_ = modelInfo["gamma"]; x.gamma = gamma_;
  arma::vec k_ = modelInfo["k"]; x.k = k_;
  arma::vec epsilon_ = modelInfo["epsilon"]; x.epsilon = epsilon_;
  arma::vec clf_epsilon_ = modelInfo["clf_epsilon"]; x.clf_epsilon = clf_epsilon_;
  double w_ = modelInfo["w"]; x.w = w_;

  arma::mat S_ = modelInfo["correlation"]; x.S = S_;
  std::string estimator = modelInfo["estimator"];
  Rcpp::Nullable<arma::mat> Target_ = modelInfo["Target"];
  x.nullable_Target = Target_;
  Rcpp::Nullable<arma::mat> Weight_ = modelInfo["Weight"];
  x.nullable_Weight = Weight_;
  Rcpp::Nullable<arma::mat> PhiTarget_ = modelInfo["PhiTarget"];
  x.nullable_PhiTarget = PhiTarget_;
  Rcpp::Nullable<arma::mat> PhiWeight_ = modelInfo["PhiWeight"];
  x.nullable_PhiWeight = PhiWeight_;
  Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks_ = modelInfo["blocks"];
  x.nullable_blocks = blocks_;
  Rcpp::Nullable<arma::vec> block_weights_ = modelInfo["block_weights"];
  x.nullable_block_weights = block_weights_;
  Rcpp::Nullable<arma::uvec> oblq_factors_ = modelInfo["oblq_factors"];
  x.nullable_oblq_factors = oblq_factors_;
  Rcpp::Nullable<arma::uvec> rot_control_ = modelInfo["oblq_factors"];
  x.nullable_rot_control = R_NilValue;

  check_rotate(x, random_starts, cores);

  rotation_manifold* manifold = choose_manifold(x.projection);
  rotation_criterion* criterion = choose_criterion(x.rotations, x.projection, x.cols_list);

  // Compute the rotation constraints:

  // manifold->param(x); // The sign of L and Phi may change (x.Inv_T is not necessary)
  criterion->F(x);
  criterion->gLP(x);
  criterion->hLP(x);
  manifold->g_constraints(x); // update x.d_constr

  // Rcpp::List xx;
  // xx["constr"] = x.d_constr;
  // xx["orth_indexes"] = x.orth_indexes;
  // xx["oblq_indexes"] = x.oblq_indexes;
  // xx["loblq_indexes"] = x.loblq_indexes;
  // return xx;

  int pq = x.p*x.q;
  int q_cor;
  int m;
  arma::uvec loblq_indexes;

  if(x.projection == "orth") {
    m = x.p*x.q + x.p;
  } else if(x.projection == "oblq") {
    q_cor = x.q*(x.q-1)/2;
    m = x.p*x.q + q_cor + x.p;
    loblq_indexes = arma::trimatl_ind(arma::size(x.Phi), -1);
  } else if(x.projection == "poblq") {
    loblq_indexes = x.loblq_indexes;
    q_cor = loblq_indexes.size();
    m = x.p*x.q + q_cor + x.p;
  } else {
    Rcpp::stop("Unkown projection");
  }

  // Compute the Hessian:

  arma::mat H, Hess;
  if(estimator == "uls") {
    Hess = hessian_uls(x.S, x.L, x.Phi, x.projection, loblq_indexes);
  } else if(estimator == "ml") {
    Hess = hessian_ml(x.S, x.L, x.Phi, x.projection, loblq_indexes);
  } else {
    Rcpp::stop("Standard errors are not implemented yet for this estimator");
  }

  // Add the constraints to the hessian matrix:
  int kk = x.d_constr.n_rows;
  H = arma::join_rows(Hess, x.d_constr.t());
  // Add zeros to make H a square matrix:
  arma::mat zeros(kk, kk);
  arma::mat C = arma::join_rows(x.d_constr, zeros);
  H = arma::join_cols(H, C);
  arma::mat H_inv = arma::inv(H);

  // arma::uvec caca = arma::trimatl_ind(arma::size(x.Phi), -1);
  // arma::uvec indexes_1(x.q);
  // for(int i=0; i < x.q; ++i) indexes_1[i] = ((i+1)*x.q) - (x.q-i);
  // Rcpp::List xx;
  // xx["Hess"] = Hess;
  // xx["H"] = H;
  // xx["H_inv"] = H_inv;
  // xx["constr"] = x.d_constr;
  // xx["orth_indexes"] = x.orth_indexes;
  // xx["oblq_indexes"] = x.oblq_indexes;
  // xx["loblq_indexes"] = x.loblq_indexes;
  // xx["projection"] = x.projection;
  // xx["estimator"] = estimator;
  // return xx;

  // Find A^{-1}BA^{-1}:
  arma::mat A_inv = H_inv(arma::span(0, m-1), arma::span(0, m-1));
  arma::mat BB = B(x.S, x.L, x.Phi, loblq_indexes, nullable_X, estimator, x.projection,
                   type, eta); // Asymptotic covariance of the correlation matrix
  arma::mat VAR = A_inv * BB * A_inv;
  arma::vec se = sqrt(arma::diagvec(VAR)/(nobs-1));

  // arma::uvec caca = arma::trimatl_ind(arma::size(x.Phi), -1);
  // arma::uvec indexes_1(x.q);
  // for(int i=0; i < x.q; ++i) indexes_1[i] = ((i+1)*x.q) - (x.q-i);
  // Rcpp::List xx;
  // xx["Hess"] = Hess;
  // xx["H"] = H;
  // xx["H_inv"] = H_inv;
  // xx["constr"] = x.d_constr;
  // xx["orth_indexes"] = x.orth_indexes;
  // xx["oblq_indexes"] = x.oblq_indexes;
  // xx["loblq_indexes"] = x.loblq_indexes;
  // xx["projection"] = x.projection;
  // xx["estimator"] = estimator;
  // xx["BB"] = BB;
  // return xx;

  arma::mat Lambda_se(x.p, x.q);
  arma::mat Phi_se(x.q, x.q);
  Phi_se.diag().zeros();
  arma::vec psi_se(x.p);

  for(int i=0; i < pq; ++i) Lambda_se[i] = se[i];
  if(x.projection == "orth") {
    int ij = 0;
    for(int i=pq; i < m; ++i) {
      psi_se[ij] = se[i];
      ++ij;
    }
  } else if(x.projection == "oblq") {
    int ij = 0;
    for(int i=pq; i < (pq+q_cor); ++i) {
      Phi_se[loblq_indexes(ij)] = se[i];
      ++ij;
    }
    Phi_se = arma::symmatl(Phi_se);
    ij = 0;
    for(int i=pq+q_cor; i < m; ++i) {
      psi_se[ij] = se[i];
      ++ij;
    }
  } else if(x.projection == "poblq") {
    int ij = 0;
    for(int i=pq; i < (pq+q_cor); ++i) {
      Phi_se[loblq_indexes(ij)] = se[i];
      ++ij;
    }
    Phi_se = arma::symmatl(Phi_se);
    ij = 0;
    for(int i=pq+q_cor; i < m; ++i) {
      psi_se[ij] = se[i];
      ++ij;
    }
  } else {
    Rcpp::stop("Unkown projection");
  }

  Rcpp::List result;
  result["Lambda_se"] = Lambda_se;
  result["Phi_se"] = Phi_se;
  result["psi_se"] = psi_se;

  return result;

}
