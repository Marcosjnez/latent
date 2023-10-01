/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 04/09/202
 *
 */

// #include "structures.h"
// #include "auxiliary_manifolds.h"
// #include "se.h"

bool contains(std::vector<std::string> chars, std::string object) {

  /*
   * Check whether 'object' is in 'chars'
   */

  bool res = std::find(chars.begin(), chars.end(), object) != chars.end();

  return res;

}

arma::mat boot_sample(arma::mat X, bool replace) {

  /*
   * Generate a bootstrap sample from matrix 'X'
   */

  int p = X.n_rows;
  int q = X.n_cols;
  arma::uvec indexes(p);
  arma::uvec seq = consecutive(0, p-1);
  arma::mat X_boot(p, q);

  if(replace) {

    for(int i=0; i < q; ++i) {
      indexes = arma::randi<arma::uvec>(p, arma::distr_param(0, p-1));
      arma::vec v = X.col(i);
      X_boot.col(i) = v(indexes);
    }

  } else {

    for(int i=0; i < q; ++i) {
      // indexes = arma::randperm(p, p);
      std::random_shuffle(seq.begin(), seq.end());
      arma::vec v = X.col(i);
      X_boot.col(i) = v(seq);
      // X_boot.col(i) = arma::shuffle(X.col(i));
    }

  }

  return X_boot;

}

arma::vec eig_PAF(arma::mat S) {

  /*
   * Compute the eigenvalues of the reduced covariance matrix S
   */

  arma::mat inv_S = arma::inv(S); // Sometimes S is not positive-definite, so use inv()
  S.diag() -= 1/arma::diagvec(inv_S);
  arma::vec eigval_PAF = eig_sym(S);

  return eigval_PAF;

}

Rcpp::CharacterVector as_character(arma::vec x){

  int n = x.size();
  Rcpp::CharacterVector x_char(n);

  for(int i = 0; i < n; i++){
    x_char[i] = std::to_string(x[i]).substr(0, std::to_string(x[i]).find(".") + 2 + 1);
  }

  return x_char;

}

Rcpp::List out_pa(arma::umat dimensions, Rcpp::Nullable<arma::vec> nullable_quantile,
                  bool PCA, bool PAF, bool mean) {

  Rcpp::List result;

  bool quants_true = nullable_quantile.isNotNull();
  arma::vec quantile;
  if(quants_true) {
    quantile = Rcpp::as<arma::vec>(nullable_quantile);
  }

  if(PCA & PAF) {

    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(dimensions));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PCA", "PAF");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;

  } else if(PCA) {

    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(dimensions));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PCA");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;

  } else if(PAF) {

    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(dimensions));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PAF");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;

  }

  return result;

}

Rcpp::List pa(arma::mat X, arma::mat S, int n_boot, std::string type, Rcpp::Nullable<arma::vec> nullable_quantile,
              bool mean, bool replace, Rcpp::Nullable<std::vector<std::string>> nullable_PA, int cores) {

  /*
   * Perform parallel analysis with either PCA or PAF
   */

  Rcpp::List result;

  std::vector<std::string> PA;
  if(nullable_PA.isNotNull()) {
    PA = Rcpp::as<std::vector<std::string>>(nullable_PA);
  } else {
    PA = {"PCA", "PAF"};
  }

  int n = X.n_rows;
  int p = X.n_cols;

  bool quants_true = nullable_quantile.isNotNull();

  if(!quants_true & !mean) {
    Rcpp::stop("Please either enter a value for quantile or set mean = TRUE");
  }

  arma::vec quantile;

  if(quants_true) {
    quantile = Rcpp::as<arma::vec>(nullable_quantile);
  }

  int s = quantile.size();
  if(mean) s += 1;

  arma::cube X_boots(n, p, n_boot);
  arma::mat PCA_boot(p, n_boot);
  arma::mat PAF_boot(p, n_boot);

  arma::vec eigval_PCA;
  arma::vec eigval_PAF;
  arma::mat PCA_cutoff;
  arma::mat PAF_cutoff;
  arma::uvec PCA_groups(s);
  arma::uvec PAF_groups(s);

  bool PCA = contains(PA, "PCA");
  bool PAF = contains(PA, "PAF");

  // Generate the bootstrap samples

  if(type == "pearson") {

    for(int i=0; i < n_boot; ++i) {

      X_boots.slice(i) = boot_sample(X, replace);
      arma::mat S_boot;

      if(X_boots.slice(i).has_nan()) {
        S_boot = pairwise_cor(X_boots.slice(i));
      } else {
        S_boot = arma::cor(X_boots.slice(i));
      }

      if(PCA) PCA_boot.col(i) = eig_sym(S_boot);
      if(PAF) PAF_boot.col(i) = eig_PAF(S_boot);

    }

  } else if(type == "poly") {

// #ifdef _OPENMP
//       omp_set_num_threads(cores);
// #pragma omp parallel for
// #endif
    for(int i=0; i < n_boot; ++i) {

      X_boots.slice(i) = boot_sample(X, replace);
      polyfast_object polychor = poly_no_cores(X_boots.slice(i), "none", 0.00);
      // Rcpp::List polychor = polyfast(X_boots.slice(i), "none", 0L, false, 1L);
      arma::mat S_boot = std::get<0>(polychor);;
      if(PCA) PCA_boot.col(i) = eig_sym(S_boot);
      if(PAF) PAF_boot.col(i) = eig_PAF(S_boot);

    }

  } else {

    Rcpp::stop("Available correlation types: 'pearson' and 'poly'");

  }

  // Compute the reference eigenvalues

  if(PCA) {

    PCA_boot = PCA_boot.t();
    eigval_PCA = eig_sym(S);
    if(quants_true) PCA_cutoff = arma::quantile(PCA_boot, quantile);
    if(mean) PCA_cutoff.insert_rows(0,  arma::mean(PCA_boot));
    PCA_cutoff = PCA_cutoff.t();

    for(int i=0; i < s; ++i) {

      arma::umat booleans = arma::reverse(eigval_PCA > PCA_cutoff.col(i));
      arma::uvec ones = arma::find(booleans == 0);
      PCA_groups[i] = ones[0];

    }

  }
  // Rcpp::stop("Well until here");

  if(PAF) {

    PAF_boot = PAF_boot.t();
    eigval_PAF = eig_PAF(S);
    if(quants_true) PAF_cutoff = arma::quantile(PAF_boot, quantile);
    if(mean) PAF_cutoff.insert_rows(0, arma::mean(PAF_boot));
    PAF_cutoff = PAF_cutoff.t();

    for(int i=0; i < s; ++i) {

      arma::umat booleans = arma::reverse(eigval_PAF > PAF_cutoff.col(i));
      arma::uvec ones = arma::find(booleans == 0);
      PAF_groups[i] = ones[0];

    }

  }

  // Output

  arma::umat dimensions;

  if(PCA & PAF) {

    arma::umat dimensions = arma::join_rows(PCA_groups, PAF_groups);
    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(dimensions));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PCA", "PAF");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;
    result["PCA_cutoff"] = PCA_cutoff;
    result["PCA_boot"] = PCA_boot;
    result["PAF_cutoff"] = PAF_cutoff;
    result["PAF_boot"] = PAF_boot;

  } else if(PCA) {

    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(PCA_groups));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PCA");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;
    result["PCA_cutoff"] = PCA_cutoff;
    result["PCA_boot"] = PCA_boot;

  } else if(PAF) {

    Rcpp::NumericMatrix dims = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(PAF_groups));
    Rcpp::colnames(dims) = Rcpp::CharacterVector::create("PAF");
    Rcpp::CharacterVector names;

    if(quants_true) {
      names = as_character(quantile);
    }
    if(mean) {
      names.push_front("mean");
    }
    Rcpp::rownames(dims) = names;

    result["dimensions"] = dims;
    result["PAF_cutoff"] = PAF_cutoff;
    result["PAF_boot"] = PAF_boot;

  }

  return result;

}

Rcpp::List parallel(arma::mat X, int nboot, std::string cor, std::string missing,
                    Rcpp::Nullable<arma::vec> nullable_quantile,
                    bool mean, bool replace, Rcpp::Nullable<std::vector<std::string>> nullable_PA,
                    bool hierarchical, Rcpp::Nullable<Rcpp::List> nullable_efa,
                    int cores) {

  arguments_cor xcor;
  xcor.X = X;
  xcor.cor = cor;
  xcor.p = X.n_cols;
  xcor.missing = missing;
  xcor.cores = cores;

  Rcpp::List efa;
  if(nullable_efa.isNotNull()) {
    efa = nullable_efa;
  }
  if (efa.containsElementNamed("estimator")) {
    std::string estimator_ = efa["estimator"];
    xcor.estimator = estimator_;
  }

  check_cor(xcor);
  Rcpp::List correlation_result = xcor.correlation_result;

  Rcpp::List first_order = pa(xcor.X, xcor.R, nboot, cor, nullable_quantile, mean, replace, nullable_PA, cores);

  first_order["correlation"] = correlation_result;
  if(!hierarchical) return first_order;

  Rcpp::List result;
  result["first_order"] = first_order;
  arma::umat dims = first_order["dimensions"];
  int s = dims.n_rows;
  int q = dims.n_cols;

  arma::umat dims_2(s, q, arma::fill::zeros);

  int n = X.n_rows;
  int p = X.n_cols;

  // Arguments to pass to efa:

  arguments_efast x;

  pass_to_efast(efa, x);

  arma::uvec unique = arma::unique(dims(arma::find(dims > 1))); // unique elements greater than 1
  if(arma::max(unique) < 2) return first_order;
  int unique_size = unique.size();
  arma::uvec groups(unique_size);

  // Rcpp::stop("Well until here");

  for(int i=0; i < unique_size; ++i) {

      Rcpp::List fit = efast(xcor.R, unique[i], x.cor, x.estimator, x.rotation, x.projection,
                             xcor.missing, x.nullable_nobs,
                             x.nullable_Target, x.nullable_Weight,
                             x.nullable_PhiTarget, x.nullable_PhiWeight,
                             x.nullable_blocks,
                             x.nullable_block_weights,
                             x.nullable_oblq_factors,
                             x.gamma, x.epsilon, x.k, x.w,
                             x.random_starts, x.cores,
                             x.nullable_init, x.nullable_efa_control,
                             x.nullable_rot_control);

      Rcpp::List rot = fit["rotation"];
      arma::mat Phi = rot["phi"];
      arma::mat loadings = rot["lambda"];
      arma::mat L = loadings * Phi;
      arma::mat W = arma::solve(xcor.R, L);
      if(xcor.X.has_nan()) xcor.X.replace(arma::datum::nan, 0); // Avoid NAs in the multiplication
      arma::mat fs = xcor.X * W;
      arma::mat Sfs = arma::cor(fs);
      // Rcpp::stop("Well until here");
      Rcpp::List second_order = pa(fs, Sfs, nboot, "pearson", nullable_quantile, mean,
                                   replace, R_NilValue, false);

      arma::uvec indexes = arma::find(dims == unique[i]);
      arma::umat temp_dims = second_order["dimensions"];
      dims_2(indexes) = temp_dims(indexes);

    }

  std::vector<std::string> PA;
  if(nullable_PA.isNotNull()) {
    PA = Rcpp::as<std::vector<std::string>>(nullable_PA);
  } else {
    PA = {"PCA", "PAF"};
  }

  bool PCA = contains(PA, "PCA");
  bool PAF = contains(PA, "PAF");

  Rcpp::List second_order = out_pa(dims_2, nullable_quantile, PCA, PAF, mean);
  result["second_order"] = second_order;

  return result;

}

// Rcpp::List cv_eigen(arma::mat X, int N, bool hierarchical,
//                     Rcpp::Nullable<Rcpp::List> nullable_efa,
//                     int cores) {
//
//   arma::mat S = arma::cor(X);
//   int q = S.n_cols;
//   arma::mat CV_eigvals(N, q);
//   int p = X.n_rows;
//   arma::uvec indexes = consecutive(0, p-1);
//   int half = p/2;
//
// #ifdef _OPENMP
//   omp_set_num_threads(cores);
// #pragma omp parallel for
// #endif
//   for(int i=0; i < N; ++i) {
//
//     arma::uvec selected = arma::randperm(p, half);
//     arma::mat A = X.rows(selected);
//     arma::mat B = X;
//     B.shed_rows(selected);
//     arma::mat cor_A = arma::cor(A);
//     arma::mat cor_B = arma::cor(B);
//     arma::vec eigval;
//     arma::mat eigvec;
//     eig_sym(eigval, eigvec, cor_A);
//
//     arma::vec cv_values = eigvec.t() * cor_B * eigvec;
//     CV_eigvals.row(i) = arma::diagvec(cv_values);
//
//   }
//
//   arma::vec avg_CV_eigvals = arma::mean(CV_eigvals, 0);
//   arma::uvec which = arma::find(avg_CV_eigvals > 1);
//   int dim = which.size();
//
//   Rcpp::List result;
//   result["CV_eigvals"] = avg_CV_eigvals;
//   result["dim"] = dim;
//
//   if(dim <= 1 || !hierarchical) return result;
//
//   Rcpp::List efa;
//
//   if(nullable_efa.isNotNull()) {
//     efa = nullable_efa;
//   }
//
//   // Arguments to pass to efa:
//
//   arguments_efast x;
//
//   pass_to_efast(efa, x);
//
//   // efa:
//
//   Rcpp::List fit = efast(S, dim, x.cor, x.estimator, x.rotation, x.projection,
//                          "none", x.nullable_nobs,
//                          x.nullable_Target, x.nullable_Weight,
//                          x.nullable_PhiTarget, x.nullable_PhiWeight,
//                          x.nullable_blocks,
//                          x.nullable_block_weights,
//                          x.nullable_oblq_factors,
//                          x.gamma, x.epsilon, x.k, x.w,
//                          x.random_starts, x.cores,
//                          x.nullable_init, x.nullable_efa_control,
//                          x.nullable_rot_control);
//
//   Rcpp::List rot = fit["rotation"];
//   arma::mat Phi = rot["lambda"];
//   arma::mat loadings = rot["lambda"];
//   arma::mat L = loadings * Phi;
//   arma::mat W = arma::solve(S, L);
//   arma::mat fs = X * W;
//
//   arma::mat CV_eigvals2(N, dim);
//
//   omp_set_num_threads(cores);
// #pragma omp parallel for
//   for(int i=0; i < N; ++i) {
//
//     arma::uvec selected = arma::randperm(p, half);
//     arma::mat A = fs.rows(selected);
//     arma::mat B = fs;
//     B.shed_rows(selected);
//     arma::mat cor_A = arma::cor(A);
//     arma::mat cor_B = arma::cor(B);
//     arma::vec eigval;
//     arma::mat eigvec;
//     eig_sym(eigval, eigvec, cor_A);
//
//     arma::vec cv_values = eigvec.t() * cor_B * eigvec;
//     CV_eigvals2.row(i) = arma::diagvec(cv_values);
//
//   }
//
//   arma::vec avg_CV_eigvals2 = arma::mean(CV_eigvals2, 0);
//   // arma::uvec which2 = arma::find(avg_CV_eigvals2 > 1);
//   // int dim2 = which2.size();
//   arma::uvec ones = arma::find((avg_CV_eigvals2 < 1) == 1);
//   int dim2 = ones[0];
//
//   result["CV_eigvals2"] = avg_CV_eigvals2;
//   result["dim2"] = dim2;
//   result["fit"] = fit;
//   result["fs"] = fs;
//
//   return result;
//
// }
