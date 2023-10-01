/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 23/09/2023
 *
 */

// [[Rcpp::export]]
arma::mat diagConc(arma::mat X1, arma::mat X2) {

  int p1 = X1.n_rows;
  int q1 = X1.n_cols;
  int p2 = X2.n_rows;
  int q2 = X2.n_cols;
  arma::mat augmented(p1+p2, q1+q2, arma::fill::zeros);
  augmented(arma::span(0, p1-1L), arma::span(0, q1-1L)) = X1;
  augmented(arma::span(p1, p1+p2-1L), arma::span(q1, q1+q2-1L)) = X2;

  return augmented;

}

Rcpp::List resizeList( const Rcpp::List& x, int n ){
  int oldsize = x.size() ;
  Rcpp::List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

Rcpp::CharacterVector resizeChar( const Rcpp::CharacterVector& x, int n ){
  int oldsize = x.size() ;
  Rcpp::CharacterVector y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

arma::vec diagcov(arma::mat X) {

  // Get the variance of variabes with missing data
  const size_t numCols = X.n_cols;

  // Initialize the correlation matrix
  arma::vec diag_cov(numCols);

  // Loop over all pairs of columns
  for (size_t i = 0; i < numCols; ++i) {
    // Get the columns for the pair (i, j)
    arma::vec col = X.col(i);

    // Find indices where both columns have non-NaN values
    arma::uvec validIndices = arma::find_finite(col);

    // Extract non-NaN values from both columns
    arma::vec validCol = col(validIndices);
    validCol -= arma::mean(validCol);

    // Calculate the correlation between the two columns
    double d = arma::accu(validCol % validCol) / (validIndices.n_elem-1);

    // Assign the correlation value to the correlation matrix
    diag_cov(i) = d;
  }

  return diag_cov;

}

arma::mat pairwise_cor(arma::mat X) {

  const size_t numCols = X.n_cols;

  // Initialize the correlation matrix
  arma::mat corrMatrix(numCols, numCols, arma::fill::eye);

  // Loop over all pairs of columns
  for (size_t i = 0; i < (numCols-1); ++i) {
    for (size_t j = (i+1); j < numCols; ++j) {
      // Get the columns for the pair (i, j)
      arma::vec col1 = X.col(i);
      arma::vec col2 = X.col(j);

      // Find indices where both columns have non-NaN values
      arma::uvec validIndices = arma::find_finite(col1 % col2);

      // Extract non-NaN values from both columns
      arma::vec validCol1 = col1(validIndices);
      arma::vec validCol2 = col2(validIndices);

      // Calculate the correlation between the two columns
      double correlation = as_scalar(arma::cor(validCol1, validCol2));

      // Assign the correlation value to the correlation matrix
      corrMatrix(i, j) = corrMatrix(j, i) = correlation;
    }
  }

  return corrMatrix;
}

void missingness(arguments_cor& x) {

  bool missing = x.X.has_nan();

  if(missing) {

    if(x.missing == "impute.mean") {
      for(int i=0; i < x.p; ++i) {
        arma::vec xvec = x.X.col(i);
        if(x.cor == "poly") {
          int number = std::round(arma::mean(xvec(arma::find_finite(xvec))));
          x.X.col(i).replace(arma::datum::nan, number);
        } else {
          double number = arma::mean(xvec(arma::find_finite(xvec)));
          x.X.col(i).replace(arma::datum::nan, number);
        }
      }
    } else if(x.missing == "impute.median") {
      for(int i=0; i < x.p; ++i) {
        arma::vec xvec = x.X.col(i);
        if(x.cor == "poly") {
          int number = std::round(arma::median(xvec(arma::find_finite(xvec))));
          x.X.col(i).replace(arma::datum::nan, number);
        } else {
          double number = arma::median(xvec(arma::find_finite(xvec)));
          x.X.col(i).replace(arma::datum::nan, number);
        }
      }
    } else if(x.missing == "complete.cases") {
      std::vector<int> indices;
      for (int i = 0; i < x.nobs; ++i) {
        // Check if the row contains any NaN values
        if(x.X.row(i).has_nan()) {
          indices.push_back(i);
        }
      }
      arma::uvec remove = arma::conv_to<arma::uvec>::from(indices);
      x.X.shed_rows(remove);
    } else if(x.missing == "pairwise.complete.cases") {
      return;
    } else {
      Rcpp::stop("Please, provide a method to handle missing data via the 'missing' argument");
    }

  }

}

double root_quad(double a, double b, double c) {

  // Solution to quadratic interpolation:

  double res = 0.5 * (- b + sqrt(b * b - 4 * a * c) ) / a;

  return res;

}

arma::mat diagg(int q) {

  arma::mat D(q, q, arma::fill::eye);
  arma::mat H(q, q*q, arma::fill::zeros);
  arma::uvec indexes = arma::find(D == 1);
  for(int i=0; i < q; ++i) H(i, indexes[i]) = 1;

  return H;

}

arma::mat zeros(arma::mat X, std::vector<arma::uvec> indexes) {

  int I = indexes.size();
  arma::mat X0 = X;

  for(int i=0; i < I; ++i) {

    X0(indexes[i], indexes[i]).zeros();

  }

  X0.diag() = X.diag();

  return X0;

}

arma::vec new_k(std::vector<std::string> x, std::string y, arma::vec k) {

  // Resize gamma, k and epsilon

  int k_size = k.size();
  int x_size = x.size();

  std::vector<int> id; // Setup storage for found IDs

  for(int i =0; i < x.size(); i++) // Loop through input
    if(x[i] == y) {// check if input matches target
      id.push_back(i);
    }

    arma::uvec indexes = arma::conv_to<arma::uvec>::from(id);
    std::vector<double> ks = arma::conv_to<std::vector<double>>::from(k);

    if(k_size < x_size) {
      ks.resize(x_size, k.back());
      k = arma::conv_to<arma::vec>::from(ks);
      k(indexes) = k(arma::span(0, indexes.size()-1));
    }

    return k; // send locations to R (c++ index shift!)
}

arma::mat tc(int g) {

  // Generate a matrix with all the combination of triads

  if(g == 0) {
    arma::mat Ng; return Ng;
  } else if(g == 1) {
    arma::mat Ng(1, 1, arma::fill::ones); return Ng;
  } else if(g == 2) {
    arma::mat Ng(2, 1, arma::fill::ones); return Ng;
  }

  int k = g*(g-1)*(g-2)/6;
  arma::mat Ng(g, k);

  int i = 0;

  for(int k=0; k < (g-2); ++k) {
    for(int j=k+1; j < (g-1); ++j) {
      for(int l=j+1; l < g; ++l) {

        Ng(k, i) = 1;
        Ng(j, i) = 1;
        Ng(l, i) = 1;
        ++i;

      }
    }
  }

  return Ng;

}

// [[Rcpp::export]]
arma::mat smoothing(arma::mat X, double min_eigval = 0.001) {

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, X);
  arma::vec new_eigval = arma::clamp(eigval, min_eigval, arma::datum::inf);

  arma::mat smoothed = eigvec * arma::diagmat(new_eigval) * eigvec.t();
  arma::mat diag = arma::diagmat(1/arma::sqrt(arma::diagvec(smoothed)));

  return diag * smoothed * diag;

}

arma::mat dxt(int p, int q) {

  /*
   * derivative of a matrix wrt its transpose
   */

  int pq = p*q;

  arma::mat res(pq, pq);
  arma::mat temp(p, q);

  for(int i=0; i < pq; ++i) {
    temp.zeros();
    temp(i) = 1;
    res.col(i) = arma::vectorise(temp.t(), 0);
  }

  return res;

}

arma::mat skew(arma::mat X) {

  // Skew-symmetric matrix

  return 0.5 * (X - X.t());

}

arma::mat symm(arma::mat X) {

  // Symmetric matrix

  return 0.5 * (X + X.t());

}

arma::mat lyap_sym(arma::mat Y, arma::mat Q) {

  // Solve the lyapunov equation YX + XY = Q with symmetric Q and X:

  int q = Y.n_cols;
  arma::vec I(q, arma::fill::ones);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Y);

  arma::mat M = eigvec.t() * Q * eigvec;
  arma::mat W1 = I * eigval.t();
  arma::mat W = W1 + W1.t();
  arma::mat YY = M / W;
  arma::mat A = eigvec * YY * eigvec.t();

  return A;

}

// [[Rcpp::export]]
arma::uvec consecutive(int lower, int upper) {

  // Generate a sequence of integers from lower to upper

  int size = upper - lower + 1;
  arma::uvec ivec(size);
  std::iota(ivec.begin(), ivec.end(), lower);

  return ivec;
}

std::vector<arma::uvec> vector_to_list(arma::uvec v){

  // Pass a vector to a list

  int n = v.size();
  std::vector<arma::uvec> lista(n);
  v.insert_rows(0, 1);

  for(int i=0; i < n; ++i) {

    lista[i] = v[i] + consecutive(1, v[i+1]);

  }

  return lista;

}

std::vector<arma::uvec> vector_to_list2(arma::uvec v){

  // Pass a vector to a list of sequential vectors

  int n = v.size();
  int add = 0;
  std::vector<arma::uvec> lista(n);

  for(int i=0; i < n; ++i) {

    if(i != 0) {
      add = lista[i-1].back() + 1;
    }

    lista[i] = add + consecutive(1, v[i]) - 1;

  }

  return lista;

}

std::vector<arma::uvec> vector_to_list3(arma::uvec v){

  // Pass a vector to a list of sequential vectors

  int n = v.size();
  int add = 0;
  std::vector<arma::uvec> lista(n);

  for(int i=0; i < n; ++i) {

    if(i != 0) {
      add = lista[i-1].back();
    }

    lista[i] = consecutive(1, v[i]);

  }

  return lista;

}

std::vector<arma::uvec> vector_to_list4(arma::uvec v){

  // Pass a vector to a list of sequential vectors

  int n = v.size();
  int add = 0;
  std::vector<arma::uvec> lista(n);

  for(int i=0; i < n; ++i) {

    if(i != 0) {
      add = lista[i-1].back();
    }

    lista[i] = add + consecutive(1, v[i]);

  }

  return lista;

}

// [[Rcpp::export]]
arma::vec orthogonalize(arma::mat X, arma::vec x) {

  // Make x orthogonal to every column of X (X must be orthogonal)

  for(int i=0; i < X.n_cols; ++i) {

    // x -= arma::accu(X.col(i) % x) / arma::accu(X.col(i) % X.col(i)) * X.col(i);
    x -= arma::accu(X.col(i) % x) * X.col(i);

  }

  // x /= sqrt(arma::accu(x % x));

  return x;

}

arma::uvec list_to_vector(std::vector<arma::uvec> X) {

  // Unlist to a vector

  arma::uvec single_vector = std::accumulate(X.begin(), X.end(),
                                             arma::uvec(), [](arma::uvec a, arma::uvec b) {
                                               a = arma::join_cols(a, b);
                                               return a;
                                             });

  return single_vector;

}

arma::uvec cum_indices(std::vector<arma::uvec> x) {

  int n = x.size();
  for(int i=1; i < n; ++i) {
    arma::uvec z = x[i-1L];
    x[i] += z.max() + 1L;
  }
  arma::uvec y = list_to_vector(x);
  return y;
}

std::vector<arma::uvec> increment(arma::uvec oblq_indexes, int p) {

  arma::uvec oblq_indexes_total = oblq_indexes;
  int n_blocks = oblq_indexes.size();
  int total = arma::accu(oblq_indexes);
  if(p != total) {
    oblq_indexes_total.insert_rows(n_blocks, 1);
    oblq_indexes_total[n_blocks] = (p - total + 0.00);
  }
  std::vector<arma::uvec> indexes_list = vector_to_list2(oblq_indexes_total);

  return indexes_list;

}

std::vector<int> subvector(std::vector<int> v, int lower, int upper) {

  std::vector<int> subv(&v[lower], &v[upper]);

  return subv;

}

std::vector<std::vector<int>> subvectors(std::vector<std::vector<int>> v, int lower, int upper) {

  std::vector<std::vector<int>> subv(&v[lower], &v[upper]);

  return subv;

}

// Gram-Schmidt process:

arma::mat gram(arma::mat X) {

  int n = X.n_rows;
  int k = X.n_cols;
  X.col(0) /= sqrt(arma::accu(X.col(0) % X.col(0)));

  for(int i=1; i < k; ++i) {

    for(int j=0; j < i; ++j) {

      X.col(i) -= arma::accu(X.col(j) % X.col(i)) / arma::accu(X.col(j) % X.col(j)) * X.col(j);

    }

    X.col(i) /= sqrt(arma::accu(X.col(i) % X.col(i)));

  }

  return X;

}

arma::mat kdiag(arma::mat X) {

  /*
   * Transform every column into a diagonal matrix and bind the results
   */

  int pq = X.n_rows;
  int q = X.n_cols;
  int p = pq/q;

  arma::mat res2(pq, 0);

  for(int j=0; j < q; ++j) {

    arma::mat res1(0, p);

    for(int i=0; i < q; ++i) {
      int index_1 = i*p;
      int index_2 = index_1 + (p-1);
      arma::mat temp = arma::diagmat(X(arma::span(index_1, index_2), j));
      res1 = arma::join_cols(res1, temp);
    }

    res2 = arma::join_rows(res2, res1);

  }

  return res2;

}

arma::mat cbind_diag(arma::mat X) {

  /*
   * Transform every column into a diagonal matrix and bind
   */

  int p = X.n_rows;
  int q = X.n_cols;
  arma::mat res(p, 0);

  for(int i=0; i < q; ++i) {
    res = arma::join_rows(res, arma::diagmat(X.col(i)));
  }

  return res;

}

arma::mat bc(int g) {

  // Generate a matrix with all the combination of pairs

  int k = g*(g-1)/2;
  arma::mat Ng(g, k);

  int i = 0;

  for(int k=0; k < (g-1); ++k) {
    for(int j=k+1; j < g; ++j) {

      Ng(k, i) = 1;
      Ng(j, i) = 1;
      ++i;

    }
  }

  return Ng;

}

arma::vec tucker_congruence(arma::mat X, arma::mat Y) {

  arma::vec YX = diagvec(Y.t() * X);
  arma::vec YY = diagvec(Y.t() * Y);
  arma::vec XX = diagvec(X.t() * X);

  // arma::vec YX = arma::sum(Y % X, 0);
  // arma::vec YY = arma::sum(Y % Y, 0);
  // arma::vec XX = arma::sum(X % X, 0);

  arma::vec congruence = YX / arma::sqrt(YY % XX);

  return arma::abs(congruence);

}

bool is_duplicate(arma::cube Targets, arma::mat Target, int length) {

  for(int i=length; i > -1; --i) {

    if(arma::approx_equal(Targets.slice(i), Target, "absdiff", 0)) return true;

  }

  return false;

}

void pass_to_efast(Rcpp::List efa_args, arguments_efast& x) {

  if (efa_args.containsElementNamed("estimator")) {
    std::string estimator_ = efa_args["estimator"]; x.estimator = estimator_;
  }
  if(efa_args.containsElementNamed("rotation")) {
    std::vector<std::string> rotation_ = efa_args["rotation"]; x.rotation = rotation_;
  }
  if(efa_args.containsElementNamed("projection")) {
    std::string projection_ = efa_args["projection"]; x.projection = projection_;
  }
  if(efa_args.containsElementNamed("nobs")) {
    Rcpp::Nullable<int> nobs_ = efa_args["nobs"]; x.nullable_nobs = nobs_;
  }
  if(efa_args.containsElementNamed("init")) {
    Rcpp::Nullable<arma::vec> init_ = efa_args["init"]; x.nullable_init = init_;
  }
  if (efa_args.containsElementNamed("Target")) {
    Rcpp::Nullable<arma::mat> Target_ = efa_args["Target"]; x.nullable_Target = Target_;
  }
  if (efa_args.containsElementNamed("Weight")) {
    Rcpp::Nullable<arma::mat> Weight_ = efa_args["Weight"]; x.nullable_Weight = Weight_;
  }
  if (efa_args.containsElementNamed("PhiTarget")) {
    Rcpp::Nullable<arma::mat> PhiTarget_ = efa_args["PhiTarget"]; x.nullable_PhiTarget = PhiTarget_;
  }
  if (efa_args.containsElementNamed("PhiWeight")) {
    Rcpp::Nullable<arma::mat> PhiWeight_ = efa_args["PhiWeight"]; x.nullable_PhiWeight = PhiWeight_;
  }
  if (efa_args.containsElementNamed("blocks")) {
    Rcpp::Nullable<std::vector<std::vector<arma::uvec>>> blocks_ = efa_args["blocks"]; x.nullable_blocks = blocks_;
  }
  if (efa_args.containsElementNamed("block_weights")) {
    Rcpp::Nullable<arma::vec> block_weights_ = efa_args["block_weights"]; x.nullable_block_weights = block_weights_;
  }
  if (efa_args.containsElementNamed("oblq_factors")) {
    Rcpp::Nullable<arma::uvec> oblq_factors_ = efa_args["oblq_factors"]; x.nullable_oblq_factors = oblq_factors_;
  }
  if (efa_args.containsElementNamed("gamma")) {
    arma::vec gamma_ = efa_args["gamma"]; x.gamma = gamma_;
  }
  if (efa_args.containsElementNamed("epsilon")) {
    arma::vec epsilon_ = efa_args["epsilon"]; x.epsilon = epsilon_;
  }
  if (efa_args.containsElementNamed("k")) {
    arma::vec k_ = efa_args["k"]; x.k = k_;
  }
  if (efa_args.containsElementNamed("w")) {
    double w_ = efa_args["w"]; x.w = w_;
  }
  if (efa_args.containsElementNamed("random_starts")) {
    int random_starts_ = efa_args["random_starts"]; x.random_starts = random_starts_;
  }
  if (efa_args.containsElementNamed("cores")) {
    int cores_ = efa_args["cores"]; x.cores = cores_;
  }
  if (efa_args.containsElementNamed("efa_control")) {
    Rcpp::Nullable<Rcpp::List> efa_control_ = efa_args["efa_control"]; x.nullable_efa_control = efa_control_;
  }
  if (efa_args.containsElementNamed("rot_control")) {
    Rcpp::Nullable<Rcpp::List> rot_control_ = efa_args["rot_control"]; x.nullable_rot_control = rot_control_;
  }

}

// Functions for updating the target matrix:

arma::mat get_target(arma::mat L, arma::mat Phi, double cutoff) {

  int I = L.n_rows;
  int J = L.n_cols;

  L.elem( arma::find_nonfinite(L) ).zeros();
  arma::mat loadings = L;

  if(cutoff > 0) {

    arma::mat A(I, J, arma::fill::ones);
    A.elem( find(abs(L) <= cutoff) ).zeros();
    return A;

  }

  /*
   * Find the squared normalized loadings.
   */

  arma::vec sqrt_communalities = sqrt(diagvec(L * Phi * L.t()));
  arma::mat norm_loadings = loadings;
  norm_loadings.each_col() /= sqrt_communalities;
  norm_loadings = pow(norm_loadings, 2);

  /*
   * Sort the squared normalized loadings by column in increasing direction and
   * compute the mean of the adyacent differences (DIFFs)
   */

  arma::mat sorted_norm_loadings = sort(norm_loadings);
  arma::mat diff_sorted_norm_loadings = diff(sorted_norm_loadings);
  arma::mat diff_means = mean(diff_sorted_norm_loadings, 0);

  // return diff_means;
  /*
   * Sort the absolute loading values by column in increasing direction and
   * find the column loading cutpoints (the smallest loading which DIFF is above the average)
   */

  arma::mat sorted_loadings = sort(abs(loadings));
  arma::vec cuts(J);

  for(int j=0; j < J; ++j) {
    for(int i=0; i < I; ++i) {
      if (diff_sorted_norm_loadings(i, j) >= diff_means(j)) {
        cuts(j) = sorted_loadings(i, j);
        // cuts(j) = sorted_norm_loadings(i, j);
        break;
      }
    }
  }

  /*
   * Create a target matrix inserting ones where squared normalized loadings are
   *  above the cutpoint
   */

  arma::mat Target(I, J, arma::fill::zeros);
  for(int j=0; j < J; ++j) {
    for(int i=0; i < I; ++i) {

      if(norm_loadings(i, j) > cuts(j)) {
        Target(i, j) = 1;
      }

    }
  }

  // return Target;

  arma::mat Target2 = Target;

  /*
   * check conditions C1 C2 C3
   */

  /*
   * C2
   * Replicate the loading matrix but with overall positive factors
   * Create submatrices for each column where the rows are 0
   * Check the rank of these submatrices
   */

  arma::mat multiplier = L;
  arma::mat a(1, J);
  double full_rank = J-1;

  for (int j=0; j < J; ++j) {

    if (mean(L.col(j)) < 0) {
      multiplier.col(j) = -L.col(j);
    }

    int size = I - accu(Target2.col(j)); // Number of 0s in column j

    arma::mat m(size, J); // submatrix of 0s in column j

    int p = 0;
    for(int i=0; i < I; ++i) {
      if(Target2(i, j) == 0) {
        m.row(p) = Target2.row(i);
        p = p+1;
      }
    }
    m.shed_col(j);

    double r = arma::rank(m);

    a(0, j) = r;
  }

  double condition = accu(full_rank - a);

  if (condition == 0) { // if all submatrices are of full rank

    return Target;

  } else {

    // Rcpp::Rcout << "Solution might not be identified" << std::endl;

    // indices de a que indican que las filas de m no son linealmente independientes o el numero de filas de m es inferior a J-1:
    int size = 0;
    for(int j=0; j < J; ++j) {
      if (a(0, j) != full_rank) {
        size = size+1;
      }
    }

    arma::uvec c(size);

    int p = 0;
    for(int j=0; j < J; ++j) {
      if (a(0, j) != full_rank) {
        c(p) = j;
        p = p+1;
      }
    }

    int h = 1;
    // Targ2[Targ2 == 0] <- NA
    for(int i=0; i < I; ++i) {
      for(int j=0; j < J; ++j) {
        if (Target2(i, j) == 0) {
          Target2(i, j) = arma::datum::nan;
        }
      }
    }

    for (int i=0; i < c.size(); ++i) {

      int h = c(i);
      // Targ2[which.min(as.matrix(multiplier[, h + 1]) * Targ2[, h]),h] <- NA
      arma::uword min_index = arma::index_min(multiplier.col(h) % Target2.col(h));
      Target2(min_index, h) = arma::datum::nan;

      // m <- Targ2[which(is.na(Targ2[, h])), -h]
      arma::uvec indexes = arma::find_nonfinite(Target2.col(h));
      arma::mat m(indexes.size(), J);

      for(int k=0; k < indexes.size(); ++k) {
        m.row(k) = Target2.row(indexes(k));
      }
      m.shed_col(h);

      // m[which(is.na(m))] <- 0
      m.elem( arma::find_nonfinite(m) ).zeros();

    }

    // Targ2[is.na(Targ2)] <- 0
    Target2.elem( arma::find_nonfinite(Target2) ).zeros();
    // Targ2[Targ2 == 1] <- NA
    // Targ[, 2:ncol(Targ)] <- Targ2

    return Target2;
  }

}

arma::mat get_target(arma::mat loadings, Rcpp::Nullable<arma::mat> nullable_Phi, double cutoff) {

  int J = loadings.n_cols;

  arma::mat Phi(J, J);

  if(nullable_Phi.isNotNull()) {
    Phi = Rcpp::as<arma::mat>(nullable_Phi);
  } else {
    Phi.eye();
  }

  return get_target(loadings, Phi, cutoff);

}

void update_target(int n_generals, int n, int nfactors,
                   arma::mat loadings, arma::mat Phi, double cutoff,
                   arma::mat& new_Target) {

  if(n_generals == 1) {

    loadings.shed_col(0);
    Phi.shed_col(0);
    Phi.shed_row(0);
    new_Target = get_target(loadings, Phi, cutoff);
    arma::vec add(n, arma::fill::ones);

    new_Target.insert_cols(0, add);

  } else {

    arma::mat loadings_g = loadings(arma::span::all, arma::span(0, n_generals-1));
    arma::mat loadings_s = loadings(arma::span::all, arma::span(n_generals, nfactors-1));

    arma::mat Phi_g = Phi(arma::span(0, n_generals-1), arma::span(0, n_generals-1));
    arma::mat Phi_s = Phi(arma::span(n_generals, nfactors-1), arma::span(n_generals, nfactors-1));

    arma::mat new_Target_g = get_target(loadings_g, Phi_g, cutoff);
    arma::mat new_Target_s = get_target(loadings_s, Phi_s, cutoff);

    new_Target = join_rows(new_Target_g, new_Target_s);

  }

}

// Derivatives wrt model correlation

arma::mat gLRhat(arma::mat Lambda, arma::mat Phi) {

  int p = Lambda.n_rows;
  int q = Lambda.n_cols;
  arma::mat I(p, p, arma::fill::eye);
  arma::mat LP = Lambda * Phi;
  arma::mat g1 = arma::kron(LP, I);
  arma::mat g21 = arma::kron(I, LP);
  arma::mat g2 = g21 * dxt(p, q);
  arma::mat g = g1 + g2;

  return g;

}

arma::mat gPRhat(arma::mat Lambda, arma::mat Phi, arma::uvec indexes) {

  int q = Phi.n_cols;
  // arma::uvec indexes = trimatl_ind(arma::size(Phi), -1);
  arma::mat g1 = arma::kron(Lambda, Lambda);
  arma::mat g2 = g1 * dxt(q, q);
  arma::mat g_temp = g1 + g2;
  arma::uvec indexes_diag_q(q);
  for(int i=0; i < q; ++i) indexes_diag_q[i] = i * q + i;
  g_temp.cols(indexes_diag_q) *= 0.5;
  // arma::mat g = g_temp.cols(indexes);

  return g_temp;
}

arma::mat gURhat(arma::mat Psi) {

  int p = Psi.n_cols;
  arma::mat gPsi = dxt(p, p);
  gPsi.diag().ones();
  // arma::uvec indexes = arma::trimatl_ind(arma::size(Psi), 0);

  // return gPsi.cols(indexes);
  return gPsi;

}
