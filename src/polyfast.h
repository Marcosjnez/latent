/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 06/09/2026
 */

Rcpp::List polyfast(arma::mat X, std::string missing, const std::string smooth,
                    double min_eigval, const int nboot, const bool fit, const int cores) {

  Rcpp::Timer timer;
  const int q = X.n_cols;
  if(q < 1 || X.n_rows < 2 || cores < 1) {
    Rcpp::stop("polyfast requires at least two rows, one column and one thread.");
  }
  if(X.has_inf()) Rcpp::stop("Ordinal values must be finite or missing.");
  if(missing == "complete.obs" || missing == "complete.cases" || missing == "listwise") {
    const arma::uvec rows = arma::find(arma::sum(X != X, 1) == 0);
    X = X.rows(rows);
  } else if(missing == "all.obs") {
    if(X.has_nan()) Rcpp::stop("Missing ordinal values with missing = 'all.obs'.");
  } else if(missing != "pairwise.complete.cases" && missing != "pairwise.complete.obs") {
    Rcpp::stop("Use pairwise or complete observations for polyfast.");
  }
  const int n = X.n_rows;
  std::vector<std::vector<double>> levels(q), taus(q), cumulative(q);
  std::vector<std::vector<int>> codes(q, std::vector<int>(n, -1));
  std::vector<int> category_count(q);
  std::vector<double> category_min(q);

  for(int j = 0; j < q; ++j) {
    // Ordinal scores need not be consecutive integers.
    const arma::vec values = X.col(j);
    const arma::uvec observed = arma::find_finite(values);
    const arma::vec categories = arma::unique(values.elem(observed));
    const int K = categories.n_elem;
    if(K < 2) Rcpp::stop("Every variable needs at least two observed categories.");
    levels[j] = arma::conv_to<std::vector<double>>::from(categories);
    category_min[j] = categories[0];
    category_count[j] = K;

    std::vector<int> frequencies(K, 0);
    for(arma::uword row : observed) {
      int k = std::lower_bound(levels[j].begin(), levels[j].end(), values[row]) - levels[j].begin();
      codes[j][row] = k;
      ++frequencies[k];
    }
    taus[j].resize(K+1);
    cumulative[j].resize(K+1);
    taus[j].front() = neg_inf;
    taus[j].back() = pos_inf;
    cumulative[j].front() = 0.0;
    cumulative[j].back() = 1.0;
    int total = 0;
    for(int k = 1; k < K; ++k) {
      total += frequencies[k-1];
      cumulative[j][k] = static_cast<double>(total)/observed.n_elem;
      taus[j][k] = R::qnorm(cumulative[j][k], 0.0, 1.0, true, false);
    }
  }

  const int npairs = q*(q-1)/2;
  arma::mat correlation(q, q, arma::fill::eye), iters(q, q, arma::fill::zeros);
  arma::vec grad(npairs), hess(npairs);
  std::vector<std::vector<std::vector<int>>> tables(npairs);
  std::vector<int> pair_counts(npairs);

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for(int i = 0; i < q-1; ++i) {
    for(int j = i+1; j < q; ++j) {
      const int k = i*(2*q-i-1)/2 + j-i-1;
      auto& table = tables[k];
      table.assign(category_count[i], std::vector<int>(category_count[j], 0));
      int count = 0;
      for(int row = 0; row < n; ++row) {
        const int a = codes[i][row], b = codes[j][row];
        if(a >= 0 && b >= 0) {
          ++table[a][b];
          ++count;
        }
      }
      pair_counts[k] = count;
      if(count < 2) continue;
      std::vector<double> rho = optimize(taus[i], taus[j], table, category_count[i],
                                         category_count[j], cumulative[i], cumulative[j],
                                         count, 0.0);
      correlation(i, j) = correlation(j, i) = rho[0];
      iters(i, j) = iters(j, i) = rho[1];
      grad(k) = rho[2];
      hess(k) = rho[3];
    }
  }
  if(std::find_if(pair_counts.begin(), pair_counts.end(), [](int x) {return x < 2;}) != pair_counts.end()) {
    Rcpp::stop("Each pair of ordinal variables needs at least two joint observations.");
  }

  timer.step("polychorics");
  return Rcpp::List::create(Rcpp::_["type"] = "polychorics",
                            Rcpp::_["correlation"] = correlation,
                            Rcpp::_["thresholds"] = taus,
                            Rcpp::_["contingency_tables"] = tables,
                            Rcpp::_["cumulative_freqs"] = cumulative,
                            Rcpp::_["category_min"] = category_min,
                            Rcpp::_["category_count"] = category_count,
                            Rcpp::_["category_levels"] = levels,
                            Rcpp::_["iters"] = iters,
                            Rcpp::_["grad"] = grad,
                            Rcpp::_["hess"] = hess,
                            Rcpp::_["elapsed"] = timer);
}
