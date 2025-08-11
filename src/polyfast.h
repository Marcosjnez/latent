/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 04/07/2025
 */

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

typedef std::tuple<arma::mat,
                   std::vector<std::vector<double>>,
                   std::vector<std::vector<double>>,
                   std::vector<std::vector<std::vector<int>>>,
                   arma::mat> polyfast_object;

polyfast_object poly(const arma::mat& X, const std::string smooth, double min_eigval,
                     const bool fit, const int cores) {

  /*
   * Function to estimate the full polychoric correlation matrix
   *                      (Multiple cores)
   */

  const int n = X.n_rows;
  const int q = X.n_cols;

  arma::mat cor;
  if(X.has_nan()) {
    cor = pairwise_cor(X);
  } else {
    cor = arma::cor(X);
  }

  std::vector<std::vector<int>> cols(q);
  std::vector<int> maxs(q);
  std::vector<int> mins(q);
  std::vector<std::vector<double>> taus(q);
  std::vector<size_t> s(q);
  std::vector<std::vector<double>> mvphi(q);
  arma::mat X2 = X;

  for(size_t i = 0; i < q; ++i) {

    mins[i] = X2.col(i).min();
    X2.col(i) -= mins[i];
    cols[i] = arma::conv_to<std::vector<int>>::from(X2.col(i));
    maxs[i] = *max_element(cols[i].begin(), cols[i].end());

    std::vector<int> frequencies = count(cols[i], n, maxs[i]);
    mvphi[i] = cumsum(frequencies); // Cumulative frequencies
    taus[i] = mvphi[i];

    for (size_t j = 0; j < maxs[i]; ++j) {
      mvphi[i][j] /= n;
      taus[i][j] = Qnorm(mvphi[i][j]);
    }
    mvphi[i].push_back(1.0);
    mvphi[i].insert(mvphi[i].begin(), 0.0);
    taus[i].push_back(pos_inf);
    taus[i].insert(taus[i].begin(), neg_inf);
    s[i] = taus[i].size() - 1L;
  }

  arma::mat polys(q, q, arma::fill::eye);
  arma::mat iters(q, q, arma::fill::zeros);

  int d = 0.5*q*(q-1);
  std::vector<std::vector<std::vector<int>>> tabs(d);
  arma::vec seq = arma::linspace(0, q-1, q);
  arma::vec I = arma::cumsum(q - seq) - 2*q;

#ifdef _OPENMP
  omp_set_num_threads(cores);
#pragma omp parallel for
#endif
  for(size_t i=0; i < (q-1L); ++i) {
    for(size_t j=(i+1L); j < q; ++j) {
      int k = I[i+1] + j;
      tabs[k] = joint_frequency_table(cols[i], n, maxs[i], cols[j], maxs[j]);
      std::vector<double> rho = optimize(taus[i], taus[j], tabs[k], s[i], s[j], mvphi[i], mvphi[j], n, cor(i, j));
      polys(i, j) = polys(j, i) = rho[0];
      iters(i, j) = iters(j, i) = rho[1];
    }
  }

  polyfast_object result = std::make_tuple(polys, taus, mvphi, tabs, iters);

  return result;

}

Rcpp::List polyfast(arma::mat X, std::string missing, std::string acov, const std::string smooth,
                    double min_eigval, const int nboot, const bool fit, const int cores) {

  /*
   * Function to estimate the full polychoric correlation matrix
   */

  Rcpp::Timer timer;

  // polyfast_object x = poly(xcor.X, smooth, min_eigval, false, cores);
  polyfast_object x = poly(X, smooth, min_eigval, false, cores);

  timer.step("polychorics");

  arma::mat polys = std::get<0>(x);
  std::vector<std::vector<double>> taus = std::get<1>(x);
  std::vector<std::vector<double>> mvphis = std::get<2>(x);
  std::vector<std::vector<std::vector<int>>> tabs = std::get<3>(x);
  arma::mat iters = std::get<4>(x);

  Rcpp::List result;
  result["type"] = "polychorics";
  result["correlation"] = polys;
  result["thresholds"] = taus;
  result["contingency_tables"] = tabs;
  result["cumulative_freqs"] = mvphis;
  result["iters"] = iters;
  result["elapsed"] = timer;

  return result;

}
