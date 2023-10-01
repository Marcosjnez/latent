// Optimization algorithms for rotation

class cor_optim {

public:

  virtual cor_NTR optim(arguments_cor x, cor_manifold *manifold,
                        cor_criterion *criterion) = 0;

};

// Riemannian gradient descent:

class cor_RGD:public cor_optim {

public:

  cor_NTR optim(arguments_cor x, cor_manifold *manifold,
                cor_criterion *criterion) {

    return cor_gd(x, manifold, criterion);

  }

};

// Riemannian Newton Trust-Region:

class cor_RNTR:public cor_optim {

public:

  cor_NTR optim(arguments_cor x, cor_manifold *manifold,
                cor_criterion *criterion) {

    return cor_ntr(x, manifold, criterion);

  }

};

// Select the optimization algorithm:

cor_optim* choose_cor_optim(std::string cor) {

  cor_optim* algorithm;
  if(cor == "gradient") {
    algorithm = new cor_RGD();
  } else if(cor == "newtonTR") {
    algorithm = new cor_RNTR();
  } else if(cor == "BFGS") {
  } else if(cor == "L-BFGS") {
  } else {

    Rcpp::stop("Available optimization rutines for rotation: \n 'gradient', 'BFGS', 'L-BFGS', 'newtonTR'. The default method is 'newtonTR'.");

  }

  return algorithm;

}

typedef std::tuple<arma::mat,
                   std::vector<std::vector<double>>,
                   std::vector<std::vector<double>>,
                   std::vector<std::vector<std::vector<int>>>,
                   arma::mat,
                   double,
                   arma::mat,
                   bool> polyfast_object;

// Rcpp::List poly2(const arma::mat& X, const int cores) {
//
//   /*
//    * Function to estimate the full polychoric correlation matrix
//    */
//
//   Rcpp::Timer timer;
//
//   Rcpp::List result;
//   const int n = X.n_rows;
//   const int q = X.n_cols;
//
//   arma::mat cor = arma::cor(X);
//   std::vector<std::vector<int>> cols(q);
//   std::vector<int> maxs(q);
//   std::vector<std::vector<double>> taus(q);
//   std::vector<size_t> s(q);
//   std::vector<std::vector<double>> mvphi(q);
//
//   omp_set_num_threads(cores);
// #pragma omp parallel for
//   for(size_t i = 0; i < q; ++i) {
//     cols[i] = arma::conv_to<std::vector<int>>::from(X.col(i));
//     maxs[i] = *max_element(cols[i].begin(), cols[i].end());
//     std::vector<int> frequencies = count(cols[i], n, maxs[i]-1L);
//     mvphi[i] = cumsum(frequencies);
//     taus[i] = mvphi[i]; // Cumulative frequencies
//     for (size_t j = 0; j < maxs[i]; ++j) {
//       mvphi[i][j] /= n;
//       taus[i][j] = Qnorm(mvphi[i][j]);
//     }
//     mvphi[i].push_back(1.0);
//     mvphi[i].insert(mvphi[i].begin(), 0.0);
//     taus[i].push_back(pos_inf);
//     taus[i].insert(taus[i].begin(), neg_inf);
//     s[i] = taus[i].size() -1L;
//   }
//
//   timer.step("thresholds");
//
//   int K = 0.5*q*(q-1);
//   std::vector<std::vector<std::vector<int>>> tabs(K);
//   int k = 0;
//
// #ifdef _OPENMP
//   omp_set_num_threads(cores);
// #pragma omp parallel for
// #endif
//   for(size_t i=0; i < (q-1L); ++i) {
//     for(int j=(i+1L); j < q; ++j) {
//       tabs[k] = joint_frequency_table(cols[i], n, maxs[i], cols[j], maxs[j]);
//       ++k;
//     }
//   }
//
//   timer.step("contingency_tables");
//
//   arguments_cor x;
//   x.nobs = n;
//   x.q = q;
//   x.taus = taus;
//   x.mvphi = mvphi;
//   x.s = s;
//   x.n = tabs;
//   x.gcor.set_size(q, q); x.gcor.zeros();
//   x.dgcor.set_size(q, q); x.dgcor.zeros();
//   x.n_pairs = 0.5*q*(q-1);
//
//   // Select one manifold:
//   cor_manifold* manifold = choose_cor_manifold("oblq");
//   // Select one specific criteria:
//   cor_criterion* criterion = choose_cor_criterion("poly");
//   // Select the optimization routine:
//   cor_optim* algorithm = choose_cor_optim("newtonTR");
//
//   // x.T = init;
//   // x.dT = dX;
//   // manifold->param(x);
//   // criterion->F(x);
//   // criterion->gcor(x);
//   // manifold->dcor(x);
//   // manifold->grad(x);
//   // criterion->dgcor(x);
//   // manifold->grad(x);
//   // manifold->dgrad(x);
//   // result["cor"] = x.cor;
//   // result["f"] = x.f;
//   // result["gcor"] = x.gcor;
//   // result["dcor"] = x.dcor;
//   // result["grad"] = x.g;
//   // result["dgcor"] = x.dgcor;
//   // result["dgrad"] = x.dg;
//   // result["taus"] = taus;
//   // result["mvphi"] = mvphi;
//   // result["tabs"] = tabs;
//   // return result;
//
//   x.T = random_oblq(q, q);
//   cor_NTR x2 = algorithm->optim(x, manifold, criterion);
//
//   timer.step("optimization");
//
//   arma::mat COR = std::get<0>(x2);
//   arma::mat T = std::get<1>(x2);
//   double obj = std::get<2>(x2);
//   int iters = std::get<3>(x2);
//   bool converge = std::get<4>(x2);
//
//   result["correlation"] = COR;
//   result["X"] = T;
//   result["f"] = obj;
//   result["convergence"] = converge;
//   result["iters"] = iters;
//   result["elapsed"] = timer;
//
//   return result;
//
// }

// double fpoly(const arma::mat& X, const arma::mat& S) {
//
//   /*
//    * Function to estimate the full polychoric correlation matrix
//    */
//
//   const int n = X.n_rows;
//   const int q = X.n_cols;
//
//   arma::mat cor = arma::cor(X);
//   std::vector<std::vector<int>> cols(q);
//   std::vector<int> maxs(q);
//   std::vector<std::vector<double>> taus(q);
//   std::vector<size_t> s(q);
//   std::vector<std::vector<double>> mvphi(q);
//
//   for(size_t i = 0; i < q; ++i) {
//     cols[i] = arma::conv_to<std::vector<int>>::from(X.col(i));
//     maxs[i] = *max_element(cols[i].begin(), cols[i].end());
//     std::vector<int> frequencies = count(cols[i], n, maxs[i]-1L);
//     mvphi[i] = cumsum(frequencies);
//     taus[i] = mvphi[i]; // Cumulative frequencies
//     for (size_t j = 0; j < maxs[i]; ++j) {
//       mvphi[i][j] /= n;
//       taus[i][j] = Qnorm(mvphi[i][j]);
//     }
//     mvphi[i].push_back(1.0);
//     mvphi[i].insert(mvphi[i].begin(), 0.0);
//     taus[i].push_back(pos_inf);
//     taus[i].insert(taus[i].begin(), neg_inf);
//     s[i] = taus[i].size() -1L;
//   }
//
//   int K = 0.5*q*(q-1);
//   std::vector<std::vector<std::vector<int>>> tabs(K);
//   int k = 0;
//
//   for(size_t i=0; i < (q-1L); ++i) {
//     for(int j=(i+1L); j < q; ++j) {
//       tabs[k] = joint_frequency_table(cols[i], n, maxs[i], cols[j], maxs[j]);
//       ++k;
//     }
//   }
//
//   arguments_cor x;
//   x.nobs = n;
//   x.q = q;
//   x.taus = taus;
//   x.mvphi = mvphi;
//   x.s = s;
//   x.n = tabs;
//   x.n_pairs = 0.5*q*(q-1);
//   x.cor = S;
//   // Select one specific criteria:
//   cor_criterion* criterion = choose_cor_criterion("poly");
//   criterion->F(x);
//
//   return x.f;
//
// }

// Rcpp::List old_poly(const arma::mat& X, const int cores) {
//
//   /*
//    * Function to estimate the full polychoric correlation matrix (fast)
//    */
//
//   Rcpp::Timer timer;
//
//   const int n = X.n_rows;
//   const int q = X.n_cols;
//
//   arma::mat cor = arma::cor(X);
//   std::vector<std::vector<int>> cols(q);
//   std::vector<int> maxs(q);
//   std::vector<std::vector<double>> taus(q);
//   std::vector<size_t> s(q);
//   std::vector<std::vector<double>> mvphi(q);
//
//   for(size_t i = 0; i < q; ++i) {
//     cols[i] = arma::conv_to<std::vector<int>>::from(X.col(i));
//     maxs[i] = *max_element(cols[i].begin(), cols[i].end());
//     std::vector<int> frequencies = count(cols[i], n, maxs[i]-1L);
//     mvphi[i] = cumsum(frequencies);
//     taus[i] = mvphi[i]; // Cumulative frequencies
//     for (size_t j = 0; j < maxs[i]; ++j) {
//       mvphi[i][j] /= n;
//       taus[i][j] = Qnorm(mvphi[i][j]);
//     }
//     mvphi[i].push_back(1.0);
//     mvphi[i].insert(mvphi[i].begin(), 0.0);
//     taus[i].push_back(pos_inf);
//     taus[i].insert(taus[i].begin(), neg_inf);
//     s[i] = taus[i].size() -1L;
//   }
//
//   timer.step("thresholds");
//
//   arma::mat polys(q, q, arma::fill::eye);
//   arma::mat iters(q, q, arma::fill::zeros);
//
// #ifdef _OPENMP
//   omp_set_num_threads(cores);
// #pragma omp parallel for
// #endif
//   for(size_t i=0; i < (q-1L); ++i) {
//     for(int j=(i+1L); j < q; ++j) {
//       std::vector<std::vector<int>> tab = joint_frequency_table(cols[i], n, maxs[i], cols[j], maxs[j]);
//       std::vector<double> rho = optimize(taus[i], taus[j], tab, s[i], s[j], mvphi[i], mvphi[j], n, cor(i, j));
//       polys(i, j) = polys(j, i) = rho[0];
//       iters(i, j) = iters(j, i) = rho[1];
//     }
//   }
//
//   timer.step("polychorics");
//
//   Rcpp::List result;
//   result["polychorics"] = polys;
//   result["thresholds"] = taus;
//   result["iters"] = iters;
//   result["elapsed"] = timer;
//
//   return result;
//
// }
