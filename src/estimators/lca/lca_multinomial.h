/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 16/05/2025
 */

/*
 * Latent class analysis (multinomial model)
 */

class lca_multinomial: public estimators {

public:

  arma::mat Y; // matrix of response patterns (lowest response should be 0)
  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  arma::uvec K; // Vector of number of categories by item
  int nclasses; // Number of latent classes
  arma::vec classes, logclasses; // log(P(X = c))
  std::vector<arma::mat> logconditionals; // Provide this filled with zeroes
  double constant_class;
  std::vector<arma::vec> constant_cond;
  std::vector<arma::mat> gconditionals; // Provide this filled with zeroes
  arma::uvec indices_classes, indices_target_classes;
  std::vector<std::vector<arma::uvec>> indices_conditionals, indices_target_conditionals;
  arma::cube loglik;
  arma::mat jointlogp;
  arma::vec logliks;
  arma::vec mid, loglik_case;
  arma::uvec targets;
  arma::uvec indices_conditionals2;
  arma::uvec indices_target_conditionals2;
  arma::mat pYXX;
  arma::mat logpYXX;
  // arma::vec PD;
  // arma::mat posterior;
  // arma::vec classes2;
  bool fix_logit;
  std::vector<arma::mat> conditionals2;

  void param() {

    latentloglik.zeros();

    classes(indices_target_classes) = transparameters(indices_classes);
    logclasses = trunc_log(classes);

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        arma::uvec indices1 = indices_conditionals[j][i];
        arma::uvec indices2 = indices_target_conditionals[j][i];
        arma::uvec col = arma::uvec{i};
        conditionals[j](indices2, col) = transparameters(indices1);
        logconditionals[j].col(i) = arma::trunc_log(conditionals[j].col(i));
        for(int s=0; s < S; ++s) {
          if(!std::isnan(Y(s, j))) {
            int category = Y(s, j);
            loglik(j, i, s) = logconditionals[j](category, i);
            latentloglik(s, i) += loglik(j, i, s);
          }
        }
      }
    }

    loglatentpars = logclasses;

  }

  void F() {

    // for(int s=0; s < S; ++s) {
    //   jointlogp.row(s) = latentloglik.row(s) + latentpars.t();
    //   double max_vector = jointlogp.row(s).max();
    //   loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
    //   logliks[s] = n[s] * loglik_case[s];
    // }
    //
    // f = -arma::accu(logliks);
    f = 0.00;

  }

  void G() {

    for(int j=0; j < J; ++j) {
      gconditionals[j].zeros();
    }

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + loglatentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      arma::vec log_posterior = -loglik_case[s] + jointlogp.row(s).t();
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          if(!std::isnan(Y(s, j))) {
            int category = Y(s, j);
            gconditionals[j](category, i) += -n[s]*arma::trunc_exp(log_posterior[i] - loglik(j, i, s));
          }
        }
      }
    }

    arma::vec gz;
    for (const auto& mat : gconditionals) {
      // Flatten each matrix as a column vector and join everything
      gz = arma::join_cols(gz, arma::vectorise(mat));
    }

    arma::vec gugz(transparameters.n_elem, arma::fill::zeros);
    gugz(indices_conditionals2) += gz(indices_target_conditionals2);
    grad = gugz;
    grad.elem( arma::find_nonfinite(grad) ).zeros();

  }

  void dG() {

    // Rcpp::stop("dG not available");
    dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H() {

    // Rcpp::stop("H not available");
    hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() { // Update the parameter estimates

    latentloglik.zeros();

    // Rprintf("Posterior:\n");
    // Rprintf("%zu ", posterior.n_rows);
    // Rprintf("%zu ", posterior.n_cols);
    // Rprintf("\n\n");
    //
    // Rprintf("n:\n");
    // Rprintf("%zu ", n.n_elem);
    // Rprintf("\n\n");

    arma::mat freqs = posterior; // S x nclasses
    freqs.each_col() %= n;
    arma::vec post_total = arma::sum(freqs, 0).t();

    for(int j=0; j < J; ++j) {
      conditionals[j].zeros();
      for(unsigned int i=0; i < nclasses; ++i) {
        for(int s=0; s < S; ++s) {
          int category = Y(s, j);
          conditionals[j](category, i) += freqs(s, i);
        }
      }
      conditionals[j].each_row() /= post_total.t();
    }

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        logconditionals[j].col(i) = arma::trunc_log(conditionals[j].col(i));
        for(int s=0; s < S; ++s) {
          int category = Y(s, j);
          loglik(j, i, s) = logconditionals[j](category, i);
          latentloglik(s, i) += loglik(j, i, s);
        }
      }
    }

    arma::vec cond_vector;
    for (const auto& mat : conditionals) {
      // Flatten each matrix as a column vector and join everything
      cond_vector = arma::join_cols(cond_vector, arma::vectorise(mat));
    }
    transparameters = arma::vectorise(arma::join_cols(latentpars, cond_vector));

    // Rf_error("170");
    // Rprintf("Multinomial nparam = %zu\n", parameters.n_elem);
    // Rprintf("Number of Multinomial indices = %zu\n", indices.n_elem);

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    for(int s=0; s < S; ++s) {
      latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
      pYXX.row(s) = arma::trunc_exp(latentloglik.row(s) + loglatentpars.t()); // P(data |X = c) P(X = c)
      mid[s] = arma::accu(pYXX.row(s)); // P(data)
    }

    posterior = pYXX; //P(X = c | data) = P(data | X = c) P(X = c) / P(data)
    posterior.each_col() /= mid;
    loglik_case = arma::trunc_log(mid);

    vectors.resize(4);
    vectors[0] = latentpars;
    vectors[1] = loglatentpars;
    vectors[2] = loglik_case;
    vectors[3] = n;

    matrices.resize(2);
    matrices[0] = posterior;
    matrices[1] = latentloglik;

    list_matrices.resize(2);
    list_matrices[0] = conditionals;
    list_matrices[1] = logconditionals;

  }

};

lca_multinomial* choose_lca_multinomial(Rcpp::List estimator_setup) {

  lca_multinomial* myestimator = new lca_multinomial();

  arma::mat Y = estimator_setup["Y"];
  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  arma::uvec K = estimator_setup["K"];
  int nclasses = estimator_setup["nclasses"];
  arma::vec n = estimator_setup["n"];
  arma::vec classes = estimator_setup["classes"];
  // arma::vec logclasses = estimator_setup["logclasses"];
  double constant_class = estimator_setup["constant_class"];
  bool fix_logit = estimator_setup["fix_logit"];
  std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
  // std::vector<arma::mat> logconditionals = estimator_setup["logconditionals"];
  std::vector<arma::vec> constant_cond = estimator_setup["constant_cond"];;
  // std::vector<arma::mat> gconditionals = estimator_setup["gconditionals"];
  arma::uvec indices_classes = estimator_setup["indices_classes"];
  arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
  std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
  std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
  arma::uvec indices = estimator_setup["indices"];
  // arma::uvec targets = estimator_setup["targets"];
  arma::uvec indices_conditionals2 = estimator_setup["indices_conditionals2"];
  arma::uvec indices_target_conditionals2 = estimator_setup["indices_target_conditionals2"];

  arma::cube loglik(J, nclasses, S, arma::fill::zeros);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec mid(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::mat pYXX(S, nclasses, arma::fill::zeros);
  // arma::mat posterior(S, nclasses, arma::fill::randu);
  // posterior = arma::normalise(posterior, 1);
  arma::mat latentloglik(S, nclasses, arma::fill::zeros);
  arma::mat jointlogp(S, nclasses, arma::fill::zeros);

  myestimator->Y = Y;
  myestimator->S = S;
  myestimator->J = J;
  myestimator->K = K;
  myestimator->nclasses = nclasses;
  myestimator->n = n;
  myestimator->classes = classes;
  myestimator->constant_class = constant_class;
  myestimator->fix_logit = fix_logit;
  myestimator->conditionals = conditionals;
  myestimator->constant_cond = constant_cond;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_target_classes = indices_target_classes;
  myestimator->indices_conditionals = indices_conditionals;
  myestimator->indices_target_conditionals = indices_target_conditionals;
  myestimator->indices = indices;
  myestimator->indices_conditionals2 = indices_conditionals2;
  myestimator->indices_target_conditionals2 = indices_target_conditionals2;

  myestimator->logconditionals = conditionals;
  myestimator->gconditionals = conditionals;
  myestimator->loglik = loglik;
  myestimator->logliks = logliks;
  myestimator->mid = mid;
  myestimator->loglik_case = loglik_case;
  myestimator->pYXX = pYXX;
  // myestimator->posterior = posterior;
  myestimator->latentloglik = latentloglik;
  myestimator->jointlogp = jointlogp;
  myestimator->evalEM = true;

 return myestimator;

}
