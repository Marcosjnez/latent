/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 16/05/2025
 */

/*
 * Latent class analysis (gaussian model)
 */

class lca_gaussian: public estimators {

public:

  arma::mat Y; // matrix of response patterns (lowest response should be 0)
  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  arma::uvec K; // Vector of number of categories by item
  int nclasses; // Number of latent classes
  arma::vec lclasses, classes, logclasses;
  std::vector<arma::mat> logconditionals; // Provide this filled with zeroes
  double constant_class;
  std::vector<arma::vec> constant_cond;
  arma::mat gmeans, gsds; // Provide this filled with zeroes
  std::vector<arma::mat> gconditionals;
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
        for(int s=0; s < S; ++s) {
          if(!std::isnan(Y(s, j))) {
            double value = Y(s, j);
            double mu = conditionals[j](0, i);
            double sigma = conditionals[j](1, i);
            loglik(j, i, s) = logdnorm2(value, mu, sigma);
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
    // logliks.elem( arma::find_nonfinite(logliks) ).zeros();
    // f = -arma::accu(logliks);
    f = 0.0;

  }

  void G() {

    gmeans.set_size(J, nclasses); gmeans.zeros();
    gsds.set_size(J, nclasses); gsds.zeros();

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + loglatentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      arma::vec term = -loglik_case[s] + jointlogp.row(s).t();
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          double constant = -n[s]*arma::trunc_exp(term[i] - loglik(j, i, s));
          if(!std::isnan(Y(s, j))) {
            double value = Y(s, j);
            double mu = conditionals[j](0, i);
            double sigma = conditionals[j](1, i);
            arma::vec gdnorm = ddnorm2(value, mu, sigma);
            gmeans(j, i) += constant*gdnorm[0];
            gsds(j, i) += constant*gdnorm[1];
          }
        }
      }
    }

    for(int j=0; j < J; ++j) {
      gconditionals[j].zeros();
      for(int i=0; i < nclasses; ++i) {
        gconditionals[j](0, i) = gmeans(j, i);
        gconditionals[j](1, i) = gsds(j, i);
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
    // dg.set_size(parameters.n_elem); dg.zeros();

  }

  void H() {

    // Rcpp::stop("H not available");
    // hessian.set_size(parameters.n_elem, parameters.n_elem); hessian.zeros();

  }

  void E() { // Update the parameter estimates and loglik of each observation

    latentloglik.zeros();

    arma::mat freqs = posterior; // S x nclasses
    freqs.each_col() %= n;
    arma::rowvec post_total = arma::sum(freqs, 0);
    freqs.each_row() /= post_total;

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        conditionals[j](0, i) = arma::accu(Y.col(j) % freqs.col(i));
        arma::vec xdiff = Y.col(j) - conditionals[j](0, i);
        conditionals[j](1, i) = arma::accu(xdiff % xdiff % freqs.col(i));
        double mu = conditionals[j](0, i);
        double sigma2 = conditionals[j](1, i);
        for(int s=0; s < S; ++s) {
          double value = Y(s, j);
          loglik(j, i, s) = logDnorm(value, mu, sigma2);
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

    // Rprintf("Gaussian nparam = %zu\n", parameters.n_elem);
    // Rprintf("Number of Gaussian indices = %zu\n", indices.n_elem);

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    for(int s=0; s < S; ++s) {
      latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
      pYXX.row(s) = arma::trunc_exp(latentloglik.row(s) + loglatentpars.t()); // P(data |X = c) P(X = c)
      mid[s] = arma::accu(pYXX.row(s)); // P(data)
    }

    posterior = pYXX;
    posterior.each_col() /= mid; //P(X = c | data) = P(data | X = c) P(X = c) / P(data)
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

lca_gaussian* choose_lca_gaussian(Rcpp::List estimator_setup) {

  lca_gaussian* myestimator = new lca_gaussian();

  arma::mat Y = estimator_setup["Y"];
  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  arma::vec n = estimator_setup["n"];
  int nclasses = estimator_setup["nclasses"];
  arma::vec classes = estimator_setup["classes"];
  double constant_class = estimator_setup["constant_class"];
  bool fix_logit = estimator_setup["fix_logit"];
  std::vector<arma::mat> conditionals = estimator_setup["conditionals"];
  arma::uvec indices_classes = estimator_setup["indices_classes"];
  arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
  std::vector<std::vector<arma::uvec>> indices_conditionals = estimator_setup["indices_conditionals"];
  std::vector<std::vector<arma::uvec>> indices_target_conditionals = estimator_setup["indices_target_conditionals"];
  arma::uvec indices = estimator_setup["indices"];
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
  myestimator->n = n;
  myestimator->nclasses = nclasses;
  myestimator->classes = classes;
  myestimator->constant_class = constant_class;
  myestimator->fix_logit = fix_logit;
  myestimator->conditionals = conditionals;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_target_classes = indices_target_classes;
  myestimator->indices_conditionals = indices_conditionals;
  myestimator->indices_target_conditionals = indices_target_conditionals;
  myestimator->indices = indices;
  myestimator->indices_conditionals2 = indices_conditionals2;
  myestimator->indices_target_conditionals2 = indices_target_conditionals2;

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
