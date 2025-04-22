/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Latent class analysis (normal) Probabilities in the raw scale
 */

class lca_gaussian: public estimators {

public:

  arma::mat Y; // matrix of response patterns (lowest response should be 0)
  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  arma::uvec K; // Vector of number of categories by item
  int nclasses; // Number of latent classes
  // arma::vec n; // number of subjects in each pattern
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

    classes(indices_target_classes) = parameters(indices_classes);
    logclasses = trunc_log(classes);

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        arma::uvec indices1 = indices_conditionals[j][i];
        arma::uvec indices2 = indices_target_conditionals[j][i];
        arma::uvec col = arma::uvec{i};
        conditionals[j](indices2, col) = parameters(indices1);
        for(int s=0; s < S; ++s) {
          double value = Y(s, j);
          double mu = conditionals[j](0, i);
          double sigma = conditionals[j](1, i);
          loglik(j, i, s) = logdnorm2(value, mu, sigma);
          latentloglik(s, i) += loglik(j, i, s);
        }
      }
    }

    latentpars = logclasses;

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

    arma::vec gclasses(nclasses, arma::fill::zeros);
    gmeans.set_size(J, nclasses); gmeans.zeros();
    gsds.set_size(J, nclasses); gsds.zeros();

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + latentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      gclasses += -n[s]*arma::trunc_exp(-loglik_case[s] + latentloglik.row(s).t());
      arma::vec term = -loglik_case[s] + jointlogp.row(s).t();
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          double constant = -n[s]*arma::trunc_exp(term[i] - loglik(j, i, s));
          double value = Y(s, j);
          double mu = conditionals[j](0, i);
          double sigma = conditionals[j](1, i);
          arma::vec gdnorm = ddnorm2(value, mu, sigma);
          gmeans(j, i) += constant*gdnorm[0];
          gsds(j, i) += constant*gdnorm[1];
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

    //     arma::vec gmeans_v = arma::vectorise(gmeans);
    //     arma::vec glog_sds_v = arma::vectorise(glog_sds);
    //     arma::vec gz = arma::join_cols(gmeans_v, glog_sds_v);

    arma::vec gugz(parameters.n_elem, arma::fill::zeros);
    // gugz(indices_classes) += gclasses(indices_target_classes);
    gugz(indices_conditionals2) += gz(indices_target_conditionals2);
    g = gugz;
    // g.replace(arma::datum::nan, 0);
    g.elem( arma::find_nonfinite(g) ).zeros();

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
    arma::vec post_total = arma::sum(freqs, 0).t();
    classes = post_total / arma::accu(post_total);
    logclasses = trunc_log(classes);

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        conditionals[j](0, i) = arma::accu(Y.col(j) % freqs.col(i)) / post_total(i);
        arma::vec xdiff = Y.col(j) - conditionals[j](0, i);
        conditionals[j](1, i) = arma::accu(xdiff % xdiff % freqs.col(i)) / post_total(i);
        double mu = conditionals[j](0, i);
        double sigma2 = conditionals[j](1, i);
        for(int s=0; s < S; ++s) {
          double value = Y(s, j);
          loglik(j, i, s) = logDnorm(value, mu, sigma2);
          latentloglik(s, i) += loglik(j, i, s);
        }
      }
    }

    // for(int s=0; s < S; ++s) {
    //   latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
    // }

    arma::vec cond_vector;
    for (const auto& mat : conditionals) {
      // Flatten each matrix as a column vector and join everything
      cond_vector = arma::join_cols(cond_vector, arma::vectorise(mat));
    }
    parameters = arma::vectorise(arma::join_cols(classes, cond_vector));

    latentpars = logclasses;
    // Rprintf("Gaussian nparam = %zu\n", parameters.n_elem);
    // Rprintf("Number of Gaussian indices = %zu\n", indices.n_elem);

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    for(int s=0; s < S; ++s) {
      latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
      pYXX.row(s) = arma::trunc_exp(latentloglik.row(s) + logclasses.t()); // P(data |X = c) P(X = c)
      mid[s] = arma::accu(pYXX.row(s)); // P(data)
    }

    posterior = pYXX;
    posterior.each_col() /= mid; //P(X = c | data) = P(data | X = c) P(X = c) / P(data)
    loglik_case = arma::trunc_log(mid);

    vectors.resize(3);
    vectors[0] = classes;
    vectors[1] = -loglik_case;
    vectors[2] = n;

    matrices.resize(2);
    matrices[0] = posterior;
    matrices[1] = latentloglik;

    list_matrices.resize(1);
    list_matrices[0] = conditionals;

  }

};
