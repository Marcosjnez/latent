/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 03/02/2025
 */

/*
 * Latent class analysis (normal) Probabilities in the log scale
 */

class lca_gaussian_softmax: public estimators {

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
  arma::mat gmeans, glog_sds; // Provide this filled with zeroes
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

    // lclasses(indices_target_classes) = parameters(indices_classes);
    // classes = softmax(lclasses, 1.00);
    // logclasses = trunc_log(classes);

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        arma::uvec indices1 = indices_conditionals[j][i];
        arma::uvec indices2 = indices_target_conditionals[j][i];
        arma::uvec col = arma::uvec{i};
        conditionals[j](indices2, col) = parameters(indices1);
        arma::vec values = Y.col(j);
        double mu = conditionals[j](0, i);
        double logsigma = conditionals[j](1, i);
        loglik.tube(j, i) = logdnorm(values, mu, logsigma);
        latentloglik.col(i) += arma::vectorise(loglik.tube(j, i));
      }
    }

    // latentpars = logclasses;

  }

  void F() {

    f = 0.0;

  }

  void G() {

    gmeans.set_size(J, nclasses); gmeans.zeros();
    glog_sds.set_size(J, nclasses); glog_sds.zeros();

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + latentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      arma::vec term = -loglik_case[s] + jointlogp.row(s).t();
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          double constant = -n[s]*arma::trunc_exp(term[i] - loglik(j, i, s));
          double value = Y(s, j);
          double mu = conditionals[j](0, i);
          double logsigma = conditionals[j](1, i);
          arma::vec gdnorm = ddnorm(value, mu, logsigma);
          gmeans(j, i) += constant*gdnorm[0];
          glog_sds(j, i) += constant*gdnorm[1];
        }
      }
    }

    for(int j=0; j < J; ++j) {
      gconditionals[j].zeros();
      for(int i=0; i < nclasses; ++i) {
        gconditionals[j](0, i) = gmeans(j, i);
        gconditionals[j](1, i) = glog_sds(j, i);
      }
    }

    arma::vec gz;
    for (const auto& mat : gconditionals) {
      // Flatten each matrix as a column vector and join everything
      gz = arma::join_cols(gz, arma::vectorise(mat));
    }

    arma::vec gugz(parameters.n_elem, arma::fill::zeros);
    gugz(indices_conditionals2) += gz(indices_target_conditionals2);
    g = gugz;
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

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    for(int s=0; s < S; ++s) {
      latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
      pYXX.row(s) = arma::trunc_exp(latentloglik.row(s) + latentpars.t()); // P(data |X = c) P(X = c)
      mid[s] = arma::accu(pYXX.row(s)); // P(data)
    }

    posterior = pYXX;
    posterior.each_col() /= mid; //P(X = c | data) = P(data | X = c) P(X = c) / P(data)
    loglik_case = arma::trunc_log(mid);

    vectors.resize(3);
    vectors[0] = softmax(latentpars, 1.00);
    vectors[1] = -loglik_case;
    vectors[2] = n;

    matrices.resize(2);
    matrices[0] = posterior;
    matrices[1] = latentloglik;

    list_matrices.resize(1);
    list_matrices[0] = conditionals;

  }

};
