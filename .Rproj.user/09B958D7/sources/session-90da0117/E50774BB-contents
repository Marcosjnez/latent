/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 03/02/2025
 */

/*
 * Latent class analysis (multinomial) Probabilities in the raw scale
 */

class lca_multinomial: public estimators {

public:

  arma::mat Y; // matrix of response patterns (lowest response should be 0)
  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  arma::uvec K; // Vector of number of categories by item
  int nclasses; // Number of latent classes
  // arma::vec n; // number of subjects in each pattern
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

    classes(indices_target_classes) = parameters(indices_classes);
    logclasses = trunc_log(classes);

    for(int j=0; j < J; ++j) {
      for(unsigned int i=0; i < nclasses; ++i) {
        arma::uvec indices1 = indices_conditionals[j][i];
        arma::uvec indices2 = indices_target_conditionals[j][i];
        arma::uvec col = arma::uvec{i};
        conditionals[j](indices2, col) = parameters(indices1);
        logconditionals[j].col(i) = arma::trunc_log(conditionals[j].col(i));
        for(int s=0; s < S; ++s) {
          int category = Y(s, j);
          loglik(j, i, s) = logconditionals[j](category, i);
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
    // f = -arma::accu(logliks);
    f = 0.00;

  }

  void G() {

    arma::vec gclasses(nclasses, arma::fill::zeros);
    for(int j=0; j < J; ++j) {
      gconditionals[j].zeros();
    }

    for(int s=0; s < S; ++s) {
      jointlogp.row(s) = latentloglik.row(s) + latentpars.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      gclasses += -n[s]*arma::trunc_exp(-loglik_case[s] + latentloglik.row(s).t());
      arma::vec term = -loglik_case[s] + jointlogp.row(s).t();
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          int category = Y(s, j);
          gconditionals[j](category, i) += -n[s]*arma::trunc_exp(term[i] - loglik(j, i, s));
        }
      }
    }

    arma::vec gu = gclasses;
    arma::vec gz;
    for (const auto& mat : gconditionals) {
      // Flatten each matrix as a column vector and join everything
      gz = arma::join_cols(gz, arma::vectorise(mat));
    }

    arma::vec gugz(parameters.n_elem, arma::fill::zeros);
    // gugz(indices_classes) += gu(indices_target_classes);
    gugz(indices_conditionals2) += gz(indices_target_conditionals2);
    g = gugz;
    // g.replace(arma::datum::nan, 0);
    g.elem( arma::find_nonfinite(g) ).zeros();

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

    arma::mat freqs = posterior; // S x nclasses
    freqs.each_col() %= n;
    arma::vec post_total = arma::sum(freqs, 0).t();
    classes = post_total / arma::accu(post_total);
    logclasses = trunc_log(classes);

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
    parameters = arma::vectorise(arma::join_cols(classes, cond_vector));

    latentpars = logclasses;

    // Rprintf("Multinomial nparam = %zu\n", parameters.n_elem);
    // Rprintf("Number of Multinomial indices = %zu\n", indices.n_elem);

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    for(int s=0; s < S; ++s) {
      latentloglik.row(s) = arma::sum(loglik.slice(s), 0);
      pYXX.row(s) = arma::trunc_exp(latentloglik.row(s) + logclasses.t()); // P(data |X = c) P(X = c)
      mid[s] = arma::accu(pYXX.row(s)); // P(data)
    }

    posterior = pYXX; //P(X = c | data) = P(data | X = c) P(X = c) / P(data)
    posterior.each_col() /= mid;
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
