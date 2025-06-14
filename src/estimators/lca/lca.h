/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 12/06/2025
 */

/*
 * Latent class analysis
 */

class lca: public estimators {

public:

  arma::vec n; // Number of repetitions of each pattern
  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  arma::uvec K; // Vector of number of categories by item
  int nclasses; // Number of latent classes
  arma::vec classes, logclasses; // log(P(X = c))
  arma::uvec indices_classes, indices_target_classes;
  arma::uvec indices_items, indices_target_items;
  arma::mat jointlogp;
  arma::cube loglik, dloglik;
  arma::mat latentloglik;
  arma::vec logliks;
  arma::vec loglik_case;
  arma::mat posterior, logposterior;

  void param() {

    classes(indices_target_classes) = transparameters(indices_classes);
    logclasses = trunc_log(classes);
    loglik.elem(indices_target_items) = transparameters(indices_items);

    latentloglik.zeros();

    for(int s=0; s < S; ++s) {
      for(int i=0; i < nclasses; ++i) {
        for(int j=0; j < J; ++j) {
          latentloglik(s, i) += loglik(i, j, s);
        }
      }
      double max_vector = jointlogp.row(s).max();
      jointlogp.row(s) = latentloglik.row(s) + logclasses.t();
      loglik_case(s) = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      logposterior.row(s) = jointlogp.row(s) - loglik_case(s);
      logliks(s) = n(s) * loglik_case(s);
    }

  }

  void F() {

    // for(int s=0; s < S; ++s) {
    //   jointlogp.row(s) = latentloglik.row(s) + logclasses.t();
    //   double max_vector = jointlogp.row(s).max();
    //   loglik_case[s] = max_vector + arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
    //   logliks[s] = n[s] * loglik_case[s];
    // }

    f = -arma::accu(logliks);

  }

  void G() {

    arma::vec dclasses(nclasses, arma::fill::zeros);

    for(int s=0; s < S; ++s) {
      for(int i=0; i < nclasses; ++i) {
        dclasses(i) += -n[s]*arma::trunc_exp(latentloglik(s, i) - loglik_case[s]);
        for(int j=0; j < J; ++j) {
          dloglik(i, j, s) += -n[s]*arma::trunc_exp(logposterior(s, i) - loglik(i, j, s));
        }
      }
    }

    // Rprintf("\n\n");
    // Rprintf("grad\n");
    // Rprintf("%zu ", grad.n_elem);
    // Rprintf("\n\n");
    // Rf_error("87");
    grad.set_size(transparameters.n_elem);
    grad.elem(indices_classes) = dclasses;
    grad.elem(indices_items) = arma::vectorise(dloglik);

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

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    vectors.resize(4);
    vectors[0] = classes;
    vectors[1] = logclasses;
    vectors[2] = loglik_case;
    vectors[3] = n;

    matrices.resize(3);
    matrices[0] = loglik;
    matrices[1] = latentloglik;
    matrices[2] = posterior;

  }

};

lca* choose_lca(const Rcpp::List& estimator_setup) {

  lca* myestimator = new lca();

  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  int nclasses = estimator_setup["nclasses"];
  arma::vec n = estimator_setup["n"];
  arma::uvec indices = estimator_setup["indices"];
  arma::uvec indices_classes = estimator_setup["indices_classes"];
  arma::uvec indices_target_classes = estimator_setup["indices_target_classes"];
  arma::uvec indices_items = estimator_setup["indices_items"];
  arma::uvec indices_target_items = estimator_setup["indices_target_items"];

  arma::vec classes(nclasses);
  arma::cube loglik(nclasses, J, S, arma::fill::zeros);
  arma::cube dloglik(nclasses, J, S, arma::fill::zeros);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::mat latentloglik(S, nclasses, arma::fill::zeros);
  arma::mat jointlogp(S, nclasses, arma::fill::zeros);
  arma::mat posterior(S, nclasses);
  arma::mat logposterior(S, nclasses);

  myestimator->S = S;
  myestimator->J = J;
  myestimator->nclasses = nclasses;
  myestimator->n = n;
  myestimator->classes = classes;
  myestimator->indices = indices;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_target_classes = indices_target_classes;
  myestimator->indices_items = indices_items;
  myestimator->indices_target_items = indices_target_items;

  myestimator->loglik = loglik;
  myestimator->dloglik = dloglik;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->latentloglik = latentloglik;
  myestimator->jointlogp = jointlogp;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;

 return myestimator;

}
