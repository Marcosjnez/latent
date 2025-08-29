/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 29/08/2025
 */

/*
 * Latent class analysis
 */

class lca2: public estimators {

public:

  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  int I; // Number of latent classes
  arma::uvec indices_classes;
  arma::uvec indices_items;
  arma::vec weights; // Number of repetitions of each response pattern
  arma::vec logweights;
  arma::vec classes;
  arma::vec logclasses;
  arma::vec logliks;
  arma::vec loglik_case;
  arma::mat jointlogp;
  arma::cube loglik, dloglik;
  arma::mat latentloglik;
  arma::mat posterior, logposterior;
  std::vector<arma::uvec> hess_indices;

  void param() {

    indices_classes = indices[1];
    indices_items = indices[2];
    classes = transparameters(indices_classes);
    logclasses = arma::trunc_log(classes);

    // loglik(transparameters(indices_items).memptr(), S, J, I);
    arma::vec values = transparameters(indices_items);
    std::memcpy(loglik.memptr(), values.memptr(), sizeof(double) * values.n_elem);

    latentloglik.zeros();

    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        for(int j=0; j < J; ++j) {
          latentloglik(s, i) += loglik(s, j, i);
        }
      }
      jointlogp.row(s) = latentloglik.row(s) + logclasses.t();
      double max_vector = jointlogp.row(s).max();
      loglik_case(s) = max_vector +
        arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      logposterior.row(s) = jointlogp.row(s) - loglik_case(s);
      logliks(s) = weights(s) * loglik_case(s);
    }

  }

  void F() {

    f = -arma::accu(logliks);

  }

  void G() {

    grad.set_size(transparameters.n_elem); grad.zeros();

    arma::vec dclasses(I, arma::fill::zeros);
    arma::mat ditemloglik(S, I, arma::fill::zeros);
    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        dclasses(i) -= arma::trunc_exp(logweights[s] + latentloglik(s, i) - loglik_case(s));
        ditemloglik(s, i) -= arma::trunc_exp(logweights[s] + logposterior(s, i));
      }
    }

    dloglik.zeros();
    for (arma::uword k = 0; k < I; ++k) {
      dloglik.slice(k).each_col() = ditemloglik.col(k);
    }

    grad.elem(indices_classes) += dclasses;
    grad.elem(indices_items) += arma::vectorise(dloglik);
    // grad.elem( arma::find_nonfinite(grad) ).zeros();
  }

  void dG() {

    dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H() {

    hess.set_size(transparameters.n_elem, transparameters.n_elem);
    hess.zeros();

    arma::mat posterior = arma::trunc_exp(logposterior);
    arma::mat M(J, J, arma::fill::ones);

    arma::mat post1 = posterior; post1.each_row() /= classes.t();
    arma::mat post2 = post1; post2.each_col() %= weights;
    hess(indices_classes, indices_classes) = post1.t() * post2;

    for (int i = 0; i < I; ++i) {
      for (int k = i; k < I; ++k) {
        arma::vec term;

        if (i == k) {

          term = -weights % posterior.col(i) % (1.0 - posterior.col(k));

          arma::vec rep_vec = arma::repmat(term / classes(k), J, 1);
          hess.submat(hess_indices[i], arma::uvec{ static_cast<arma::uword>(k) }) = rep_vec;
          hess.submat(arma::uvec{ static_cast<arma::uword>(k) }, hess_indices[i]) = rep_vec.t();

          arma::mat block = arma::kron(M, arma::diagmat(term));
          hess.submat(hess_indices[i], hess_indices[k]) = block;

        } else {

          term = weights % posterior.col(i) % posterior.col(k);

          arma::vec rep_vec = arma::repmat(term / classes(k), J, 1);
          hess.submat(hess_indices[i], arma::uvec{ static_cast<arma::uword>(k) }) = rep_vec;
          hess.submat(arma::uvec{ static_cast<arma::uword>(k) }, hess_indices[i]) = rep_vec.t();

          arma::vec rep_vec2 = arma::repmat(term / classes(i), J, 1);
          hess.submat(hess_indices[k], arma::uvec{ static_cast<arma::uword>(i) }) = rep_vec2;
          hess.submat(arma::uvec{ static_cast<arma::uword>(i) }, hess_indices[k]) = rep_vec2.t();

          arma::mat block = arma::kron(M, arma::diagmat(term));
          hess.submat(hess_indices[i], hess_indices[k]) = block;
          hess.submat(hess_indices[k], hess_indices[i]) = block;

        }
      }
    }

  }

  void E() { // Update the parameter estimates

  }

  void M() { // Update the posterior probabilities

  }

  void outcomes() {

    doubles.resize(1);
    doubles[0] = f;

    vectors.resize(5);
    vectors[0] = classes;
    vectors[1] = logclasses;
    vectors[2] = loglik_case;
    vectors[3] = logliks;
    vectors[4] = weights;

    matrices.resize(3);
    matrices[0] = latentloglik;
    matrices[1] = logposterior;
    matrices[2] = jointlogp;

    cubes.resize(2);
    cubes[0] = loglik;
    cubes[1] = dloglik;

  }

};

lca2* choose_lca2(const Rcpp::List& estimator_setup) {

  lca2* myestimator = new lca2();

  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  int I = estimator_setup["I"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::vec weights = estimator_setup["weights"];
  std::vector<arma::uvec> hess_indices = estimator_setup["hess_indices"];

  arma::vec classes(I);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::cube loglik(S, J, I, arma::fill::zeros);
  arma::cube dloglik(S, J, I, arma::fill::zeros);
  arma::mat latentloglik(S, I, arma::fill::zeros);
  arma::mat jointlogp(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);

  myestimator->S = S;
  myestimator->J = J;
  myestimator->I = I;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices = indices;
  myestimator->hess_indices = hess_indices;

  myestimator->classes = classes;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->loglik = loglik;
  myestimator->dloglik = dloglik;
  myestimator->latentloglik = latentloglik;
  myestimator->jointlogp = jointlogp;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;

  return myestimator;

}
