/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 07/10/2025
 */

/*
 * Latent class analysis with expanded class probabilities
 */

class lca_cov: public estimators {

public:

  int S; // rows of Y (number of response patterns)
  int J; // cols of Y (number of items)
  int I; // Number of latent classes
  arma::uvec indices_classes;
  arma::uvec indices_items;
  arma::uvec indices_theta;
  arma::vec weights; // Number of repetitions of each response pattern
  arma::vec logweights;
  arma::mat classes;
  arma::mat logclasses;
  arma::vec logliks;
  arma::vec loglik_case;
  arma::mat jointlogp;
  arma::cube loglik, dloglik;
  arma::mat latentloglik;
  arma::mat posterior, logposterior;
  std::vector<arma::uvec> hess_indices;

  void param(arguments_optim& x) {

    arma::vec cl = x.transparameters.elem(indices_classes);
    classes = arma::reshape(cl, S, I);
    logclasses = arma::trunc_log(classes);

    arma::vec values = x.transparameters.elem(indices_items);
    std::memcpy(loglik.memptr(), values.memptr(), sizeof(double) * values.n_elem);

    latentloglik.zeros();

    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        for(int j=0; j < J; ++j) {
          latentloglik(s, i) += loglik(s, j, i);
        }
      }
      jointlogp.row(s) = latentloglik.row(s) + logclasses.row(s);
      double max_vector = jointlogp.row(s).max();
      loglik_case(s) = max_vector +
        arma::trunc_log(arma::accu(arma::trunc_exp(jointlogp.row(s) - max_vector)));
      logposterior.row(s) = jointlogp.row(s) - loglik_case(s);
      logliks(s) = weights(s) * loglik_case(s);
    }

  }

  void F(arguments_optim& x) {

    f = arma::accu(logliks);
    x.f -= f;

  }

  void G(arguments_optim& x) {

    // grad.set_size(transparameters.n_elem); grad.zeros();

    arma::mat dclasses(S, I, arma::fill::zeros);
    arma::mat ditemloglik(S, I, arma::fill::zeros);
    for(int s=0; s < S; ++s) {
      for(int i=0; i < I; ++i) {
        dclasses(s, i) -= arma::trunc_exp(logweights[s] + latentloglik(s, i) - loglik_case(s));
        ditemloglik(s, i) -= arma::trunc_exp(logweights[s] + logposterior(s, i));
      }
    }

    // Replicate ditemloglik across slices:
    dloglik.zeros();
    for (arma::uword k = 0; k < I; ++k) {
      dloglik.slice(k).each_col() = ditemloglik.col(k);
    }

    x.grad.elem(indices_classes) += arma::vectorise(dclasses);
    x.grad.elem(indices_items) += arma::vectorise(dloglik);
    // grad.elem( arma::find_nonfinite(grad) ).zeros();
    // x.grad.elem(indices[0]) += grad;

  }

  void dG(arguments_optim& x) {

    dg.set_size(transparameters.n_elem); dg.zeros();

  }

  void H(arguments_optim& x) {

    arma::mat posterior = arma::trunc_exp(logposterior);
    arma::mat M(J, J, arma::fill::ones);

    arma::mat post1 = posterior / classes;
    arma::mat post2 = post1; post2.each_col() %= weights;

    // FIX this
    x.hess(indices_classes, indices_classes) += post1.t() * post2;

    std::vector<arma::uvec> hess_indices2(hess_indices.size());
    for(int i=0; i < hess_indices.size(); ++i) {
      hess_indices2[i] = indices[0](hess_indices[i]);
    }

    int n = indices[0][0];
    for (int i = 0; i < I; ++i) {
      for (int k = i; k < I; ++k) {

        if (i == k) {

          arma::vec term = -weights % posterior.col(i) % (1.0 - posterior.col(k));

          arma::vec rep_vec = arma::repmat(term / classes(k), J, 1);
          x.hess.submat(hess_indices2[i], arma::uvec{ n+k }) += rep_vec;
          x.hess.submat(arma::uvec{ n+k }, hess_indices2[i]) += rep_vec.t();

          arma::mat block = arma::kron(M, arma::diagmat(term));
          x.hess.submat(hess_indices2[i], hess_indices2[k]) = block;

        } else {

          arma::vec term = weights % posterior.col(i) % posterior.col(k);

          arma::vec rep_vec = arma::repmat(term / classes(k), J, 1);
          x.hess.submat(hess_indices2[i], arma::uvec{ n+k }) += rep_vec;
          x.hess.submat(arma::uvec{ n+k }, hess_indices2[i]) += rep_vec.t();

          arma::vec rep_vec2 = arma::repmat(term / classes(i), J, 1);
          x.hess.submat(hess_indices2[k], arma::uvec{ n+i }) += rep_vec2;
          x.hess.submat(arma::uvec{ n+i }, hess_indices2[k]) += rep_vec2.t();

          arma::mat block = arma::kron(M, arma::diagmat(term));
          x.hess.submat(hess_indices2[i], hess_indices2[k]) += block;
          x.hess.submat(hess_indices2[k], hess_indices2[i]) += block;

        }
      }
    }

  }

  void M(arguments_optim& x) { // Update the parameter estimates

  }

  void E(arguments_optim& x) { // Update the loglik and posterior

    // Estimated posterior:
    x.posterior = arma::trunc_exp(logposterior);
    // New hypothetical frequencies:
    x.freqs = x.posterior;
    x.freqs.each_col() %= weights;
    x.loglik = f;

    // Update the class probabilities:
    // New number of subjects in each class:
    arma::vec freqs_i = arma::sum(x.freqs, 0).t();
    // New proportion of subjects in each class:
    classes = freqs_i / arma::accu(freqs_i);
    logclasses = arma::trunc_log(classes);
    // Put a zero in the first element:
    logclasses -= logclasses(0);

    x.transparameters(indices_classes) = classes;
    x.transparameters(indices_theta) = logclasses;

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(1);
    doubles[0] = f;

    vectors.resize(3);
    vectors[0] = loglik_case;
    vectors[1] = logliks;
    vectors[2] = weights;

    matrices.resize(5);
    matrices[0] = latentloglik;
    matrices[1] = logposterior;
    matrices[2] = jointlogp;
    matrices[3] = classes;
    matrices[4] = logclasses;

    cubes.resize(2);
    cubes[0] = loglik;
    cubes[1] = dloglik;

  }

};

lca_cov* choose_lca_cov(const Rcpp::List& estimator_setup) {

  lca_cov* myestimator = new lca_cov();

  int S = estimator_setup["S"];
  int J = estimator_setup["J"];
  int I = estimator_setup["I"];
  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::vec weights = estimator_setup["weights"];
  std::vector<arma::uvec> hess_indices = estimator_setup["hess_indices"];

  arma::mat classes(S, I);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::cube loglik(S, J, I, arma::fill::zeros);
  arma::cube dloglik(S, J, I, arma::fill::zeros);
  arma::mat latentloglik(S, I, arma::fill::zeros);
  arma::mat jointlogp(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);
  arma::uvec indices_classes = indices[1];
  arma::uvec indices_items = indices[2];
  arma::uvec indices_theta = indices[3];

  myestimator->S = S;
  myestimator->J = J;
  myestimator->I = I;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices = indices;
  myestimator->indices_classes = indices_classes;
  myestimator->indices_items = indices_items;
  myestimator->indices_theta = indices_theta;
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
