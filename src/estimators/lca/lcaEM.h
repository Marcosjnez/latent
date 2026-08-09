/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 09/08/2026
 */

// Latent class analysis for Expectation-Maximization

class lcaEM: public estimators {

public:

  int S;
  int I;

  arma::vec weights;
  arma::vec logweights;

  arma::uvec indices_classes, indices_classloglik;

  arma::mat classes;
  arma::mat logclasses;
  arma::mat classloglik;
  arma::mat joint_classloglik;

  arma::vec loglik_case;
  arma::vec logliks;

  arma::mat posterior;
  arma::mat logposterior;

  // Posterior probabilities frozen during the E step:
  arma::mat posterior_em;

  double loss;
  double q_loss;
  double N;


  void param(arguments_optim& x) override {

    classes = arma::reshape(x.transparameters.elem(indices_classes), S, I);
    classloglik = arma::reshape(x.transparameters.elem(indices_classloglik), S, I);

    logclasses = arma::trunc_log(classes);

    for(int s = 0; s < S; ++s) {

      joint_classloglik.row(s) =
        classloglik.row(s) + logclasses.row(s);

      double max_vector = joint_classloglik.row(s).max();

      loglik_case(s) =
        max_vector +
        arma::trunc_log(
          arma::accu(
            arma::trunc_exp(
              joint_classloglik.row(s) - max_vector
            )
          )
        );

      logposterior.row(s) =
        joint_classloglik.row(s) - loglik_case(s);

      posterior.row(s) =
        arma::trunc_exp(logposterior.row(s));

      logliks(s) =
        weights(s) * loglik_case(s);

    }

  }


  // E step: freeze posterior probabilities
  void E(arguments_optim& x) override {

    posterior_em = posterior;
    x.posterior = posterior_em;

  }


  // Observed-data negative log-likelihood
  void observed_F(arguments_optim& x) override {

    loss = -arma::accu(logliks);
    x.f += loss;

  }


  // Gradient of the observed-data negative log-likelihood
  void observed_G(arguments_optim& x) override {

    arma::mat df_dclasses(S, I, arma::fill::zeros);
    arma::mat df_dclassloglik(S, I, arma::fill::zeros);

    for(int s = 0; s < S; ++s) {
      for(int i = 0; i < I; ++i) {

        df_dclasses(s, i) -=
          arma::trunc_exp(
            logweights(s) +
              classloglik(s, i) -
              loglik_case(s)
          );

        df_dclassloglik(s, i) -=
          arma::trunc_exp(
            logweights(s) +
              logposterior(s, i)
          );

      }
    }

    x.grad.elem(indices_classes) +=
      arma::vectorise(df_dclasses);

    x.grad.elem(indices_classloglik) +=
      arma::vectorise(df_dclassloglik);

  }


  // Expected complete-data negative log-likelihood:
  //
  // Q = - sum_s w_s sum_i tau_si *
  //       [log(class_si) + classloglik_si]
  //
  // posterior_em is kept fixed during the complete M step.
  void F(arguments_optim& x) override {

    q_loss = 0.0;

    for(int s = 0; s < S; ++s) {
      for(int i = 0; i < I; ++i) {

        q_loss -=
          weights(s) *
          posterior_em(s, i) *
          (logclasses(s, i) + classloglik(s, i));

      }
    }

    x.f += q_loss;

  }


  // Gradient of the Q-function
  void G(arguments_optim& x) override {

    arma::mat df_dclasses(S, I, arma::fill::zeros);
    arma::mat df_dclassloglik(S, I, arma::fill::zeros);

    for(int s = 0; s < S; ++s) {
      for(int i = 0; i < I; ++i) {

        double weighted_posterior =
          weights(s) * posterior_em(s, i);

        df_dclasses(s, i) -=
          weighted_posterior / classes(s, i);

        df_dclassloglik(s, i) -=
          weighted_posterior;

      }
    }

    x.grad.elem(indices_classes) +=
      arma::vectorise(df_dclasses);

    x.grad.elem(indices_classloglik) +=
      arma::vectorise(df_dclassloglik);

  }


  // Differential of the gradient of the Q-function
  void dG(arguments_optim& x) override {

    arma::mat dclasses =
      arma::reshape(
        x.dtransparameters.elem(indices_classes), S, I
      );

    arma::mat ddf_dclasses(S, I, arma::fill::zeros);
    arma::mat ddf_dclassloglik(S, I, arma::fill::zeros);

    for(int s = 0; s < S; ++s) {
      for(int i = 0; i < I; ++i) {

        double weighted_posterior =
          weights(s) * posterior_em(s, i);

        ddf_dclasses(s, i) +=
          weighted_posterior *
          dclasses(s, i) /
            (classes(s, i) * classes(s, i));

        // The Q-function is linear in classloglik, so:
        //
        // d/dtheta [dQ/dclassloglik] = 0
        //
        // The second derivatives with respect to the conditional
        // parameters are subsequently obtained by the likelihood
        // transformations in update_dgrad().

      }
    }

    x.dgrad.elem(indices_classes) +=
      arma::vectorise(ddf_dclasses);

    x.dgrad.elem(indices_classloglik) +=
      arma::vectorise(ddf_dclassloglik);

  }


  void outcomes(arguments_optim& x) override {

    doubles.resize(2);
    doubles[0] = loss;
    doubles[1] = -loss;

    names_doubles.resize(2);
    names_doubles[0] = "loss";
    names_doubles[1] = "loglik";

    vectors.resize(1);
    vectors[0] = loglik_case;

    names_vectors.resize(1);
    names_vectors[0] = "loglik_case";

    matrices.resize(1);
    matrices[0] = logposterior;

    names_matrices.resize(1);
    names_matrices[0] = "logposterior";

  }

};

lcaEM* choose_lcaEM(const Rcpp::List& estimator_setup) {

  lcaEM* myestimator = new lcaEM();

  std::vector<arma::uvec> indices = estimator_setup["indices"];

  int S = estimator_setup["S"];
  int I = estimator_setup["I"];

  arma::vec weights = estimator_setup["weights"];

  double N = arma::accu(weights);

  arma::mat classes(S, I, arma::fill::zeros);
  arma::mat classloglik(S, I, arma::fill::zeros);
  arma::mat logclasses(S, I, arma::fill::zeros);
  arma::mat joint_classloglik(S, I, arma::fill::zeros);

  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::mat posterior_em(S, I, arma::fill::zeros);

  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);

  arma::vec logweights = arma::trunc_log(weights);

  myestimator->S = S;
  myestimator->I = I;

  myestimator->weights = weights;
  myestimator->logweights = logweights;

  myestimator->indices_classes = indices[0];
  myestimator->indices_classloglik = indices[1];

  myestimator->classes = classes;
  myestimator->classloglik = classloglik;
  myestimator->logclasses = logclasses;
  myestimator->joint_classloglik = joint_classloglik;

  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;
  myestimator->posterior_em = posterior_em;

  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;

  myestimator->N = N;
  myestimator->loss = 0.0;
  myestimator->q_loss = 0.0;

  return myestimator;

}
