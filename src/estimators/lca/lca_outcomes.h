/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 09/06/2026
 */

/*
 * Latent class analysis with expanded class probabilities and distal outcomes
 */

/*
 * Latent class analysis with expanded class probabilities and distal outcomes
 */

class lca_outcomes: public estimators {

public:

  int S;
  int I;
  int P;
  arma::vec weights;
  arma::vec logweights;
  arma::uvec indices_classes, indices_classloglik, indices_betas;
  arma::mat Y, res;
  arma::mat classes;
  arma::mat logclasses;
  arma::mat classloglik;
  arma::mat joint_classloglik;
  arma::mat beta;
  arma::vec loglik_case;
  arma::vec logliks;
  arma::mat posterior, logposterior;
  double loss, N, loss_lca, loss_outcomes;

  void param(arguments_optim& x) {

    classes = arma::reshape(x.transparameters.elem(indices_classes), S, I);
    classloglik = arma::reshape(x.transparameters.elem(indices_classloglik), S, I);
    beta = arma::reshape(x.transparameters.elem(indices_betas), I, P);

    logclasses = arma::trunc_log(classes);

    for (int s = 0; s < S; ++s) {

      joint_classloglik.row(s) = classloglik.row(s) + logclasses.row(s);

      double max_vector = joint_classloglik.row(s).max();

      loglik_case(s) = max_vector +
        arma::trunc_log(
          arma::accu(
            arma::trunc_exp(joint_classloglik.row(s) - max_vector)
          )
        );

      logposterior.row(s) = joint_classloglik.row(s) - loglik_case(s);
      posterior.row(s) = arma::trunc_exp(logposterior.row(s));
      logliks(s) = weights(s) * loglik_case(s);

    }

    res = Y - logposterior * beta;

  }

  void F(arguments_optim& x) {

    loss_lca = -arma::accu(logliks);
    loss_outcomes = 0.5 * arma::accu(res % res);
    loss = loss_lca + loss_outcomes;

    x.f += loss;

  }

  void G(arguments_optim& x) {

    arma::mat df_dclasses(S, I, arma::fill::zeros);
    arma::mat df_dclassloglik(S, I, arma::fill::zeros);

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {

        df_dclasses(s, i) -=
          arma::trunc_exp(logweights(s) + classloglik(s, i) - loglik_case(s));

        df_dclassloglik(s, i) -=
          arma::trunc_exp(logweights(s) + logposterior(s, i));

      }
    }

    arma::mat df_dlogposterior = -res * beta.t();

    for (int s = 0; s < S; ++s) {

      double acc = arma::accu(df_dlogposterior.row(s));

      for (int i = 0; i < I; ++i) {

        double q_si = df_dlogposterior(s, i) - posterior(s, i) * acc;

        df_dclasses(s, i) += q_si / classes(s, i);
        df_dclassloglik(s, i) += q_si;

      }
    }

    arma::mat df_dbeta = -logposterior.t() * res;

    x.grad.elem(indices_classes) += arma::vectorise(df_dclasses);
    x.grad.elem(indices_classloglik) += arma::vectorise(df_dclassloglik);
    x.grad.elem(indices_betas) += arma::vectorise(df_dbeta);

  }

  void dG(arguments_optim& x) {

    arma::mat dclasses = arma::reshape(x.dtransparameters.elem(indices_classes), S, I);
    arma::mat dclassloglik = arma::reshape(x.dtransparameters.elem(indices_classloglik), S, I);
    arma::mat dbeta = arma::reshape(x.dtransparameters.elem(indices_betas), I, P);

    arma::mat ddf_dclasses(S, I, arma::fill::zeros);
    arma::mat ddf_dclassloglik(S, I, arma::fill::zeros);

    arma::mat dlogclasses(S, I, arma::fill::zeros);
    arma::mat djoint_classloglik(S, I, arma::fill::zeros);
    arma::vec dloglik_case(S, arma::fill::zeros);
    arma::mat dlogposterior(S, I, arma::fill::zeros);
    arma::mat dposterior(S, I, arma::fill::zeros);

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        dlogclasses(s, i) = dclasses(s, i) / classes(s, i);
        djoint_classloglik(s, i) = dclassloglik(s, i) + dlogclasses(s, i);
      }
    }

    for (int s = 0; s < S; ++s) {
      double acc = 0.0;
      for (int i = 0; i < I; ++i) {
        acc += posterior(s, i) * djoint_classloglik(s, i);
      }
      dloglik_case(s) = acc;
    }

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {
        dlogposterior(s, i) = djoint_classloglik(s, i) - dloglik_case(s);
        dposterior(s, i) = posterior(s, i) * dlogposterior(s, i);
      }
    }

    for (int s = 0; s < S; ++s) {
      for (int i = 0; i < I; ++i) {

        double term_classes =
          arma::trunc_exp(logweights(s) + classloglik(s, i) - loglik_case(s));

        double term_classloglik =
          arma::trunc_exp(logweights(s) + logposterior(s, i));

        ddf_dclasses(s, i) -=
          term_classes * (dclassloglik(s, i) - dloglik_case(s));

        ddf_dclassloglik(s, i) -=
          term_classloglik * dlogposterior(s, i);

      }
    }

    arma::mat dres = -(dlogposterior * beta + logposterior * dbeta);

    arma::mat df_dlogposterior = -res * beta.t();
    arma::mat ddf_dlogposterior =
      -(dres * beta.t() + res * dbeta.t());

      for (int s = 0; s < S; ++s) {

        double acc = arma::accu(df_dlogposterior.row(s));
        double dacc = arma::accu(ddf_dlogposterior.row(s));

        for (int i = 0; i < I; ++i) {

          double q_si =
            df_dlogposterior(s, i) - posterior(s, i) * acc;

          double dq_si =
            ddf_dlogposterior(s, i) -
            dposterior(s, i) * acc -
            posterior(s, i) * dacc;

          ddf_dclasses(s, i) +=
            dq_si / classes(s, i) -
            q_si * dclasses(s, i) / (classes(s, i) * classes(s, i));

          ddf_dclassloglik(s, i) += dq_si;

        }
      }

      arma::mat ddf_dbeta =
        -dlogposterior.t() * res -
        logposterior.t() * dres;

      x.dgrad.elem(indices_classes) += arma::vectorise(ddf_dclasses);
      x.dgrad.elem(indices_classloglik) += arma::vectorise(ddf_dclassloglik);
      x.dgrad.elem(indices_betas) += arma::vectorise(ddf_dbeta);

  }

  void outcomes(arguments_optim& x) {

    doubles.resize(7);
    doubles[0] =  loss;
    doubles[1] =  0.00;
    doubles[2] =  0.00;
    doubles[3] = -loss_lca;
    doubles[4] =  0.00;
    doubles[5] =  0.00;
    doubles[6] =  loss_outcomes;

    vectors.resize(1);
    vectors[0] = loglik_case;

    matrices.resize(3);
    matrices[0] = logposterior;
    matrices[1] = beta;
    matrices[2] = res;

  }

};

lca_outcomes* choose_lca_outcomes(const Rcpp::List& estimator_setup) {

  lca_outcomes* myestimator = new lca_outcomes();

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  arma::mat Y = estimator_setup["Y"];
  int S = estimator_setup["S"];
  int I = estimator_setup["I"];
  arma::vec weights = estimator_setup["weights"];

  int P = Y.n_cols;
  double N = arma::accu(weights);

  arma::mat classes(S, I, arma::fill::zeros);
  arma::mat classloglik(S, I, arma::fill::zeros);
  arma::mat logclasses(S, I, arma::fill::zeros);
  arma::mat joint_classloglik(S, I, arma::fill::zeros);
  arma::mat posterior(S, I, arma::fill::zeros);
  arma::mat logposterior(S, I, arma::fill::zeros);
  arma::mat beta(I, P, arma::fill::zeros);
  arma::mat res(S, P, arma::fill::zeros);
  arma::vec logliks(S, arma::fill::zeros);
  arma::vec loglik_case(S, arma::fill::zeros);
  arma::vec logweights = arma::trunc_log(weights);

  myestimator->S = S;
  myestimator->I = I;
  myestimator->P = P;
  myestimator->weights = weights;
  myestimator->logweights = logweights;
  myestimator->indices_classes = indices[0];
  myestimator->indices_classloglik = indices[1];
  myestimator->indices_betas = indices[2];
  myestimator->Y = Y;
  myestimator->classes = classes;
  myestimator->classloglik = classloglik;
  myestimator->logclasses = logclasses;
  myestimator->joint_classloglik = joint_classloglik;
  myestimator->posterior = posterior;
  myestimator->logposterior = logposterior;
  myestimator->beta = beta;
  myestimator->res = res;
  myestimator->logliks = logliks;
  myestimator->loglik_case = loglik_case;
  myestimator->N = N;
  myestimator->loss = 0.0;
  myestimator->loss_lca = 0.0;
  myestimator->loss_outcomes = 0.0;

  return myestimator;

}
