#ifndef LATENT_LHMM_EM_H
#define LATENT_LHMM_EM_H

#include "lhmm.h"

// Hidden Markov model for Expectation-Maximization

class lhmmEM: public lhmm {

public:

  arma::uvec indices_probabilities;
  arma::vec weighted_states, weighted_counts;

  void param(arguments_optim& x) override {

    // Before the first E step, ordinary evaluations use the observed likelihood.
    // The numerical M step reads the trial parameters directly from x.
    if(weighted_states.is_empty()) lhmm::param(x);

  }

  void observed_F(arguments_optim& x) override {

    if(!weighted_states.is_empty()) lhmm::param(x);
    lhmm::F(x);

  }

  void observed_G(arguments_optim& x) override {

    lhmm::G(x);

  }

  void E(arguments_optim& x) override {

    if(weighted_states.is_empty()) {
      indices_probabilities = arma::join_cols(indices_classes, indices_transition);
      weighted_states.set_size(indices_classloglik.n_elem);
    }

    arma::mat initial_scores, transition_scores;
    smooth(initial_scores, transition_scores, [&](int k, double count) {
      weighted_states(k) = count;
    });

    // Expected initial and transition counts, including trajectory weights.
    weighted_counts = arma::join_cols(
      arma::vectorise(initial_scores % classes),
      arma::vectorise(transition_scores % transition));

  }

  // Expected complete-data negative log-likelihood with frozen counts
  void F(arguments_optim& x) override {

    if(weighted_states.is_empty()) return lhmm::F(x);

    double q_loss = 0.0;
    for(arma::uword k = 0; k < weighted_counts.n_elem; ++k) {
      if(weighted_counts(k) != 0.0) {
        q_loss -= weighted_counts(k) *
          std::log(x.transparameters(indices_probabilities(k)));
      }
    }
    for(arma::uword k = 0; k < weighted_states.n_elem; ++k) {
      if(weighted_states(k) != 0.0) {
        q_loss -= weighted_states(k) * x.transparameters(indices_classloglik(k));
      }
    }
    x.f += q_loss;

  }

  void G(arguments_optim& x) override {

    if(weighted_states.is_empty()) return lhmm::G(x);

    for(arma::uword k = 0; k < weighted_counts.n_elem; ++k) {
      if(weighted_counts(k) != 0.0) {
        const arma::uword j = indices_probabilities(k);
        x.grad(j) -= weighted_counts(k) / x.transparameters(j);
      }
    }
    x.grad.elem(indices_classloglik) -= weighted_states;

  }

  void dG(arguments_optim& x) override {

    if(weighted_states.is_empty()) return lhmm::dG(x);

    for(arma::uword k = 0; k < weighted_counts.n_elem; ++k) {
      if(weighted_counts(k) != 0.0) {
        const arma::uword j = indices_probabilities(k);
        const double p = x.transparameters(j);
        x.dgrad(j) += (weighted_counts(k) / p) * (x.dtransparameters(j) / p);
      }
    }
    // The emission log-likelihood enters the Q-function linearly.

  }

};

#endif
