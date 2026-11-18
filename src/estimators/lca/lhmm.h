#ifndef LATENT_LHMM_H
#define LATENT_LHMM_H

#include <cmath>

/*
 * Hidden Markov model with trajectory-specific initial probabilities and a
 * shared transition matrix.
 * Follows the estimator interface in lca.h (latent commit e45add43c87b).
 *
 * indices: initial probabilities (Q x I), emission log-likelihoods (S x I),
 *          transition probabilities (I x I), in this order.
 * Q is the number of trajectories; initial probabilities can depend on covariates.
 * transition(i, j) = P(next state = i | previous state = j).
 * lengths: sizes of contiguous, independent, time-ordered trajectories.
 * weights: one positive weight per trajectory (not per observation).
 *
 * Softmax and emission transformations, including their derivatives, remain
 * outside this estimator. Internally, states are rows and times are columns
 * for contiguous access. No time-specific transformation objects are created.
 * Finite emission logs and interior softmax probabilities are the intended
 * inputs. A failed forward scale gives infinite loss; derivatives then error.
 */

class lhmm: public estimators {

public:

  int S, I, Q;
  arma::uvec lengths;
  arma::uvec indices_classes, indices_classloglik, indices_transition;
  arma::vec weights, scale, loglik_case;
  arma::mat classes, transition, emission, forward;
  double loss = 0.0;

  void param(arguments_optim& x) {

    classes = arma::reshape(x.transparameters.elem(indices_classes), Q, I);
    transition = arma::reshape(
      x.transparameters.elem(indices_transition), I, I);
    emission = arma::reshape(
      x.transparameters.elem(indices_classloglik), S, I).t();

    loss = 0.0;
    loglik_case.zeros();
    int start = 0;
    arma::vec predicted(I);

    for(arma::uword q = 0; q < lengths.n_elem; ++q) {

      const int end = start + lengths(q);

      for(int s = start; s < end; ++s) {

        if(s == start) predicted = classes.row(q).t();
        else predicted = transition * forward.col(s - 1);

        // Shift emission logs and normalize at every step to prevent underflow
        // from multiplying probabilities over a long trajectory.
        const double shift = emission.col(s).max();
        if(!std::isfinite(shift)) {
          loglik_case(q) = -arma::datum::inf;
          loss = arma::datum::inf;
          return;
        }
        emission.col(s) = arma::exp(emission.col(s) - shift);
        forward.col(s) = emission.col(s) % predicted;
        scale(s) = arma::accu(forward.col(s));

        if(!(scale(s) > 0.0) || !std::isfinite(scale(s))) {
          loglik_case(q) = -arma::datum::inf;
          loss = arma::datum::inf;
          return;
        }

        forward.col(s) /= scale(s);
        loglik_case(q) += shift + std::log(scale(s));

      }

      loss -= weights(q) * loglik_case(q);
      start = end;

    }

  }

  void F(arguments_optim& x) {

    x.f += loss;

  }

  void G(arguments_optim& x) {

    if(!std::isfinite(loss)) Rcpp::stop("lhmm: non-finite likelihood.");

    arma::mat df_dclasses(Q, I, arma::fill::zeros);
    arma::mat df_dtransition(I, I, arma::fill::zeros);
    arma::vec backward(I), u(I);
    int start = 0;

    for(arma::uword q = 0; q < lengths.n_elem; ++q) {

      const int end = start + lengths(q);
      const double w = weights(q);
      backward.ones();

      for(int s = end - 1; s >= start; --s) {

        u = emission.col(s) % backward / scale(s);

        for(int i = 0; i < I; ++i) {
          // Derivative w.r.t. emission log-likelihood = -w * smoothed probability.
          x.grad(indices_classloglik(s + S*i)) -=
            w * forward(i, s) * backward(i);
        }

        if(s == start) {
          df_dclasses.row(q) -= w * u.t();
        } else {
          df_dtransition -= w * u * forward.col(s - 1).t();
          backward = transition.t() * u;
        }

      }

      start = end;

    }

    x.grad.elem(indices_classes) += arma::vectorise(df_dclasses);
    x.grad.elem(indices_transition) += arma::vectorise(df_dtransition);

  }

  void dG(arguments_optim& x) {

    if(!std::isfinite(loss)) Rcpp::stop("lhmm: non-finite likelihood.");

    arma::mat dclasses = arma::reshape(x.dtransparameters.elem(indices_classes), Q, I);
    arma::mat dtransition = arma::reshape(
      x.dtransparameters.elem(indices_transition), I, I);
    arma::mat dclassloglik = arma::reshape(
      x.dtransparameters.elem(indices_classloglik), S, I).t();

    arma::mat dforward(I, S);
    arma::vec dscale(S);
    arma::mat ddf_dclasses(Q, I, arma::fill::zeros);
    arma::mat ddf_dtransition(I, I, arma::fill::zeros);
    arma::vec predicted(I), dpredicted(I), dvalue(I);
    arma::vec backward(I), dbackward(I), u(I), du(I);
    int start = 0;

    for(arma::uword q = 0; q < lengths.n_elem; ++q) {

      const int end = start + lengths(q);
      const double w = weights(q);

      // Differentiate the scaled forward recursion in the supplied direction.
      // Emission shifts are held constant: their differentials cancel.
      for(int s = start; s < end; ++s) {

        if(s == start) {
          predicted = classes.row(q).t();
          dpredicted = dclasses.row(q).t();
        } else {
          predicted = transition * forward.col(s - 1);
          dpredicted = dtransition * forward.col(s - 1) +
            transition * dforward.col(s - 1);
        }

        dvalue = emission.col(s) %
          (dpredicted + predicted % dclassloglik.col(s));
        dscale(s) = arma::accu(dvalue);
        dforward.col(s) =
          (dvalue - forward.col(s) * dscale(s)) / scale(s);

      }

      // Differentiate the backward recursion and accumulate d(gradient).
      backward.ones();
      dbackward.zeros();

      for(int s = end - 1; s >= start; --s) {

        u = emission.col(s) % backward / scale(s);
        du = (emission.col(s) %
          (dbackward + backward % dclassloglik.col(s)) -
          u * dscale(s)) / scale(s);

        for(int i = 0; i < I; ++i) {
          x.dgrad(indices_classloglik(s + S*i)) -= w *
            (dforward(i, s) * backward(i) + forward(i, s) * dbackward(i));
        }

        if(s == start) {
          ddf_dclasses.row(q) -= w * du.t();
        } else {
          ddf_dtransition -= w *
            (du * forward.col(s - 1).t() + u * dforward.col(s - 1).t());
          dbackward = dtransition.t() * u + transition.t() * du;
          backward = transition.t() * u;
        }

      }

      start = end;

    }

    x.dgrad.elem(indices_classes) += arma::vectorise(ddf_dclasses);
    x.dgrad.elem(indices_transition) += arma::vectorise(ddf_dtransition);

  }

  void outcomes(arguments_optim& x) {

    doubles = {loss, -loss};
    names_doubles = {"loss", "loglik"};
    vectors = {loglik_case};
    names_vectors = {"loglik_case"};
    matrices = {forward.t()};
    names_matrices = {"filtered"};

  }

};

inline lhmm* choose_lhmm(const Rcpp::List& estimator_setup) {

  std::vector<arma::uvec> indices = estimator_setup["indices"];
  const int S = estimator_setup["S"];
  const int I = estimator_setup["I"];

  arma::vec lengths(1);
  lengths(0) = S;
  if(estimator_setup.containsElementNamed("lengths")) {
    lengths = Rcpp::as<arma::vec>(estimator_setup["lengths"]);
  }
  if(lengths.is_empty() || !lengths.is_finite() ||
     arma::any(lengths < 1.0) || arma::any(lengths != arma::floor(lengths)) ||
     arma::accu(lengths) != S) {
    Rcpp::stop("lhmm: lengths must be positive integers summing to S.");
  }

  const int Q = lengths.n_elem;
  if(S < 1 || I < 1 || indices.size() != 3 ||
     indices[0].n_elem != arma::uword(Q)*I ||
     indices[1].n_elem != arma::uword(S)*I ||
     indices[2].n_elem != arma::uword(I)*I) {
    Rcpp::stop("lhmm: expected indices of sizes Q*I, S*I, and I*I.");
  }

  arma::vec weights(lengths.n_elem, arma::fill::ones);
  if(estimator_setup.containsElementNamed("weights")) {
    weights = Rcpp::as<arma::vec>(estimator_setup["weights"]);
  }
  if(weights.n_elem != lengths.n_elem || !weights.is_finite() ||
     arma::any(weights <= 0.0)) {
    Rcpp::stop("lhmm: supply one positive finite weight per trajectory.");
  }

  lhmm* myestimator = new lhmm();
  myestimator->S = S;
  myestimator->I = I;
  myestimator->Q = Q;
  myestimator->lengths = arma::conv_to<arma::uvec>::from(lengths);
  myestimator->weights = weights;
  myestimator->indices_classes = indices[0];
  myestimator->indices_classloglik = indices[1];
  myestimator->indices_transition = indices[2];
  myestimator->forward.set_size(I, S);
  myestimator->scale.set_size(S);
  myestimator->loglik_case.set_size(lengths.n_elem);

  return myestimator;

}

#endif
