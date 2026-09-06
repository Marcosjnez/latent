/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 06/09/2026
 *
 */

#include <algorithm>
#include <map>
#include "poly_acov_parallel.h"

arma::vec diagcov(arma::mat X) {

  const size_t numCols = X.n_cols;
  arma::vec diag_cov(numCols);

  for(size_t i = 0; i < numCols; ++i) {
    arma::vec col = X.col(i);
    arma::uvec validIndices = arma::find_finite(col);
    arma::vec validCol = col(validIndices);
    validCol -= arma::mean(validCol);
    double d = arma::accu(validCol % validCol)/(validIndices.n_elem-1);
    diag_cov(i) = d;
  }

  return diag_cov;

}

arma::mat center_finite(arma::mat X) {

  // Center each column ignoring NaNs

  for(arma::uword j = 0; j < X.n_cols; ++j) {
    arma::vec col = X.col(j);
    arma::uvec idx = arma::find_finite(col);
    if(idx.n_elem > 0) {
      double mu = arma::mean(col(idx));
      col(idx) -= mu;
      X.col(j) = col;
    }
  }

  return X;

}

arma::mat asymptotic_general(arma::mat X, bool cov, bool diag) {

  /*
   * Browne and Shapiro (Equation 3.2; 1986)
   *
   * cov = false : asymptotic covariance matrix of item correlations
   *               (off-diagonal only)
   * cov = true  : asymptotic covariance matrix of item covariances
   *               (including variances)
   */

  int q = X.n_cols;
  int qq = q*q;

  // Center variables first (important for fourth moments)
  arma::mat Xc = center_finite(X);

  arma::vec d;
  arma::mat S, P;

  if(Xc.has_nan()) {
    d = arma::sqrt(diagcov(Xc));
    P = pairwise_cor(Xc);
    arma::mat D = arma::diagmat(d);
    S = D*P*D;
  } else {
    S = (Xc.t()*Xc)/Xc.n_rows;
    d = arma::sqrt(arma::diagvec(S));
    arma::mat Dinv = arma::diagmat(1.0/d);
    P = Dinv*S*Dinv;
  }

  arma::mat Theta(qq, qq, arma::fill::zeros);

  int ij = 0;
  int kh = 0;

  for(int j = 0; j < q; ++j) {
    for(int i = 0; i < q; ++i) {

      kh = 0;

      for(int h = 0; h < q; ++h) {
        for(int k = 0; k < q; ++k) {

          arma::vec m = Xc.col(i) % Xc.col(j) % Xc.col(k) % Xc.col(h);
          arma::uvec validIndices = arma::find_finite(m);
          arma::vec v = m(validIndices);

          double val = arma::mean(v);

          if(cov) {
            Theta(ij, kh) = val;
          } else {
            Theta(ij, kh) = val/(d[i]*d[j]*d[k]*d[h]);
          }

          ++kh;
        }
      }

      ++ij;
    }
  }

  if(cov) {
    arma::vec s = arma::vectorise(S);
    arma::mat Gamma = Theta-s*s.t();

    arma::uvec lower_indices = arma::trimatl_ind(arma::size(S), 0);
    return Gamma(lower_indices, lower_indices);
  }

  arma::vec p = arma::vectorise(P);
  arma::mat Gamma = Theta-p*p.t();

  arma::mat Ms = dxt(q, q)*0.5;
  Ms.diag() += 0.5;

  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);

  for(int i = 0; i < q; ++i) {
    int ii = i*q+i;
    Kd(ii, i) = 1.0;
  }

  arma::mat A = Ms*arma::kron(I, P)*Kd;
  arma::mat B = Gamma*Kd;
  arma::mat G = Kd.t()*Gamma*Kd;

  arma::mat asymptotic = Gamma-A*B.t()-B*A.t()+A*G*A.t();

  arma::uvec lower_indices;
  if(diag) {
    lower_indices = arma::trimatl_ind(arma::size(P), 0);
  } else {
    lower_indices = arma::trimatl_ind(arma::size(P), -1);
  }

  return asymptotic(lower_indices, lower_indices);

}

arma::mat asymptotic_normal(const arma::mat& S, bool cov, bool diag) {

  /*
   * Browne and Shapiro (Equation 4.1; 1986)
   *
   * Input is the covariance matrix S.
   *
   * cov = false : asymptotic covariance matrix of item correlations
   * cov = true  : asymptotic covariance matrix of item covariances
   */

  int q = S.n_rows;
  int qq = q*q;

  arma::mat Ms = dxt(q, q)*0.5;
  Ms.diag() += 0.5;

  if(cov) {
    arma::mat Gamma = 2.0*Ms*arma::kron(S, S);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(S), 0);
    return Gamma(lower_indices, lower_indices);
  }

  arma::vec d = arma::sqrt(arma::diagvec(S));
  arma::mat Dinv = arma::diagmat(1.0/d);
  arma::mat P = Dinv*S*Dinv;

  arma::mat I(q, q, arma::fill::eye);
  arma::mat Kd(qq, q, arma::fill::zeros);

  for(int i = 0; i < q; ++i) {
    int ii = i*q+i;
    Kd(ii, i) = 1.0;
  }

  arma::mat A = Ms*arma::kron(I, P)*Kd;
  arma::mat Gamma = 2.0*Ms*arma::kron(P, P);
  arma::mat B = Gamma*Kd;
  arma::mat G = 2.0*P % P;

  arma::mat asymptotic = Gamma-A*B.t()-B*A.t()+A*G*A.t();

  arma::uvec lower_indices;
  if(diag) {
    lower_indices = arma::trimatl_ind(arma::size(P), 0);
  } else {
    lower_indices = arma::trimatl_ind(arma::size(P), -1);
  }

  return asymptotic(lower_indices, lower_indices);

}

arma::mat asymptotic_elliptical(const arma::mat& S, double eta, bool cov,
                                bool diag) {

  /*
   * Browne and Shapiro style elliptical correction.
   *
   * Input is the covariance matrix S.
   * Here eta = 1+kappa.
   */

  if(!cov) {
    return eta*asymptotic_normal(S, false, false);
  }

  int q = S.n_rows;
  int qq = q*q;

  arma::mat Ms = dxt(q, q)*0.5;
  Ms.diag() += 0.5;

  arma::mat Gamma_normal = 2.0*Ms*arma::kron(S, S);
  arma::vec s = arma::vectorise(S);

  arma::mat Gamma = eta*Gamma_normal+(eta-1.0)*(s*s.t());

  arma::uvec lower_indices;
  if(diag) {
    lower_indices = arma::trimatl_ind(arma::size(S), 0);
  } else {
    lower_indices = arma::trimatl_ind(arma::size(S), -1);
  }

  return Gamma(lower_indices, lower_indices);

}


namespace latent_asymptotic_poly {

inline double normal_cdf_bound(const double x) {

  if(x == neg_inf) return 0.0;
  if(x == pos_inf) return 1.0;

  return Pnorm(x);

}

inline double bvn_pdf_dx(const double rho,
                         const double x,
                         const double y) {

  if(!std::isfinite(x) || !std::isfinite(y)) return 0.0;

  const double one_minus_rho2 = 1.0-rho*rho;

  return -((x-rho*y)/one_minus_rho2)*dbinorm(rho, x, y);

}

inline double bvn_pdf_dy(const double rho,
                         const double x,
                         const double y) {

  if(!std::isfinite(x) || !std::isfinite(y)) return 0.0;

  const double one_minus_rho2 = 1.0-rho*rho;

  return -((y-rho*x)/one_minus_rho2)*dbinorm(rho, x, y);

}

inline void first_threshold_derivatives(const double rho,
                                        const double threshold,
                                        const double sign,
                                        const double lower_other,
                                        const double upper_other,
                                        double& dprobability,
                                        double& drho_threshold) {

  const double denominator = std::sqrt(1.0-rho*rho);
  const double upper = (upper_other-rho*threshold)/denominator;
  const double lower = (lower_other-rho*threshold)/denominator;

  dprobability = sign*Dnorm(threshold)*(Pnorm(upper)-Pnorm(lower));
  drho_threshold = sign*(
    bvn_pdf_dx(rho, threshold, upper_other)-
    bvn_pdf_dx(rho, threshold, lower_other)
  );

}

inline void second_threshold_derivatives(const double rho,
                                         const double threshold,
                                         const double sign,
                                         const double lower_other,
                                         const double upper_other,
                                         double& dprobability,
                                         double& drho_threshold) {

  const double denominator = std::sqrt(1.0-rho*rho);
  const double upper = (upper_other-rho*threshold)/denominator;
  const double lower = (lower_other-rho*threshold)/denominator;

  dprobability = sign*Dnorm(threshold)*(Pnorm(upper)-Pnorm(lower));
  drho_threshold = sign*(
    bvn_pdf_dy(rho, upper_other, threshold)-
    bvn_pdf_dy(rho, lower_other, threshold)
  );

}

inline arma::vec finite_thresholds(const Rcpp::RObject& object,
                                   const arma::uword variable) {

  Rcpp::NumericVector input(object);
  arma::vec result(input.size());

  for(R_xlen_t i = 0; i < input.size(); ++i) {
    result[static_cast<arma::uword>(i)] = input[i];
  }

  arma::uword first = 0L;
  arma::uword last = result.n_elem;

  if(last > 0L && result[0L] == neg_inf) ++first;
  if(last > first && result[last-1L] == pos_inf) --last;

  if(last <= first) {
    Rcpp::stop("Variable " + std::to_string(variable+1L) +
      " does not contain any finite threshold.");
  }

  result = result.subvec(first, last-1L);

  if(!result.is_finite()) {
    Rcpp::stop("All internal thresholds must be finite.");
  }

  if(result.n_elem > 1L &&
     arma::any(arma::diff(result) <= 0.0)) {
    Rcpp::stop("Thresholds must be strictly increasing within variables.");
  }

  return result;

}

inline arma::uword category_index(const arma::vec& levels,
                                  const double value) {

  const double* first = levels.begin();
  const double* last = levels.end();
  const double* position = std::lower_bound(first, last, value);

  if(position == last || *position != value) {
    Rcpp::stop("An observed category could not be matched to its levels.");
  }

  return static_cast<arma::uword>(position-first);

}

} // namespace latent_asymptotic_poly

arma::mat composite_poly_scores(const arma::mat& data,
                                const Rcpp::List& thresholds,
                                const arma::mat& correlation) {

  using namespace latent_asymptotic_poly;

  const arma::uword nobs = data.n_rows;
  const arma::uword nitems = data.n_cols;
  const double probability_floor = 1e-16;

  if(nobs < 1L || nitems < 2L) {
    Rcpp::stop("data must contain at least one observation and two variables.");
  }

  if(!data.is_finite()) {
    Rcpp::stop("composite_poly_scores() currently requires complete ordinal data.");
  }

  if(correlation.n_rows != nitems ||
     correlation.n_cols != nitems ||
     !correlation.is_finite()) {
    Rcpp::stop("correlation must be a finite square matrix matching data.");
  }

  if(!arma::approx_equal(correlation, correlation.t(),
                         "absdiff", 1e-08)) {
    Rcpp::stop("correlation must be symmetric.");
  }

  if(arma::any(arma::abs(correlation.diag()-1.0) > 1e-08)) {
    Rcpp::stop("correlation must have a unit diagonal.");
  }

  if(thresholds.size() != static_cast<R_xlen_t>(nitems)) {
    Rcpp::stop("thresholds must contain one vector per observed variable.");
  }

  std::vector<arma::vec> tau(nitems);
  std::vector<arma::vec> levels(nitems);
  std::vector<arma::vec> bounds(nitems);
  std::vector<arma::vec> bounds_cdf(nitems);

  arma::uvec threshold_offsets(nitems, arma::fill::zeros);
  arma::umat categories(nobs, nitems, arma::fill::zeros);

  arma::uword nthresholds = 0L;

  for(arma::uword j = 0L; j < nitems; ++j) {

    tau[j] = finite_thresholds(thresholds[j], j);
    threshold_offsets[j] = nthresholds;
    nthresholds += tau[j].n_elem;

    levels[j] = arma::sort(arma::unique(data.col(j)));

    if(levels[j].n_elem != tau[j].n_elem+1L) {
      Rcpp::stop("The number of thresholds for variable " +
        std::to_string(j+1L) +
        " does not match its observed categories. The threshold representation "
        "cannot contain empty internal categories.");
    }

    for(arma::uword i = 0L; i < nobs; ++i) {
      categories(i, j) = category_index(levels[j], data(i, j));
    }

    bounds[j].set_size(tau[j].n_elem+2L);
    bounds[j][0L] = neg_inf;
    bounds[j].subvec(1L, tau[j].n_elem) = tau[j];
    bounds[j][bounds[j].n_elem-1L] = pos_inf;

    bounds_cdf[j].set_size(bounds[j].n_elem);

    for(arma::uword h = 0L; h < bounds[j].n_elem; ++h) {
      bounds_cdf[j][h] = normal_cdf_bound(bounds[j][h]);
    }

  }

  const arma::uword ncorrelations = nitems*(nitems-1L)/2L;
  const arma::uword nparameters = nthresholds+ncorrelations;
  arma::mat scores(nobs, nparameters, arma::fill::zeros);
  arma::uword floored_probabilities = 0L;

  arma::uword q = 0L;

  for(arma::uword j = 0L; j < nitems-1L; ++j) {

    for(arma::uword k = j+1L; k < nitems; ++k, ++q) {

      const double rho = correlation(j, k);

      if(std::abs(rho) >= 1.0-1e-10) {
        Rcpp::stop("Every polychoric correlation must lie strictly inside (-1, 1).");
      }

      for(arma::uword i = 0L; i < nobs; ++i) {

        const arma::uword category_j = categories(i, j);
        const arma::uword category_k = categories(i, k);

        const double lower_j = bounds[j][category_j];
        const double upper_j = bounds[j][category_j+1L];
        const double lower_k = bounds[k][category_k];
        const double upper_k = bounds[k][category_k+1L];

        double probability = pbinorm(
          rho,
          lower_j, lower_k,
          upper_j, upper_k,
          bounds_cdf[j][category_j],
          bounds_cdf[k][category_k],
          bounds_cdf[j][category_j+1L],
          bounds_cdf[k][category_k+1L]
        );

        if(!std::isfinite(probability) || probability < probability_floor) {
          probability = probability_floor;
          ++floored_probabilities;
        }

        const double probability_rho =
          dbinorm(rho, upper_j, upper_k)-
          dbinorm(rho, lower_j, upper_k)-
          dbinorm(rho, upper_j, lower_k)+
          dbinorm(rho, lower_j, lower_k);

        scores(i, nthresholds+q) += probability_rho/probability;

        auto add_first_score = [&](const arma::uword threshold,
                                   const double sign) {

          double probability_threshold = 0.0;
          double rho_threshold = 0.0;

          first_threshold_derivatives(
            rho, tau[j][threshold], sign,
            lower_k, upper_k,
            probability_threshold, rho_threshold
          );

          scores(i, threshold_offsets[j]+threshold) +=
            probability_threshold/probability;

        };

        auto add_second_score = [&](const arma::uword threshold,
                                    const double sign) {

          double probability_threshold = 0.0;
          double rho_threshold = 0.0;

          second_threshold_derivatives(
            rho, tau[k][threshold], sign,
            lower_j, upper_j,
            probability_threshold, rho_threshold
          );

          scores(i, threshold_offsets[k]+threshold) +=
            probability_threshold/probability;

        };

        if(category_j > 0L) {
          add_first_score(category_j-1L, -1.0);
        }

        if(category_j < tau[j].n_elem) {
          add_first_score(category_j, 1.0);
        }

        if(category_k > 0L) {
          add_second_score(category_k-1L, -1.0);
        }

        if(category_k < tau[k].n_elem) {
          add_second_score(category_k, 1.0);
        }

      }

    }

  }

  if(floored_probabilities > 0L) {
    Rcpp::warning(std::to_string(floored_probabilities) +
      " bivariate probabilities were replaced by 1e-16.");
  }

  return scores;

}

Rcpp::List asymptotic_poly(const arma::mat& data,
                           const arma::mat& correlation,
                           const Rcpp::List& thresholds,
                           bool return_scores,
                           double probability_floor,
                           double inversion_tolerance,
                           Rcpp::Nullable<Rcpp::List> polyfast_object,
                           const int cores) {

  using namespace latent_asymptotic_poly;

  // Check inputs

  if(cores < 1) {
    Rcpp::stop("cores must be a positive integer.");
  }

#ifndef _OPENMP
  if(cores > 1) {
    Rcpp::warning("OpenMP is not available in this build; ACOV OpenMP regions will run serially.");
  }
#endif

  const arma::uword nobs = data.n_rows;
  const arma::uword nitems = data.n_cols;

  if(nobs < 2L || nitems < 2L) {
    Rcpp::stop("data must contain at least two observations and two variables.");
  }

  if(!data.is_finite()) {
    Rcpp::stop("asymptotic_poly() currently requires complete ordinal data.");
  }

  if(correlation.n_rows != nitems ||
     correlation.n_cols != nitems ||
     !correlation.is_finite()) {
    Rcpp::stop("correlation must be a finite square matrix matching data.");
  }

  if(!arma::approx_equal(correlation, correlation.t(),
                         "absdiff", 1e-08)) {
    Rcpp::stop("correlation must be symmetric.");
  }

  if(arma::any(arma::abs(correlation.diag()-1.0) > 1e-08)) {
    Rcpp::stop("correlation must have a unit diagonal.");
  }

  if(thresholds.size() != static_cast<R_xlen_t>(nitems)) {
    Rcpp::stop("thresholds must contain one vector per observed variable.");
  }

  if(!std::isfinite(probability_floor) ||
     probability_floor <= 0.0 || probability_floor >= 1.0) {
    Rcpp::stop("probability_floor must be strictly between zero and one.");
  }

  if(!std::isfinite(inversion_tolerance) || inversion_tolerance <= 0.0) {
    Rcpp::stop("inversion_tolerance must be a positive finite number.");
  }

  // Reusable polyfast information

  const bool use_polyfast = polyfast_object.isNotNull();
  std::vector<double> category_min;
  std::vector<int> category_count;
  std::vector<std::vector<std::vector<int>>> cached_tables;

  if(use_polyfast) {

    Rcpp::List cache(polyfast_object.get());

    if(!cache.containsElementNamed("category_min") ||
       !cache.containsElementNamed("category_count") ||
       !cache.containsElementNamed("contingency_tables")) {
      Rcpp::stop("polyfast_object must contain category_min, category_count, and contingency_tables.");
    }

    category_min = Rcpp::as<std::vector<double>>(cache["category_min"]);
    category_count = Rcpp::as<std::vector<int>>(cache["category_count"]);
    cached_tables =
      Rcpp::as<std::vector<std::vector<std::vector<int>>>>(
        cache["contingency_tables"]
      );

    if(category_min.size() != nitems || category_count.size() != nitems) {
      Rcpp::stop("The polyfast category metadata does not match data.");
    }

  }

  // Thresholds and ordinal categories

  std::vector<arma::vec> tau(nitems);
  std::vector<arma::vec> levels(nitems);
  std::vector<arma::vec> bounds(nitems);
  std::vector<arma::vec> bounds_cdf(nitems);
  std::vector<arma::vec> category_probabilities(nitems);
  std::vector<arma::mat> marginal_scores(nitems);
  std::vector<arma::uvec> marginal_floored(nitems);

  arma::uvec threshold_counts(nitems, arma::fill::zeros);
  arma::uvec threshold_offsets(nitems, arma::fill::zeros);
  arma::umat categories(nobs, nitems, arma::fill::zeros);

  arma::uword nthresholds = 0L;

  for(arma::uword j = 0L; j < nitems; ++j) {

    tau[j] = finite_thresholds(thresholds[j], j);
    threshold_counts[j] = tau[j].n_elem;
    threshold_offsets[j] = nthresholds;
    nthresholds += tau[j].n_elem;

    if(use_polyfast) {

      if(category_count[j] < 2) {
        Rcpp::stop("Every variable in polyfast_object must contain at least two categories.");
      }

      const arma::uword ncategories =
        static_cast<arma::uword>(category_count[j]);

      if(ncategories != tau[j].n_elem+1L) {
        Rcpp::stop("The number of thresholds for variable " +
          std::to_string(j+1L) +
          " does not match the polyfast category metadata.");
      }

      for(arma::uword i = 0L; i < nobs; ++i) {

        const double shifted = data(i, j)-category_min[j];
        const double rounded = std::round(shifted);

        if(std::abs(shifted-rounded) > 1e-08 ||
           rounded < 0.0 || rounded >= category_count[j]) {
          Rcpp::stop("The data do not match the categories stored in polyfast_object.");
        }

        categories(i, j) = static_cast<arma::uword>(rounded);

      }

    } else {

      levels[j] = arma::sort(arma::unique(data.col(j)));

      if(levels[j].n_elem != tau[j].n_elem+1L) {
        Rcpp::stop("The number of thresholds for variable " +
          std::to_string(j+1L) +
          " does not match its observed categories. The threshold representation "
          "cannot contain empty internal categories.");
      }

      for(arma::uword i = 0L; i < nobs; ++i) {
        categories(i, j) = category_index(levels[j], data(i, j));
      }

    }

    bounds[j].set_size(tau[j].n_elem+2L);
    bounds[j][0L] = neg_inf;
    bounds[j].subvec(1L, tau[j].n_elem) = tau[j];
    bounds[j][bounds[j].n_elem-1L] = pos_inf;

    bounds_cdf[j].set_size(bounds[j].n_elem);

    for(arma::uword h = 0L; h < bounds[j].n_elem; ++h) {
      bounds_cdf[j][h] = normal_cdf_bound(bounds[j][h]);
    }

    category_probabilities[j] = arma::diff(bounds_cdf[j]);

    const arma::uword ncategories = tau[j].n_elem+1L;
    marginal_scores[j].zeros(ncategories, tau[j].n_elem);
    marginal_floored[j].zeros(ncategories);

    for(arma::uword category = 0L; category < ncategories; ++category) {

      double probability = category_probabilities[j][category];

      if(!std::isfinite(probability) || probability < probability_floor) {
        probability = probability_floor;
        marginal_floored[j][category] = 1L;
      }

      if(category > 0L) {
        marginal_scores[j](category, category-1L) -=
          Dnorm(tau[j][category-1L])/probability;
      }

      if(category < tau[j].n_elem) {
        marginal_scores[j](category, category) +=
          Dnorm(tau[j][category])/probability;
      }

    }

  }

  // Correlation pairs and response patterns

  const arma::uword ncorrelations = nitems*(nitems-1L)/2L;
  const arma::uword nparameters = nthresholds+ncorrelations;
  arma::umat pairs(ncorrelations, 2L, arma::fill::zeros);

  arma::uword pair_index = 0L;

  for(arma::uword j = 0L; j < nitems-1L; ++j) {
    for(arma::uword k = j+1L; k < nitems; ++k) {
      pairs(pair_index, 0L) = j;
      pairs(pair_index, 1L) = k;
      ++pair_index;
    }
  }

  if(use_polyfast && cached_tables.size() != ncorrelations) {
    Rcpp::stop("The number of polyfast contingency tables does not match data.");
  }

  std::map<std::vector<arma::uword>, arma::uword> pattern_map;

  for(arma::uword i = 0L; i < nobs; ++i) {

    std::vector<arma::uword> key(nitems);

    for(arma::uword j = 0L; j < nitems; ++j) {
      key[j] = categories(i, j);
    }

    ++pattern_map[key];

  }

  const arma::uword npatterns = pattern_map.size();
  arma::umat pattern_values(npatterns, nitems, arma::fill::zeros);
  arma::uvec pattern_weights(npatterns, arma::fill::zeros);

  arma::uword pattern_index = 0L;

  for(const auto& pattern : pattern_map) {

    for(arma::uword j = 0L; j < nitems; ++j) {
      pattern_values(pattern_index, j) = pattern.first[j];
    }

    pattern_weights[pattern_index] = pattern.second;
    ++pattern_index;

  }

  // Validate the entire pair cache before starting any worker. Rcpp error
  // reporting and access to R objects must remain on the main thread.

  for(arma::uword q = 0L; q < ncorrelations; ++q) {

    const arma::uword j = pairs(q, 0L);
    const arma::uword k = pairs(q, 1L);
    const arma::uword ncategories_j = tau[j].n_elem+1L;
    const arma::uword ncategories_k = tau[k].n_elem+1L;

    if(std::abs(correlation(j, k)) >= 1.0-1e-10) {
      Rcpp::stop("Every polychoric correlation must lie strictly inside (-1, 1).");
    }

    if(use_polyfast) {

      const std::vector<std::vector<int>>& table = cached_tables[q];

      if(table.size() != ncategories_j) {
        Rcpp::stop("A polyfast contingency table has an invalid number of rows.");
      }

      arma::uword total = 0L;

      for(arma::uword category_j = 0L;
          category_j < ncategories_j; ++category_j) {

        if(table[category_j].size() != ncategories_k) {
          Rcpp::stop("A polyfast contingency table has an invalid number of columns.");
        }

        for(arma::uword category_k = 0L;
            category_k < ncategories_k; ++category_k) {
          const int count = table[category_j][category_k];
          if(count < 0) {
            Rcpp::stop("polyfast contingency tables cannot contain negative counts.");
          }
          total += static_cast<arma::uword>(count);
        }

      }

      if(total != nobs) {
        Rcpp::stop("A polyfast contingency table does not match the sample size.");
      }

    }

  }

  // Pairwise cell scores and threshold cross-derivatives

  arma::mat A21(ncorrelations, nthresholds, arma::fill::zeros);
  std::vector<arma::mat> correlation_scores(ncorrelations);
  std::vector<arma::umat> correlation_floored(ncorrelations);

  int cores_used = acov_parallel_for(ncorrelations, cores,
    [&](const arma::uword q) {

    const arma::uword j = pairs(q, 0L);
    const arma::uword k = pairs(q, 1L);
    const arma::uword ncategories_j = tau[j].n_elem+1L;
    const arma::uword ncategories_k = tau[k].n_elem+1L;
    const double rho = correlation(j, k);

    arma::umat counts(ncategories_j, ncategories_k, arma::fill::zeros);

    if(use_polyfast) {

      const std::vector<std::vector<int>>& table = cached_tables[q];

      for(arma::uword category_j = 0L;
          category_j < ncategories_j; ++category_j) {
        for(arma::uword category_k = 0L;
            category_k < ncategories_k; ++category_k) {
          counts(category_j, category_k) =
            static_cast<arma::uword>(table[category_j][category_k]);
        }
      }

    } else {

      for(arma::uword i = 0L; i < nobs; ++i) {
        ++counts(categories(i, j), categories(i, k));
      }

    }

    correlation_scores[q].zeros(ncategories_j, ncategories_k);
    correlation_floored[q].zeros(ncategories_j, ncategories_k);

    for(arma::uword category_j = 0L;
        category_j < ncategories_j; ++category_j) {

      const double lower_j = bounds[j][category_j];
      const double upper_j = bounds[j][category_j+1L];

      for(arma::uword category_k = 0L;
          category_k < ncategories_k; ++category_k) {

        const arma::uword count = counts(category_j, category_k);
        if(count == 0L) continue;

        const double lower_k = bounds[k][category_k];
        const double upper_k = bounds[k][category_k+1L];

        double probability = pbinorm(
          rho,
          lower_j, lower_k,
          upper_j, upper_k,
          bounds_cdf[j][category_j],
          bounds_cdf[k][category_k],
          bounds_cdf[j][category_j+1L],
          bounds_cdf[k][category_k+1L]
        );

        if(!std::isfinite(probability) || probability < probability_floor) {
          probability = probability_floor;
          correlation_floored[q](category_j, category_k) = 1L;
        }

        const double probability_rho =
          dbinorm(rho, upper_j, upper_k)-
          dbinorm(rho, lower_j, upper_k)-
          dbinorm(rho, upper_j, lower_k)+
          dbinorm(rho, lower_j, lower_k);

        correlation_scores[q](category_j, category_k) =
          probability_rho/probability;

        const double weight = static_cast<double>(count)/
          static_cast<double>(nobs);

        auto add_first_derivative = [&](const arma::uword threshold,
                                        const double sign) {

          double probability_threshold = 0.0;
          double rho_threshold = 0.0;

          first_threshold_derivatives(
            rho, tau[j][threshold], sign,
            lower_k, upper_k,
            probability_threshold, rho_threshold
          );

          const double score_derivative =
            rho_threshold/probability-
            probability_rho*probability_threshold/
              (probability*probability);

          A21(q, threshold_offsets[j]+threshold) -=
            weight*score_derivative;

        };

        auto add_second_derivative = [&](const arma::uword threshold,
                                         const double sign) {

          double probability_threshold = 0.0;
          double rho_threshold = 0.0;

          second_threshold_derivatives(
            rho, tau[k][threshold], sign,
            lower_j, upper_j,
            probability_threshold, rho_threshold
          );

          const double score_derivative =
            rho_threshold/probability-
            probability_rho*probability_threshold/
              (probability*probability);

          A21(q, threshold_offsets[k]+threshold) -=
            weight*score_derivative;

        };

        if(category_j > 0L) {
          add_first_derivative(category_j-1L, -1.0);
        }

        if(category_j < tau[j].n_elem) {
          add_first_derivative(category_j, 1.0);
        }

        if(category_k > 0L) {
          add_second_derivative(category_k-1L, -1.0);
        }

        if(category_k < tau[k].n_elem) {
          add_second_derivative(category_k, 1.0);
        }

      }

    }

  });

  // Casewise estimating functions by unique response pattern

  arma::mat pattern_scores(npatterns, nparameters, arma::fill::zeros);
  arma::uvec marginal_floored_counts(nitems, arma::fill::zeros);
  arma::uvec correlation_floored_counts(ncorrelations, arma::fill::zeros);

  // Fill threshold-score columns contiguously. Each worker owns all columns
  // for one variable, so there is no shared write region.
  const int threshold_score_cores = acov_parallel_for(nitems, cores,
    [&](const arma::uword j) {

    const arma::uword first = threshold_offsets[j];

    for(arma::uword h = 0L; h < threshold_counts[j]; ++h) {

      double* score_column = pattern_scores.colptr(first+h);

      for(arma::uword r = 0L; r < npatterns; ++r) {
        const arma::uword category = pattern_values(r, j);
        score_column[r] = marginal_scores[j](category, h);
      }

    }

    arma::uword floored = 0L;

    for(arma::uword r = 0L; r < npatterns; ++r) {
      floored += marginal_floored[j][pattern_values(r, j)];
    }

    marginal_floored_counts[j] = floored;

  });

  cores_used = std::max(cores_used, threshold_score_cores);

  // Correlation scores are naturally independent by pair. Writing one full
  // column per worker follows Armadillo's column-major storage and avoids the
  // cache-unfriendly row-wise traversal used previously.
  const int correlation_score_cores = acov_parallel_for(ncorrelations, cores,
    [&](const arma::uword q) {

    const arma::uword j = pairs(q, 0L);
    const arma::uword k = pairs(q, 1L);
    double* score_column = pattern_scores.colptr(nthresholds+q);
    arma::uword floored = 0L;

    for(arma::uword r = 0L; r < npatterns; ++r) {

      const arma::uword category_j = pattern_values(r, j);
      const arma::uword category_k = pattern_values(r, k);

      score_column[r] =
        correlation_scores[q](category_j, category_k);

      floored += correlation_floored[q](category_j, category_k);

    }

    correlation_floored_counts[q] = floored;

  });

  cores_used = std::max(cores_used, correlation_score_cores);

  const arma::uword floored_probabilities =
    arma::accu(marginal_floored_counts)+
    arma::accu(correlation_floored_counts);

  arma::vec pattern_probabilities =
    arma::conv_to<arma::vec>::from(pattern_weights);
  pattern_probabilities /= static_cast<double>(nobs);

  arma::vec score_mean = pattern_scores.t()*pattern_probabilities;

  arma::mat weighted_scores = pattern_scores;
  weighted_scores.each_col() %= arma::sqrt(pattern_probabilities);

  // Dense crossproducts are intentionally left to Armadillo/BLAS. A tuned
  // GEMM is substantially more cache-efficient and vectorized than a custom
  // OpenMP dot-product kernel, and may use its own BLAS thread pool.
  arma::mat INNER = weighted_scores.t()*weighted_scores;
  INNER = 0.5*(INNER+INNER.t());

  // Lower-triangular sensitivity matrix

  arma::mat A11(nthresholds, nthresholds, arma::fill::zeros);

  for(arma::uword j = 0L; j < nitems; ++j) {

    const arma::uword first = threshold_offsets[j];
    const arma::uword last = first+threshold_counts[j]-1L;

    A11.submat(first, first, last, last) =
      INNER.submat(first, first, last, last);

  }

  A11 = 0.5*(A11+A11.t());

  arma::mat A22(ncorrelations, ncorrelations, arma::fill::zeros);

  for(arma::uword q = 0L; q < ncorrelations; ++q) {
    A22(q, q) = INNER(nthresholds+q, nthresholds+q);
  }

  arma::mat B(nparameters, nparameters, arma::fill::zeros);
  B.submat(0L, 0L, nthresholds-1L, nthresholds-1L) = A11;
  B.submat(nthresholds, 0L, nparameters-1L, nthresholds-1L) = A21;
  B.submat(nthresholds, nthresholds,
           nparameters-1L, nparameters-1L) = A22;

  // Block inverse

  arma::mat A11_inverse;
  bool generalized_A11 = !arma::inv_sympd(A11_inverse, A11);

  if(generalized_A11) {
    A11_inverse = arma::pinv(A11, inversion_tolerance);
  }

  arma::vec A22_diagonal = A22.diag();
  bool generalized_A22 = !A22_diagonal.is_finite() ||
    arma::any(arma::abs(A22_diagonal) <= inversion_tolerance);

  arma::mat A22_inverse(ncorrelations, ncorrelations, arma::fill::zeros);

  if(generalized_A22) {
    A22_inverse = arma::pinv(A22, inversion_tolerance);
  } else {
    A22_inverse.diag() = 1.0/A22_diagonal;
  }

  if(generalized_A11) {
    Rcpp::warning("A generalized inverse was used for the threshold-information block A11.");
  }

  if(generalized_A22) {
    Rcpp::warning("A generalized inverse was used for the correlation-information block A22.");
  }

  arma::mat B_inverse(nparameters, nparameters, arma::fill::zeros);
  B_inverse.submat(0L, 0L,
                   nthresholds-1L, nthresholds-1L) = A11_inverse;
  B_inverse.submat(nthresholds, nthresholds,
                   nparameters-1L, nparameters-1L) = A22_inverse;

  arma::mat B_inverse_lower = A21*A11_inverse;

  if(generalized_A22) {
    B_inverse_lower = -A22_inverse*B_inverse_lower;
  } else {
    B_inverse_lower.each_col() %= -A22_inverse.diag();
  }

  B_inverse.submat(nthresholds, 0L,
                   nparameters-1L, nthresholds-1L) = B_inverse_lower;

  // Influence scores and two-step sandwich

  arma::mat influence_thresholds =
    weighted_scores.cols(0L, nthresholds-1L)*A11_inverse.t();

  arma::mat influence_correlations =
    weighted_scores.cols(nthresholds, nparameters-1L)-
    influence_thresholds*A21.t();

  if(generalized_A22) {
    influence_correlations = influence_correlations*A22_inverse.t();
  } else {
    influence_correlations.each_row() %= A22_inverse.diag().t();
  }

  arma::mat influence_scores =
    arma::join_rows(influence_thresholds, influence_correlations);

  arma::mat NVCOV = influence_scores.t()*influence_scores;
  NVCOV = 0.5*(NVCOV+NVCOV.t());

  arma::mat VCOV = NVCOV/static_cast<double>(nobs);
  VCOV = 0.5*(VCOV+VCOV.t());

  if(floored_probabilities > 0L) {
    Rcpp::warning(std::to_string(floored_probabilities) +
      " univariate or bivariate response-pattern probabilities were replaced by probability_floor.");
  }

  // Result

  arma::uvec threshold_offsets_output = threshold_offsets+1L;
  arma::umat pairs_output = pairs+1L;

  Rcpp::List result;
  result["VCOV"] = VCOV;
  result["NVCOV"] = NVCOV;
  result["B"] = B;
  result["B_inverse"] = B_inverse;
  result["INNER"] = INNER;
  result["A11"] = A11;
  result["A21"] = A21;
  result["A22"] = A22;
  result["score_mean"] = score_mean;
  result["threshold_counts"] = threshold_counts;
  result["threshold_offsets"] = threshold_offsets_output;
  result["pairs"] = pairs_output;
  result["nobs"] = nobs;
  result["npatterns"] = npatterns;
  result["floored_probabilities"] = floored_probabilities;
  result["generalized_A11"] = generalized_A11;
  result["generalized_A22"] = generalized_A22;
  result["polyfast_reused"] = use_polyfast;
  // Maximum actual OpenMP team size (not the BLAS/LAPACK thread count).
  result["cores_requested"] = cores;
  result["cores_used"] = cores_used;
#ifdef _OPENMP
  result["openmp_available"] = true;
#else
  result["openmp_available"] = false;
#endif
  result["parameter_order"] =
    "finite thresholds by variable, followed by strict-lower-triangle correlations";

  if(return_scores) {
    arma::umat patterns_output = pattern_values+1L;
    result["patterns"] = patterns_output;
    result["pattern_weights"] = pattern_weights;
    result["pattern_scores"] = pattern_scores;
  }

  return result;

}
