/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 25/08/2026
 */

Rcpp::List get_jacob(Rcpp::List control_manifold,
                     Rcpp::List control_transform,
                     Rcpp::List control_estimator,
                     Rcpp::List control_optimizer) {

  (void)control_manifold;
  (void)control_estimator;

  arguments_optim x;
  x.ntransforms = control_transform.size();

  // Parameter coordinates

  if(!control_optimizer.containsElementNamed("parameters") ||
     !control_optimizer.containsElementNamed("transparameters") ||
     !control_optimizer.containsElementNamed("transparam2param")) {
    Rcpp::stop("The optimizer control does not contain the parameter mappings "
               "required to compute the Jacobian.");
  }

  std::vector<arma::vec> parameters_list =
    control_optimizer["parameters"];
  std::vector<arma::vec> transparameters_list =
    control_optimizer["transparameters"];

  if(parameters_list.size() == 0L ||
     transparameters_list.size() == 0L) {
    Rcpp::stop("The fitted parameter vectors are missing.");
  }

  x.parameters = parameters_list[0L];
  x.transparameters = transparameters_list[0L];
  x.transparameters_init = x.transparameters;
  x.transparam2param =
    Rcpp::as<arma::uvec>(control_optimizer["transparam2param"]);
  x.nparam = x.parameters.n_elem;
  x.ntransparam = x.transparameters.n_elem;

  if(x.transparam2param.n_elem !=
     static_cast<arma::uword>(x.nparam)) {
    Rcpp::stop("transparam2param has incompatible dimensions.");
  }

  if(x.transparam2param.n_elem > 0L &&
     x.transparam2param.max() >=
       static_cast<arma::uword>(x.ntransparam)) {
    Rcpp::stop("transparam2param contains an index outside the transformed "
               "parameter vector.");
  }

  x.transparameters(x.transparam2param) =
    x.parameters;
  x.transparameters_init = x.transparameters;

  if(control_optimizer.containsElementNamed("idx_transforms")) {
    x.idx_transforms =
      Rcpp::as<arma::uvec>(control_optimizer["idx_transforms"]);
  } else {
    x.idx_transforms =
      arma::regspace<arma::uvec>(
        0L,
        x.ntransforms > 0L ?
          static_cast<arma::uword>(x.ntransforms-1L) :
          0L
      );

    if(x.ntransforms == 0L) {
      x.idx_transforms.reset();
    }
  }

  // Transformation objects

  std::vector<transformations*> xtransforms(x.ntransforms);

  for(int i = 0L; i < x.ntransforms; ++i) {
    xtransforms[i] =
      choose_transform(control_transform[i]);
  }

  product_transform final_transform;

  // Evaluate every transformation once so each local transformation object
  // contains the fitted quantities needed by its jacobian() method.
  final_transform.transform(x, xtransforms);

  // Full transparameter dependency Jacobian

  // Every transformed parameter is represented as a coordinate, so start from
  // the identity. Transformation rows are then replaced by their direct and
  // cumulative dependencies while retaining their own diagonal identity.
  x.jacob.eye(x.ntransparam, x.ntransparam);

  for(arma::uword t : x.idx_transforms) {

    if(t >= xtransforms.size()) {
      Rcpp::stop("A requested transformation index is outside the "
                 "transformation list.");
    }

    arma::uvec indices_in =
      xtransforms[t]->indices_in;
    arma::uvec indices_out =
      xtransforms[t]->indices_out;

    if((indices_in.n_elem > 0L &&
        indices_in.max() >=
          static_cast<arma::uword>(x.ntransparam)) ||
       (indices_out.n_elem > 0L &&
        indices_out.max() >=
          static_cast<arma::uword>(x.ntransparam))) {
      Rcpp::stop("A transformation contains parameter indices outside the "
                 "transformed-parameter vector.");
    }

    xtransforms[t]->jacobian(x);
    const arma::mat& local_jacobian =
      xtransforms[t]->jacob;

    if(local_jacobian.n_rows != indices_out.n_elem ||
       local_jacobian.n_cols != indices_in.n_elem) {
      Rcpp::stop(
        "The Jacobian dimensions of transformation " +
        std::to_string(t+1L) +
        " do not match its input and output indices."
      );
    }

    if(!local_jacobian.is_finite()) {
      Rcpp::stop(
        "Transformation " +
        std::to_string(t+1L) +
        " returned a non-finite Jacobian."
      );
    }

    // The current rows of the inputs already contain their direct and
    // transitive dependencies.
    arma::mat input_jacobian(
      indices_in.n_elem,
      x.ntransparam,
      arma::fill::zeros
    );

    for(arma::uword i = 0L;
        i < indices_in.n_elem;
        ++i) {

      for(arma::sp_mat::const_row_iterator it =
            x.jacob.begin_row(indices_in[i]);
          it != x.jacob.end_row(indices_in[i]);
          ++it) {
        input_jacobian(i, it.col()) = *it;
      }

    }

    arma::mat output_jacobian =
      local_jacobian*input_jacobian;

    // Replace each output row with its dependency row. The diagonal identity is
    // retained because the full matrix is indexed by all transformed
    // parameters, not only by the freely estimated parameters.
    for(arma::uword i = 0L;
        i < indices_out.n_elem;
        ++i) {

      arma::uword output_index =
        indices_out[i];

      x.jacob.row(output_index).zeros();
      x.jacob(output_index, output_index) = 1.0;

      for(arma::uword j = 0L;
          j < static_cast<arma::uword>(x.ntransparam);
          ++j) {

        if(j == output_index) {
          continue;
        }

        double value = output_jacobian(i, j);

        if(value != 0.0) {
          x.jacob(output_index, j) = value;
        }

      }

    }

  }

  // Result

  Rcpp::List result;
  result["jacob"] = x.jacob;

  return result;

}
