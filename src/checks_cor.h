void check_cor(arguments_cor& x) {

  if(x.X.is_square()) {

    x.R = x.X;
    x.correlation_result["correlation"] = x.R;
    if(x.cor != "poly" & x.cor != "pearson") {
      Rcpp::stop("Unkown correlation");
    }

  } else {

    x.nobs = x.X.n_rows;
    missingness(x); // Handle missing values

    // ALLOW DWLS for continuous data
    if(x.cor == "poly") {
      if(x.estimator == "dwls") {
        x.correlation_result = polyfast(x.X, x.missing, "var", "none", 0.00, 0L, false, x.cores);
        arma::mat W = x.correlation_result["acov"];
        x.W = 1/W; x.W.diag().ones();
      } else {
        x.correlation_result = polyfast(x.X, x.missing, "none", "none", 0.00, 0L, false, x.cores);
      }
      arma::mat polys = x.correlation_result["correlation"];
      x.R = polys;
    } else if(x.cor == "pearson") {
      // Rcpp::Rcout << "1" << std::endl;
      if(x.X.has_nan()) {
        x.R = pairwise_cor(x.X);
      } else {
        x.R = arma::cor(x.X);
      }
      if(x.estimator == "dwls") {
        arma::vec asymp_diag;
        if(x.std_error == "normal") {
          asymp_diag = arma::diagvec(asymptotic_normal(x.R));
          x.correlation_result["std_error"] = "normal";
        } else {
          asymp_diag = arma::diagvec(asymptotic_general(x.X));
          x.correlation_result["std_error"] = "general";
        }
        arma::mat W = arma::reshape(asymp_diag, x.p, x.p);
        x.W = 1/W; x.W.diag().ones();
      }
      x.correlation_result["type"] = "pearson";
      x.correlation_result["correlation"] = x.R;
    } else {
      Rcpp::stop("Unkown correlation");
    }

  }

}
