void check_cfa(arguments_cfa& x) {

  Rcpp::List cfa_control;
  // Check control parameters:
  if (x.nullable_control.isNotNull()) {
    cfa_control = Rcpp::as<Rcpp::List>(x.nullable_control);
  }

  if(x.p < x.q) Rcpp::stop("Too many factors");
  x.df = arma::accu(x.p*(x.p+1)/2) - x.parameters.size();
  x.df_null = x.p*(x.p-1)/2;

  // Choose custom weight matrix for the dwls estimator:
  if(cfa_control.containsElementNamed("W")) {
    arma::mat W = cfa_control["W"];
    if(W.n_cols != x.p | W.n_rows != x.p) {
      Rcpp::stop("W must be a matrix with the same dimensions as the correlation matrix");
    }
    x.W = W;
  }

  // If the raw data was not provided and the dwls estimator was selected, then
  // compute the normal approximation to the asymptotic correlation matrix:
  if(x.X.is_square() & x.estimator == "dwls" & !cfa_control.containsElementNamed("W")) {
    arma::vec asymp_diag;
    if(x.std_error == "normal") {
      asymp_diag = arma::diagvec(asymptotic_normal(x.R));
      x.correlation_result["std_error"] = "normal";
    } else {
      Rcpp::stop("estimator = 'dwls' requires either the raw data or a weight matrix W of the same dimension as the correlation matrix (control = list(W = ...))");
    }
    arma::mat W = arma::reshape(asymp_diag, x.p, x.p);
    x.W = 1/W; x.W.diag().ones();
  }

  if(x.estimator == "dwls") {
    arma::mat residuals = x.R; residuals.diag().zeros();
    x.f_null = 0.5*arma::accu(residuals % residuals % x.W);
  }

  // If estimator = "uls", then W is just filled with ones:
  if(x.estimator == "uls") {
    arma::mat W(x.p, x.p, arma::fill::ones);
    x.W = W;
    x.f_null = 0.5*(arma::accu(x.R % x.R) - x.p); // Assumes diag of x.R is identity
  }

  // estimator = "gls" requires a weight matrix W:
  if(x.estimator == "gls") {
    if(x.W.is_empty()) {
      Rcpp::stop("For the gls estimator, please provide a weight matrix in efa_control = list(W = ...)");
    }
    x.f_null = 0.5*(arma::accu(x.R % x.R) - x.p);
  }

  if(x.estimator == "ml") {
    x.logdetR = arma::log_det_sympd(x.R);
    x.f_null = -x.logdetR;
  }

  if(cfa_control.containsElementNamed("upper")) {
    arma::vec upper = cfa_control["upper"];
    x.upper = upper;
  }

  if(cfa_control.containsElementNamed("lower")) {
    arma::vec lower = cfa_control["lower"];
    x.lower = lower;
  }

  if(cfa_control.containsElementNamed("target_positive")) {
    arma::uvec target_positive = cfa_control["target_positive"];
    x.target_positive = target_positive-1;
  }

  // Arrange stuff for CFA:
  x.Ip.set_size(x.p, x.p); x.Ip.eye();
  x.Iq.set_size(x.q, x.q); x.Iq.eye();
  x.w = arma::vectorise(x.W);
  x.indexes_diag_q.set_size(x.q);
  for(int i=0; i < x.q; ++i) x.indexes_diag_q[i] = i * x.q + i;
  x.indexes_diag_p.set_size(x.p);
  for(int i=0; i < x.p; ++i) x.indexes_diag_p[i] = i * x.p + i;
  x.indexes_diag_q2.set_size(x.q);
  x.indexes_diag_q2[0] = 0;
  for(int i=1; i < x.q; ++i) x.indexes_diag_q2[i] = x.indexes_diag_q2[i-1] + (x.q-i+2);
  x.indexes_p = arma::trimatl_ind(arma::size(x.psi));
  x.indexes_q = arma::trimatl_ind(arma::size(x.phi));
  x.S_indexes = x.indexes_p;
  x.Rhat_inv = x.Ip;
  x.targetT_indexes = consecutive(0L, x.q*x.q-1L);
  if(x.projection == "positive") x.positive = false;

  x.dlambda.set_size(x.p, x.q); x.dlambda.zeros();
  x.dphi.set_size(x.q, x.q); x.dphi.zeros();
  x.dpsi.set_size(x.p, x.p); x.dpsi.zeros();
  x.T.set_size(x.q, x.q);
  x.dT.set_size(x.q, x.q);
  x.Phi_Target = x.phi;
  x.oblq_indexes = arma::find(x.Phi_Target == 1);

  x.U.set_size(x.p, x.p);
  x.dU.set_size(x.p, x.p);
  x.Psi_Target = x.psi;
  x.psi_oblq_indexes = arma::find(x.Psi_Target == 1);

}

void check_opt(arguments_optim& x) {

  Rcpp::List cfa_control;
  // Check control parameters:
  if (x.nullable_control.isNotNull()) {
    cfa_control = Rcpp::as<Rcpp::List>(x.nullable_control);
  }

  if(cfa_control.containsElementNamed("optim")) {
    std::string optim = cfa_control["optim"];
    x.optim = optim;
  }

  if(cfa_control.containsElementNamed("search")) {
    std::string search = cfa_control["search"];
    x.search = search;
    if(x.search != "back" & x.search != "wolfe") {
      Rcpp::stop("Unkown line-search method. Available methods: back and wolfe");
    }
  } else {
    x.search = "back";
  }

  if(cfa_control.containsElementNamed("maxit") ){
    x.maxit = cfa_control["maxit"];
  } else {
    x.maxit = 1e3;
  }

  // if(cfa_control.containsElementNamed("upper")) {
  //   arma::vec upper = cfa_control["upper"];
  //   x.upper = upper;
  // }
  //
  // if(cfa_control.containsElementNamed("lower")) {
  //   arma::vec lower = cfa_control["lower"];
  //   x.lower = lower;
  // }

  if(cfa_control.containsElementNamed("c1")) {
    x.c1 = cfa_control["c1"];
    if(x.c1 <= 0 | x.c1 >= 1) Rcpp::stop("The parameter c1 must be a positive number lower than 1.");
  } else {
    if(x.optim == "gradient") {
      x.c1 = 0.5;
    } else {
      x.c1 = 10e-04;
    }
  }

  if(cfa_control.containsElementNamed("c2")) {
    x.c2 = cfa_control["c2"];
    if(x.c2 < x.c1 | x.c2 >= 1) Rcpp::stop("The parameter c2 must be a positive number between c1 and 1.");
  } else {
    x.c2 = std::max(x.c1, 0.5);
  }

  if(cfa_control.containsElementNamed("rho")) {
    x.rho = cfa_control["rho"];
    if(x.rho <= 0 | x.rho >= 1) Rcpp::stop("The parameter rho must be a positive number between 0 and 1.");
  }

  if(cfa_control.containsElementNamed("M")) {
    x.M = cfa_control["M"];
  }

  if(cfa_control.containsElementNamed("armijo_maxit")) {
    x.armijo_maxit = cfa_control["armijo_maxit"];
  }

  if(x.random_starts < 1) {
    Rcpp::stop("The number of random starts should be a positive integer");
  }
  if(x.cores < 1) {
    Rcpp::stop("The number of cores should be a positive integer");
  }

}
