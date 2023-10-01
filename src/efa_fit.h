/*
 * Author: Marcos Jimenez
 * email: marcosjnezhquez@gmail.com
 * Modification date: 18/03/2022
 *
 */

arma::vec sdp_cpp(arma::mat Cov) {

  Rcpp::Environment Rcsdp("package:Rcsdp");
  Rcpp::Function csdp = Rcsdp["csdp"];

  int p = Cov.n_rows;
  arma::vec LoBounds(p, arma::fill::zeros);
  arma::vec UpBounds(p, arma::fill::ones);
  arma::vec Var = UpBounds;
  arma::vec opt = UpBounds;

  arma::mat Cov2 = -Cov;
  Cov2.diag().zeros();

  Rcpp::List C = Rcpp::List::create(Cov2, -UpBounds, LoBounds);
  Rcpp::List A(p);
  for (int i=0; i < p; ++i) {
    arma::vec b(p, arma::fill::zeros);
    b(i) = 1;
    A[i] = Rcpp::List::create(diagmat(b), -b, b);
  }

  arma::vec size_K(3, arma::fill::ones);
  size_K *= p;
  Rcpp::CharacterVector type(3);
  type[0] = "s";
  type[1] = "l";
  type[2] = "l";
  Rcpp::List K = Rcpp::List::create(Rcpp::_["type"] = type, Rcpp::_["size"] = size_K);

  Rcpp::Function csdp_control = Rcsdp["csdp.control"];
  Rcpp::List control = csdp_control();
  control["printlevel"] = 0;

  Rcpp::List result = csdp(C, A, opt, K, control);
  arma::vec item_diag = result["y"];

  return item_diag;

}

Rcpp::List principal_axis(arma::vec psi, arma::mat R, int nfactors,
                    double rel_tol, int efa_max_iter) {

  Rcpp::List result;

  double criteria;
  arma::mat w, ww;
  int iteration = 0;

  do {

    iteration = iteration + 1;

    arma::mat reduced_R = R - diagmat(psi);

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, reduced_R);

    arma::vec eigval2 = reverse(eigval);
    arma::mat eigvec2 = reverse(eigvec, 1);

    arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
    arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
    for(int i=0; i < nfactors; ++i) {
      if(eigenvalues[i] < 0) eigenvalues[i] = 0;
    }
    arma::mat D = diagmat(sqrt(eigenvalues));

    w = A * D;
    ww = w * w.t();

    arma::vec new_psi = 1 - diagvec(ww);

    criteria = arma::accu( abs(psi - new_psi) );
    psi = new_psi;

  } while (criteria > rel_tol && iteration < efa_max_iter);

  bool convergence = false;
  if(criteria < rel_tol) {
    convergence = true;
  }
  arma::mat Rhat = ww;
  Rhat.diag() = R.diag();
  arma::mat residuals = R - Rhat;
  double f = arma::accu(residuals % residuals);

  result["convergence"] = convergence;
  result["f"] = f;
  result["loadings"] = w;
  result["uniquenesses"] = psi;
  result["Rhat"] = Rhat;
  result["residuals"] = residuals;
  result["iterations"] = iteration;

  return result;

}

double ml_objective(arma::vec psi, arma::mat R, int nfactors, int n_items) {

  arma::mat sc = diagmat(1/sqrt(psi));
  arma::mat Sstar = sc * R * sc;

  arma::vec eigval;
  eig_sym(eigval, Sstar);

  arma::vec e = eigval(arma::span(0, n_items - nfactors - 1));

  // double objective = -arma::accu(log(e) + 1/e - 1);
  double objective = arma::accu(log(e) - e) + n_items - nfactors;

  return -objective;

}

arma::vec ml_gradient(arma::vec psi, arma::mat R, int nfactors, int n_items) {

  arma::vec sqrt_psi = sqrt(psi);
  arma::mat sc = diagmat(1/sqrt_psi);
  arma::mat Sstar = sc * R * sc;

  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, Sstar);

  arma::vec eigval2 = reverse(eigval);
  arma::mat eigvec2 = reverse(eigvec, 1);

  arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
  arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
  // for(int i=0; i < nfactors; ++i) {
  //   if(eigenvalues[i] < 0) eigenvalues[i] = 0;
  // }
  // arma::mat D = diagmat(sqrt(eigenvalues - 1));
  // arma::mat w = A * D;
  // w = diagmat(sqrt_psi) * w;
  // arma::mat ww = w * w.t();
  // arma::mat residuals = R - ww - diagmat(psi);
  // arma::mat gradient = -diagvec(residuals) / (psi % psi);

  arma::mat gradient = ((A % A) * (eigenvalues - 1) + 1 -
    arma::diagvec(R)/psi)/psi;

  return gradient;
}

double uls_objective(arma::vec psi, arma::mat R, int nfactors) {

  int n_items = psi.size();
  arma::mat reduced_R = R - diagmat(psi);

  arma::vec eigval;
  // arma::mat eigvec;
  // eig_sym(eigval, eigvec, reduced_R);
  //
  // arma::vec eigval2 = reverse(eigval);
  // arma::mat eigvec2 = reverse(eigvec, 1);
  //
  // arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
  // arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
  // for(int i=0; i < nfactors; ++i) {
  //   if(eigenvalues[i] < 0) eigenvalues[i] = 0;
  // }
  // arma::mat D = diagmat(sqrt(eigenvalues));
  //
  // arma::mat w = A * D;
  // arma::mat ww = w * w.t();
  // arma::mat residuals = R - ww - diagmat(psi);
  //
  // double objective = 0.5*arma::accu(residuals % residuals);

  eig_sym(eigval, reduced_R);
  arma::vec e = eigval(arma::span(0, n_items - nfactors - 1));
  double objective = arma::accu(e % e);

  return objective;

}

arma::vec uls_gradient(arma::vec psi, arma::mat R, int nfactors) {

  int n_items = psi.size();
  arma::mat reduced_R = R - diagmat(psi);

  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, reduced_R);

  // arma::vec eigval2 = reverse(eigval);
  // arma::mat eigvec2 = reverse(eigvec, 1);
  //
  // arma::mat A = eigvec2(arma::span::all, arma::span(0, nfactors-1));
  // arma::vec eigenvalues = eigval2(arma::span(0, nfactors-1));
  // for(int i=0; i < nfactors; ++i) {
  //   if(eigenvalues[i] < 0) eigenvalues[i] = 0;
  // }
  // arma::mat D = diagmat(sqrt(eigenvalues));
  //
  // arma::mat w = A * D;
  // arma::mat ww = w * w.t();
  // arma::mat residuals = R - ww - diagmat(psi);
  // arma::mat gradient = -diagvec(residuals);

  arma::vec e_values = eigval(arma::span(0, n_items - nfactors - 1));
  arma::mat e_vectors = eigvec(arma::span::all, arma::span(0, n_items - nfactors - 1));
  arma::mat gradient = -2*arma::diagvec(e_vectors * arma::diagmat(e_values) * e_vectors.t());

  return gradient;

}

Rcpp::List optim_rcpp(arma::vec psi, arma::mat R, int nfactors, std::string method,
                int efa_max_iter, double efa_factr, int m) {

  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  int n_items = psi.size();

  Rcpp::List results;
  Rcpp::List control;
  control["maxit"] = efa_max_iter;
  control["factr"] = efa_factr;
  control["lmm"] = m;
  arma::vec parscale(n_items, arma::fill::ones);
  parscale *= 0.01;
  control["parscale"] = parscale;

  if (method == "uls") {
    results = optim(Rcpp::_["par"] = psi,
                    Rcpp::_["fn"] = Rcpp::InternalFunction(&uls_objective),
                    Rcpp::_["gr"] = Rcpp::InternalFunction(&uls_gradient),
                    Rcpp::_["method"] = "L-BFGS-B",
                    Rcpp::_["lower"] = 0.005,
                    Rcpp::_["upper"] = arma::diagvec(R),
                    Rcpp::_["control"] = control,
                    Rcpp::_["R"] = R,
                    Rcpp::_["nfactors"] = nfactors);
  } else if (method == "ml") {

    results = optim(Rcpp::_["par"] = psi,
                    Rcpp::_["fn"] = Rcpp::InternalFunction(&ml_objective),
                    Rcpp::_["gr"] = Rcpp::InternalFunction(&ml_gradient),
                    Rcpp::_["method"] = "L-BFGS-B",
                    Rcpp::_["lower"] = 0.005,
                    Rcpp::_["upper"] = arma::diagvec(R),
                    Rcpp::_["control"] = control,
                    Rcpp::_["R"] = R,
                    Rcpp::_["nfactors"] = nfactors,
                    Rcpp::_["n_items"] = n_items);

  }

  return results;

}
