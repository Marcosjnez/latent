/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 01/09/2026
 */

arma::uvec check_poblq_blocks_oblique(SEXP oblique_input,
                                      arma::uword p) {

  Rcpp::NumericVector input(oblique_input);

  if(input.size() == 0L) {
    Rcpp::stop("The poblq_blocks manifold requires at least one oblique block.");
  }

  arma::uvec oblique(input.size());
  arma::uword total = 0L;

  for(R_xlen_t i=0L; i < input.size(); ++i) {

    double value = input[i];

    if(!std::isfinite(value) ||
       value <= 0.00 ||
       value != std::floor(value)) {
      Rcpp::stop("The poblq_blocks oblique block sizes must be positive integers.");
    }

    if(value > static_cast<double>(p-total)) {
      Rcpp::stop("The poblq_blocks oblique block sizes cannot sum to more than p.");
    }

    oblique[i] = static_cast<arma::uword>(value);
    total += oblique[i];

  }

  return oblique;

}

arma::mat poblq_blocks_constraints(const arma::uvec& oblique,
                                    arma::uword p) {

  arma::mat constraints(p, p, arma::fill::zeros);
  arma::uword start = 0L;

  for(arma::uword block=0L; block < oblique.n_elem; ++block) {

    arma::uword end = start+oblique[block]-1L;
    constraints.submat(start, start, end, end).ones();

    for(arma::uword i=start; i <= end; ++i) {
      constraints(i, i) = 0.00;
    }

    start = end+1L;

  }

  return constraints;

}

arma::mat poblq_blocks_retract(arma::mat X,
                               const arma::uvec& oblique) {

  if(X.n_rows == 0L ||
     X.n_rows != X.n_cols ||
     !X.is_finite()) {
    Rcpp::stop("The poblq_blocks manifold requires a non-empty finite square matrix.");
  }

  const arma::uword p = X.n_cols;
  const double tolerance = std::sqrt(arma::datum::eps);
  arma::mat result(p, p, arma::fill::zeros);
  arma::mat projector(p, p, arma::fill::zeros);
  arma::uword start = 0L;

  for(arma::uword block=0L; block < oblique.n_elem; ++block) {

    arma::uword end = start+oblique[block]-1L;
    arma::mat R = X.cols(start, end)-projector*X.cols(start, end);
    arma::rowvec norms = arma::sqrt(arma::sum(R % R, 0));

    if(!norms.is_finite() || arma::any(norms <= tolerance)) {
      Rcpp::stop("The poblq_blocks retraction produced a zero column in an oblique block.");
    }

    R.each_row() /= norms;
    result.cols(start, end) = R;

    arma::mat gram_inverse;
    bool inverted = arma::inv_sympd(gram_inverse, R.t()*R);

    if(!inverted || !gram_inverse.is_finite()) {
      Rcpp::stop("The poblq_blocks retraction requires full-rank oblique blocks.");
    }

    projector += R*gram_inverse*R.t();
    projector = 0.50*(projector+projector.t());
    start = end+1L;

  }

  if(start < p) {

    arma::mat R = X.cols(start, p-1L)-projector*X.cols(start, p-1L);
    arma::mat eigvec;
    arma::vec eigval;
    bool decomposed = arma::eig_sym(eigval, eigvec, R.t()*R);

    if(!decomposed ||
       !eigval.is_finite() ||
       !eigvec.is_finite() ||
       eigval.min() <= tolerance*std::max(1.00, eigval.max())) {
      Rcpp::stop("The poblq_blocks retraction requires a full-rank orthogonal block.");
    }

    R *= eigvec*arma::diagmat(1.00/arma::sqrt(eigval))*eigvec.t();
    result.cols(start, p-1L) = R;

  }

  if(!result.is_finite()) {
    Rcpp::stop("The poblq_blocks retraction produced non-finite values.");
  }

  return result;

}

// Partially oblique manifold with blocks:

class poblq_blocks:public poblq {

public:

  arma::uvec oblique, block_starts;
  std::vector<arma::vec> block_eigval;
  std::vector<arma::mat> block_eigvec;

  arma::mat solve_blocks_sylvester(const arma::mat& C,
                                   arma::uword left,
                                   arma::uword right) {

    const arma::vec& left_values = block_eigval[left];
    const arma::vec& right_values = block_eigval[right];
    const arma::mat& left_vectors = block_eigvec[left];
    const arma::mat& right_vectors = block_eigvec[right];

    arma::mat denominator =
      arma::repmat(left_values, 1L, right_values.n_elem)+
      arma::repmat(right_values.t(), left_values.n_elem, 1L);

    const double tolerance = std::sqrt(arma::datum::eps);

    if(!denominator.is_finite() ||
       arma::any(arma::vectorise(denominator) <= tolerance)) {
      Rcpp::stop("The poblq_blocks Sylvester equation is not uniquely defined.");
    }

    arma::mat solution =
      left_vectors*
      ((left_vectors.t()*C*right_vectors)/denominator)*
      right_vectors.t();

    if(!solution.is_finite()) {
      Rcpp::stop("The poblq_blocks Sylvester equation could not be solved.");
    }

    return solution;

  }

  arma::mat normal_coefficients(const arma::mat& G) {

    arma::mat coefficients(p, p, arma::fill::zeros);
    const arma::uword nblocks = oblique.n_elem;

    for(arma::uword i=0L; i < nblocks; ++i) {

      arma::uword start_i = block_starts[i];
      arma::uword end_i = start_i+oblique[i]-1L;
      arma::mat Xi = X.cols(start_i, end_i);
      arma::mat Gi = G.cols(start_i, end_i);

      coefficients.submat(start_i, start_i, end_i, end_i) =
        arma::diagmat(Xi.t()*Gi);

      for(arma::uword j=i+1L; j < nblocks; ++j) {

        arma::uword start_j = block_starts[j];
        arma::uword end_j = start_j+oblique[j]-1L;
        arma::mat Xj = X.cols(start_j, end_j);
        arma::mat Gj = G.cols(start_j, end_j);
        arma::mat C = Xj.t()*Gi+Gj.t()*Xi;
        arma::mat B = solve_blocks_sylvester(C, j, i);

        coefficients.submat(start_j, start_i, end_j, end_i) = B;
        coefficients.submat(start_i, start_j, end_i, end_j) = B.t();

      }

    }

    arma::uword orthogonal_start = arma::sum(oblique);

    if(orthogonal_start < p) {

      arma::mat Xo = X.cols(orthogonal_start, p-1L);
      arma::mat Go = G.cols(orthogonal_start, p-1L);

      coefficients.submat(orthogonal_start, orthogonal_start,
                          p-1L, p-1L) =
        0.50*(Xo.t()*Go+Go.t()*Xo);

      for(arma::uword i=0L; i < nblocks; ++i) {

        arma::uword start_i = block_starts[i];
        arma::uword end_i = start_i+oblique[i]-1L;
        arma::mat Xi = X.cols(start_i, end_i);
        arma::mat Gi = G.cols(start_i, end_i);
        arma::mat C = Xo.t()*Gi+Go.t()*Xi;
        arma::mat B = C*block_eigvec[i];
        arma::rowvec inverse =
          (1.00/(block_eigval[i]+1.00)).t();

        B.each_row() %= inverse;
        B *= block_eigvec[i].t();

        coefficients.submat(orthogonal_start, start_i,
                            p-1L, end_i) = B;
        coefficients.submat(start_i, orthogonal_start,
                            end_i, p-1L) = B.t();

      }

    }

    return coefficients;

  }

  arma::mat project(const arma::mat& G) {

    arma::mat coefficients = normal_coefficients(G);

    return G-X*coefficients;

  }

  void param(arguments_optim& x) {

    poblq::param(x);

    const arma::uword nblocks = oblique.n_elem;
    const double tolerance = std::sqrt(arma::datum::eps);
    block_eigval.resize(nblocks);
    block_eigvec.resize(nblocks);

    for(arma::uword i=0L; i < nblocks; ++i) {

      arma::uword start = block_starts[i];
      arma::uword end = start+oblique[i]-1L;
      arma::mat Y = psi.submat(start, start, end, end);
      bool decomposed = arma::eig_sym(block_eigval[i],
                                      block_eigvec[i], Y);

      if(!decomposed ||
         !block_eigval[i].is_finite() ||
         !block_eigvec[i].is_finite() ||
         block_eigval[i].min() <=
         tolerance*std::max(1.00, block_eigval[i].max())) {
        Rcpp::stop("The poblq_blocks projection requires full-rank oblique blocks.");
      }

    }

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, p);
    A = normal_coefficients(g);

    x.rg.elem(indices) =
      arma::vectorise(g-X*A);

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, p);
    dg = arma::reshape(x.dg.elem(indices), p, p);
    dX = arma::reshape(x.dparameters.elem(indices), p, p);
    A = normal_coefficients(g);

    arma::mat h = dg-dX*A;

    x.dH.elem(indices) =
      arma::vectorise(project(h));

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) =
      arma::vectorise(poblq_blocks_retract(X, oblique));

  }

};

poblq_blocks* choose_poblq_blocks(Rcpp::List manifold_setup) {

  poblq_blocks* mymanifold = new poblq_blocks();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];

  if(p == 0L) {
    Rcpp::stop("The poblq_blocks manifold requires a positive p dimension.");
  }

  arma::uvec oblique = check_poblq_blocks_oblique(
    manifold_setup["oblique"],
    static_cast<arma::uword>(p)
  );
  arma::mat constraints = poblq_blocks_constraints(
    oblique,
    static_cast<arma::uword>(p)
  );
  arma::umat constraint_pairs =
    poblq_constraint_pairs(constraints);
  arma::mat constraint_basis =
    poblq_constraint_basis(
      constraint_pairs,
      static_cast<arma::uword>(p)
    );
  arma::uvec block_starts(oblique.n_elem);
  arma::uword start = 0L;

  for(arma::uword i=0L; i < oblique.n_elem; ++i) {
    block_starts[i] = start;
    start += oblique[i];
  }

  mymanifold->indices = indices;
  mymanifold->constraints = constraints;
  mymanifold->constraint_pairs = constraint_pairs;
  mymanifold->constraint_basis = constraint_basis;
  mymanifold->p = p;
  mymanifold->q = p;
  mymanifold->oblq_indices = arma::find(constraints == 1.00);
  mymanifold->oblique = oblique;
  mymanifold->block_starts = block_starts;

  return mymanifold;

}
