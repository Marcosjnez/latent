/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 27/08/2026
 */

arma::mat lyap_sym(arma::mat Y, arma::mat Q) {

  // Solve the lyapunov equation YX + XY = Q with symmetric Q and X:

  int q = Y.n_cols;
  arma::vec I(q, arma::fill::ones);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Y);

  arma::mat M = eigvec.t() * Q * eigvec;
  arma::mat W1 = I * eigval.t();
  arma::mat W = W1 + W1.t();
  arma::mat YY = M / W;
  arma::mat A = eigvec * YY * eigvec.t();

  return A;

}

arma::mat check_poblq_constraints(arma::mat constraints,
                                  arma::uword q) {

  if(constraints.n_rows != q || constraints.n_cols != q) {
    Rcpp::stop("The poblq constraints must be a q by q matrix.");
  }

  if(!constraints.is_finite()) {
    Rcpp::stop("The poblq constraints cannot contain missing or non-finite values.");
  }

  const double tolerance = std::sqrt(arma::datum::eps);

  if(!arma::approx_equal(constraints, constraints.t(),
                         "absdiff", tolerance)) {
    Rcpp::stop("The poblq constraints must be symmetric.");
  }

  for(arma::uword i=0L; i < constraints.n_elem; ++i) {

    if(std::abs(constraints[i]) <= tolerance) {
      constraints[i] = 0.00;
    } else if(std::abs(constraints[i]-1.00) <= tolerance) {
      constraints[i] = 1.00;
    } else {
      Rcpp::stop("The poblq constraints must contain only zeros and ones.");
    }

  }

  return constraints;

}

arma::umat poblq_constraint_pairs(const arma::mat& constraints) {

  const arma::uword q = constraints.n_cols;
  arma::uword nconstraints = 0L;

  for(arma::uword column=0L; column < q; ++column) {
    for(arma::uword row=column; row < q; ++row) {
      if(constraints(column, row) == 0.00) {
        nconstraints++;
      }
    }
  }

  arma::umat constraint_pairs(2L, nconstraints);
  arma::uword constraint = 0L;

  for(arma::uword column=0L; column < q; ++column) {
    for(arma::uword row=column; row < q; ++row) {

      if(constraints(column, row) == 0.00) {
        constraint_pairs(0L, constraint) = column;
        constraint_pairs(1L, constraint) = row;
        constraint++;
      }

    }
  }

  return constraint_pairs;

}

arma::mat poblq_constraint_basis(const arma::umat& constraint_pairs,
                                 arma::uword q) {

  const arma::uword nconstraints = constraint_pairs.n_cols;
  arma::mat constraint_basis(q*q, nconstraints,
                             arma::fill::zeros);
  const double inv_sqrt2 = 1.00/std::sqrt(2.00);

  for(arma::uword constraint=0L;
      constraint < nconstraints; ++constraint) {

    arma::uword column = constraint_pairs(0L, constraint);
    arma::uword row = constraint_pairs(1L, constraint);

    if(row == column) {
      constraint_basis(column+column*q, constraint) = 1.00;
    } else {
      constraint_basis(column+row*q, constraint) = inv_sqrt2;
      constraint_basis(row+column*q, constraint) = inv_sqrt2;
    }

  }

  return constraint_basis;

}

arma::mat poblq_lyap(const arma::mat& psi,
                     const arma::mat& Q,
                     const arma::mat& constraint_basis) {

  const arma::uword q = psi.n_cols;
  const arma::uword nconstraints = constraint_basis.n_cols;

  if(psi.n_rows != q || Q.n_rows != q || Q.n_cols != q ||
     constraint_basis.n_rows != q*q) {
    Rcpp::stop("The matrices supplied to the poblq Lyapunov equation have incompatible dimensions.");
  }

  if(nconstraints == 0L) {
    return arma::mat(q, q, arma::fill::zeros);
  }

  if(!psi.is_finite() || !Q.is_finite() ||
     !constraint_basis.is_finite()) {
    Rcpp::stop("The poblq Lyapunov equation contains non-finite values.");
  }

  arma::mat system(nconstraints, nconstraints,
                   arma::fill::zeros);

  for(arma::uword i=0L; i < nconstraints; ++i) {

    arma::mat B = arma::reshape(constraint_basis.col(i), q, q);
    arma::mat LB = psi*B+B*psi;
    system.col(i) =
      constraint_basis.t()*arma::vectorise(LB);

  }

  system = 0.50*(system+system.t());
  arma::mat Qsym = 0.50*(Q+Q.t());
  arma::vec rhs =
    constraint_basis.t()*arma::vectorise(Qsym);

  arma::mat R;
  bool decomposed = arma::chol(R, system);

  if(!decomposed || !R.is_finite()) {
    Rcpp::stop("The poblq projection is not uniquely defined because the constraint derivatives are linearly dependent.");
  }

  arma::vec intermediate;
  arma::vec coefficients;
  bool solved_lower = arma::solve(intermediate,
                                  arma::trimatl(R.t()), rhs);
  bool solved_upper = arma::solve(coefficients,
                                  arma::trimatu(R), intermediate);

  if(!solved_lower || !solved_upper ||
     !coefficients.is_finite()) {
    Rcpp::stop("The projected poblq Lyapunov equation could not be solved.");
  }

  arma::mat A = arma::reshape(
    constraint_basis*coefficients, q, q
  );
  A = 0.50*(A+A.t());

  arma::vec residual = constraint_basis.t()*arma::vectorise(
    psi*A+A*psi-Qsym
  );
  double residual_scale = 1.00+arma::norm(rhs, "inf");

  if(!residual.is_finite() ||
     arma::norm(residual, "inf") > 1e-08*residual_scale) {
    Rcpp::stop("The projected poblq Lyapunov equation has an excessive numerical residual.");
  }

  return A;

}

arma::mat poblq_retract(arma::mat X,
                        const arma::mat& constraints) {

  if(X.n_cols != constraints.n_cols ||
     constraints.n_rows != constraints.n_cols) {
    Rcpp::stop("The matrix and poblq constraints have incompatible dimensions.");
  }

  if(X.n_rows == 0L || X.n_cols == 0L || !X.is_finite()) {
    Rcpp::stop("The poblq manifold requires a non-empty finite matrix.");
  }

  const arma::uword q = X.n_cols;
  const double tolerance = std::sqrt(arma::datum::eps);

  for(arma::uword i=0L; i < q; ++i) {

    if(i > 0L) {

      arma::vec upper_column = constraints.col(i).head(i);
      arma::uvec zeros = arma::find(upper_column == 0.00);

      if(zeros.n_elem > 0L) {

        arma::mat Q = arma::orth(X.cols(zeros));

        if(!Q.is_finite()) {
          Rcpp::stop("The selective Gram-Schmidt retraction could not be computed.");
        }

        if(Q.n_cols > 0L) {
          X.col(i) -= Q*(Q.t()*X.col(i));
        }

      }

    }

    if(constraints(i, i) == 0.00) {

      double column_norm = arma::norm(X.col(i), 2);

      if(!std::isfinite(column_norm) ||
         column_norm <= tolerance) {
        Rcpp::stop("The selective Gram-Schmidt retraction produced a zero column that must have unit norm.");
      }

      X.col(i) /= column_norm;

    }

  }

  if(!X.is_finite()) {
    Rcpp::stop("The poblq retraction produced non-finite values.");
  }

  return X;

}

// Partially Oblique manifold:

class poblq:public manifolds {

public:

  arma::uvec indices, oblq_indices;
  arma::umat constraint_pairs;
  std::size_t p, q;
  arma::mat X, dX, A, psi, constraints, constraint_basis, g, dg;

  arma::mat normal_coefficients(const arma::mat& G) {

    arma::mat Q = X.t()*G+G.t()*X;

    return poblq_lyap(psi, Q, constraint_basis);

  }

  arma::mat project(const arma::mat& G) {

    arma::mat coefficients = normal_coefficients(G);

    return G-X*coefficients;

  }

  void param(arguments_optim& x) {

    if(p == 0L || q == 0L) {
      Rcpp::stop("The poblq manifold requires positive p and q dimensions.");
    }

    if(indices.n_elem != p*q) {
      Rcpp::stop("The poblq manifold requires p*q parameters.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

    if(!X.is_finite()) {
      Rcpp::stop("The poblq parameters cannot contain non-finite values.");
    }

    psi = X.t()*X;

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    A = normal_coefficients(g);

    x.rg.elem(indices) =
      arma::vectorise(g-X*A);

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    dg = arma::reshape(x.dg.elem(indices), p, q);
    dX = arma::reshape(x.dparameters.elem(indices), p, q);
    A = normal_coefficients(g);

    arma::mat h = dg-dX*A;

    x.dH.elem(indices) =
      arma::vectorise(project(h));

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) =
      arma::vectorise(poblq_retract(X, constraints));

  }

  void tangent_basis(arguments_optim& x) {

    X = arma::reshape(x.parameters(indices), p, q);

    const arma::uword nlocal = p*q;
    const arma::uword nconstraints =
      constraint_basis.n_cols;

    if(nconstraints > nlocal) {
      Rcpp::stop("The poblq manifold has more constraints than parameters.");
    }

    arma::mat T_local;

    if(nconstraints == 0L) {

      T_local = arma::eye(nlocal, nlocal);

    } else {

      arma::mat normal_basis(nlocal, nconstraints,
                             arma::fill::zeros);

      for(arma::uword i=0L; i < nconstraints; ++i) {

        arma::mat B = arma::reshape(
          constraint_basis.col(i), q, q
        );
        normal_basis.col(i) =
          arma::vectorise(X*B);

      }

      bool computed = arma::null(T_local,
                                 normal_basis.t());
      const arma::uword expected_dimension =
        nlocal-nconstraints;

      if(!computed || !T_local.is_finite() ||
         T_local.n_cols != expected_dimension) {
        Rcpp::stop("The poblq tangent-space basis is not defined because the constraint derivatives are linearly dependent.");
      }

    }

    T.zeros(x.nparam, T_local.n_cols);

    for(arma::uword column=0L;
        column < T_local.n_cols; ++column) {
      for(arma::uword row=0L; row < indices.n_elem; ++row) {
        T(indices[row], column) = T_local(row, column);
      }
    }

  }

  void dconstraints(arguments_optim& x) {

    const arma::uword nconstraints =
      constraint_pairs.n_cols;

    if(nconstraints == 0L) {
      return;
    }

    X = arma::reshape(x.parameters(indices), p, q);

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    arma::sp_mat first_derivatives(
      x.transparameters.n_elem,
      nconstraints
    );
    std::vector<arma::sp_mat> second_derivatives;
    second_derivatives.reserve(nconstraints);

    for(arma::uword constraint=0L;
        constraint < nconstraints; ++constraint) {

      arma::uword column =
        constraint_pairs(0L, constraint);
      arma::uword row =
        constraint_pairs(1L, constraint);

      arma::sp_mat second_derivative(
        x.transparameters.n_elem,
        x.transparameters.n_elem
      );

      if(row == column) {

        for(arma::uword i=0L; i < p; ++i) {

          arma::uword local_index = i+column*p;
          arma::uword trans_index =
            trans_indices[local_index];

          first_derivatives(trans_index,
                            constraint) =
            2.00*X(i, column);
          second_derivative(trans_index,
                            trans_index) = 2.00;

        }

      } else {

        for(arma::uword i=0L; i < p; ++i) {

          arma::uword column_index = i+column*p;
          arma::uword row_index = i+row*p;
          arma::uword trans_column =
            trans_indices[column_index];
          arma::uword trans_row =
            trans_indices[row_index];

          first_derivatives(trans_column,
                            constraint) =
            X(i, row);
          first_derivatives(trans_row,
                            constraint) =
            X(i, column);

          second_derivative(trans_column,
                            trans_row) = 1.00;
          second_derivative(trans_row,
                            trans_column) = 1.00;

        }

      }

      second_derivatives.push_back(
        second_derivative
      );

    }

    append_constraint_derivatives(
      x,
      first_derivatives,
      second_derivatives
    );

  }

  void outcomes(arguments_optim& x) {

    tangent_basis(x);

    matrices.resize(1);
    matrices[0] = T;
    names_matrices.resize(1);
    names_matrices[0] = "tangent_basis";

  }

};

poblq* choose_poblq(Rcpp::List manifold_setup) {

  poblq* mymanifold = new poblq();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];
  std::size_t q = manifold_setup["q"];
  arma::mat constraints = manifold_setup["constraints"];

  constraints = check_poblq_constraints(
    constraints, static_cast<arma::uword>(q)
  );
  arma::umat constraint_pairs =
    poblq_constraint_pairs(constraints);
  arma::mat constraint_basis =
    poblq_constraint_basis(
      constraint_pairs,
      static_cast<arma::uword>(q)
    );
  arma::uvec oblq_indices =
    arma::find(constraints == 1.00);

  mymanifold->indices = indices;
  mymanifold->constraints = constraints;
  mymanifold->constraint_pairs = constraint_pairs;
  mymanifold->constraint_basis = constraint_basis;
  mymanifold->p = p;
  mymanifold->q = q;
  mymanifold->oblq_indices = oblq_indices;

  return mymanifold;

}

arma::mat poblq(arma::mat X, arma::mat constraints) {

  constraints = check_poblq_constraints(
    constraints, X.n_cols
  );
  X = poblq_retract(X, constraints);

  return X;

}

arma::mat rpoblq(int p, int q, arma::mat constraints) {

  if(p < 1L || q < 1L) {
    Rcpp::stop("rpoblq requires positive p and q dimensions.");
  }

  constraints = check_poblq_constraints(
    constraints, static_cast<arma::uword>(q)
  );

  arma::mat X(p, q, arma::fill::randn);
  X = poblq_retract(X, constraints);

  return X;

}
