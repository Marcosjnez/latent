/*
 * Author: Marcos Jimenez
 * email: m.j.jimenezhenriquez@vu.nl
 * Modification date: 26/08/2026
 */

// Orthogonal-column manifold with unrestricted column norms:

class orthog:public manifolds {

public:

  arma::uvec indices;
  std::size_t p, q;
  arma::mat X, dX, g, dg, A;

  arma::mat normal_coefficients(const arma::mat& X,
                                const arma::mat& G) {

    arma::vec norms2 = arma::sum(X % X, 0).t();

    if(!norms2.is_finite() ||
       arma::any(norms2 <= arma::datum::eps)) {
      Rcpp::stop("The orthog manifold requires nonzero finite columns.");
    }

    arma::mat B = X.t()*G+G.t()*X;
    arma::mat coefficients(q, q, arma::fill::zeros);

    for(arma::uword i=0L; i+1L < q; ++i) {
      for(arma::uword j=i+1L; j < q; ++j) {

        double value = B(i, j)/(norms2[i]+norms2[j]);
        coefficients(i, j) = value;
        coefficients(j, i) = value;

      }
    }

    return coefficients;

  }

  arma::mat project(const arma::mat& X,
                    const arma::mat& G) {

    arma::mat coefficients =
      normal_coefficients(X, G);

    return G-X*coefficients;

  }

  arma::mat retract(const arma::mat& X) {

    if(X.n_rows < X.n_cols) {
      Rcpp::stop("The orthog manifold requires p >= q.");
    }

    arma::vec norms =
      arma::sqrt(arma::sum(X % X, 0)).t();

    if(!norms.is_finite() ||
       arma::any(norms <= arma::datum::eps)) {
      Rcpp::stop("The orthog manifold requires nonzero finite columns.");
    }

    arma::mat Q, R;
    arma::qr_econ(Q, R, X);

    arma::rowvec signs = arma::sign(R.diag()).t();

    for(arma::uword j=0L; j < signs.n_elem; ++j) {
      if(signs[j] == 0.00) {
        signs[j] = 1.00;
      }
    }

    Q.each_row() %= signs;

    return Q*arma::diagmat(norms);

  }

  void param(arguments_optim& x) {

    if(p < q) {
      Rcpp::stop("The orthog manifold requires p >= q.");
    }

    if(indices.n_elem != p*q) {
      Rcpp::stop("The orthog manifold requires p*q parameters.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

  }

  void proj(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    A = normal_coefficients(X, g);

    x.rg.elem(indices) =
      arma::vectorise(g-X*A);

  }

  void hess(arguments_optim& x) {

    g = arma::reshape(x.g.elem(indices), p, q);
    dg = arma::reshape(x.dg.elem(indices), p, q);
    dX = arma::reshape(x.dparameters.elem(indices), p, q);
    A = normal_coefficients(X, g);

    arma::mat h = dg-dX*A;

    x.dH.elem(indices) =
      arma::vectorise(project(X, h));

  }

  void retr(arguments_optim& x) {

    x.parameters(indices) =
      arma::vectorise(retract(X));

  }

  void tangent_basis(arguments_optim& x) {

    if(p < q) {
      Rcpp::stop("The orthog manifold requires p >= q.");
    }

    X = arma::reshape(x.parameters(indices), p, q);

    arma::vec norms =
      arma::sqrt(arma::sum(X % X, 0)).t();

    if(!norms.is_finite() ||
       arma::any(norms <= arma::datum::eps)) {
      Rcpp::stop("The orthog manifold requires nonzero finite columns.");
    }

    arma::mat Q = X*arma::diagmat(1.00/norms);
    const arma::uword nconstraints =
      q*(q-1L)/2L;
    const arma::uword dimension =
      p*q-nconstraints;

    T.zeros(x.nparam, dimension);
    arma::uword basis_column = 0L;

    // Free changes in the column norms.
    for(arma::uword j=0L; j < q; ++j) {

      arma::mat direction(p, q, arma::fill::zeros);
      direction.col(j) = Q.col(j);
      arma::vec direction_vec =
        arma::vectorise(direction);

      for(arma::uword k=0L; k < indices.n_elem; ++k) {
        T(indices[k], basis_column) = direction_vec[k];
      }

      basis_column++;

    }

    // Pairwise rotations that preserve orthogonality to first order.
    for(arma::uword i=0L; i+1L < q; ++i) {
      for(arma::uword j=i+1L; j < q; ++j) {

        double normalizer =
          std::sqrt(norms[i]*norms[i]+norms[j]*norms[j]);
        arma::mat direction(p, q, arma::fill::zeros);
        direction.col(i) =
          norms[i]/normalizer*Q.col(j);
        direction.col(j) =
          -norms[j]/normalizer*Q.col(i);
        arma::vec direction_vec =
          arma::vectorise(direction);

        for(arma::uword k=0L; k < indices.n_elem; ++k) {
          T(indices[k], basis_column) = direction_vec[k];
        }

        basis_column++;

      }
    }

    // Changes orthogonal to the current column space.
    if(p > q) {

      arma::mat Q_full, R;
      arma::qr(Q_full, R, Q);
      arma::mat Q_perp = Q_full.cols(q, p-1L);

      for(arma::uword j=0L; j < q; ++j) {
        for(arma::uword i=0L; i < p-q; ++i) {

          arma::mat direction(p, q, arma::fill::zeros);
          direction.col(j) = Q_perp.col(i);
          arma::vec direction_vec =
            arma::vectorise(direction);

          for(arma::uword k=0L; k < indices.n_elem; ++k) {
            T(indices[k], basis_column) = direction_vec[k];
          }

          basis_column++;

        }
      }

    }

    if(basis_column != dimension) {
      Rcpp::stop("The orthog tangent-space basis has an invalid dimension.");
    }

  }

  void dconstraints(arguments_optim& x) {

    if(q < 2L) {
      return;
    }

    X = arma::reshape(x.parameters(indices), p, q);

    const arma::uword nconstraints =
      q*(q-1L)/2L;

    arma::uvec trans_indices =
      x.transparam2param.elem(indices);

    arma::sp_mat first_derivatives(
      x.transparameters.n_elem,
      nconstraints
    );
    std::vector<arma::sp_mat> second_derivatives;
    second_derivatives.reserve(nconstraints);

    arma::uword constraint = 0L;

    for(arma::uword column=0L; column+1L < q; ++column) {
      for(arma::uword row=column+1L; row < q; ++row) {

        arma::sp_mat second_derivative(
          x.transparameters.n_elem,
          x.transparameters.n_elem
        );

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

        second_derivatives.push_back(
          second_derivative
        );
        constraint++;

      }
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

orthog* choose_orthog(Rcpp::List manifold_setup) {

  orthog* mymanifold = new orthog();

  arma::uvec indices = manifold_setup["indices"];
  std::size_t p = manifold_setup["p"];
  std::size_t q = manifold_setup["q"];

  mymanifold->indices = indices;
  mymanifold->p = p;
  mymanifold->q = q;

  return mymanifold;

}
