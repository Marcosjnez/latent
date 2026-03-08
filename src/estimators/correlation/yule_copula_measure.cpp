// yule_copula_measure.cpp
#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------
// Utilities
// -----------------------------

// Safe kron (Armadillo provides kron)
static inline arma::vec kron_vec(const arma::vec& a, const arma::vec& b) {
    // return vec(kron(a,b)) but ensure column vectors
    return arma::vectorise(arma::kron(a, b)); // column-major
}

// SVD-based pseudoinverse wrapper (Armadillo pinv)
static inline arma::mat pinv_mat(const arma::mat& X, double tol = 1e-12) {
    return arma::pinv(X, tol);
}

// Position mapping: 0-based (i,j) -> pos = j*R + i
static inline arma::uword pos_idx(unsigned int i, unsigned int j, unsigned int R) {
    return (arma::uword)j * (arma::uword)R + (arma::uword)i;
}

// -----------------------------
// IPFP (Armadaillo) -- counts scale, uniform rectangular targets
// targets: floor(total / R) repeated R times ; floor(total / S) repeated S times
// returns a matrix that sums to 1 (copula pmf)
// -----------------------------
// [[Rcpp::export]]
arma::mat copula_pmf_armadillo(const arma::mat& H,
    int max_iter = 2000,
    double tol = 1e-12) {
    const unsigned int R = H.n_rows;
    const unsigned int S = H.n_cols;

    // Work in counts scale (clone H)
    arma::mat K = H;

    double total = arma::accu(H);
    if (total <= 0.0) {
        // empty table -> return uniform copula pmf if possible
        arma::mat out = arma::ones<arma::mat>(R, S);
        out /= arma::accu(out);
        return out;
    }

    double restriction_x = std::floor(total / (double)R);
    double restriction_y = std::floor(total / (double)S);

    arma::vec target_rows(R);
    target_rows.fill(restriction_x);
    arma::vec target_cols(S);
    target_cols.fill(restriction_y);

    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        // Adjust rows
        for (unsigned int i = 0; i < R; ++i) {
            double sum_row = arma::accu(K.row(i));
            if (sum_row <= 0.0) continue;
            double factor = target_rows(i) / sum_row;
            K.row(i) *= factor;
            double change = std::fabs(sum_row * factor - sum_row);
            if (change > max_diff) max_diff = change;
        }

        // Adjust cols
        for (unsigned int j = 0; j < S; ++j) {
            double sum_col = arma::accu(K.col(j));
            if (sum_col <= 0.0) continue;
            double factor = target_cols(j) / sum_col;
            K.col(j) *= factor;
            double change = std::fabs(sum_col * factor - sum_col);
            if (change > max_diff) max_diff = change;
        }

        if (max_diff < tol) break;
    }

    double tsum = arma::accu(K);
    if (tsum > 0.0) {
        K /= tsum;
    }
    else {
        // fallback
        K = arma::ones<arma::mat>(R, S);
        K /= arma::accu(K);
    }

    return K;
}

// -----------------------------
// Build A_Sp (exact R-style):
// If no structural zeros -> A_Sp = kron(CS, CR) where CR = I_R - 1/R J_R
// Else -> build support Sp, Z = [1 | row indicators | col indicators] (m x (1+R+S))
// compute P_Z = Z (Z'Z)^+ Z' and A_support = I_m - P_Z
// embed A_support into RS x RS at the support indices
// -----------------------------
static inline arma::mat build_A_Sp_from_H(const arma::mat& H) {
    unsigned int R = H.n_rows;
    unsigned int S = H.n_cols;

    arma::uvec support_idx; // indices in 0..R*S-1 that belong to support
    arma::uvec all_idx(R * S);
    for (unsigned int k = 0; k < R * S; ++k) all_idx(k) = k;

    // Detect structural zeros: structural zeros are those cells that are impossible.
    // Here we consider H == 0 as structural if any zero exists. (User may refine.)
    bool has_structural = false;
    for (unsigned int i = 0; i < R; ++i)
        for (unsigned int j = 0; j < S; ++j)
            if (H(i, j) == 0.0) { has_structural = true; break; }
    // If all cells positive, no structural zeros
    if (!has_structural) {
        // Rectangular canonical basis
        arma::mat J_R = arma::ones<arma::mat>(R, R);
        arma::mat J_S = arma::ones<arma::mat>(S, S);
        arma::mat CR = arma::eye<arma::mat>(R, R) - (1.0 / (double)R) * J_R;
        arma::mat CS = arma::eye<arma::mat>(S, S) - (1.0 / (double)S) * J_S;
        arma::mat A_R = arma::kron(CS, CR); // size (R*S) x (R*S)
        return A_R;
    }
    else {
        // Build Sp as positions with H > 0 (strictly positive)
        std::vector<arma::uword> idx_vec;
        for (unsigned int j = 0; j < S; ++j) {
            for (unsigned int i = 0; i < R; ++i) {
                if (H(i, j) > 0.0) {
                    idx_vec.push_back(pos_idx(i, j, R));
                }
            }
        }
        unsigned int m = idx_vec.size();
        if (m == 0) {
            // fallback: no support -> return rectangular basis
            arma::mat J_R = arma::ones<arma::mat>(R, R);
            arma::mat J_S = arma::ones<arma::mat>(S, S);
            arma::mat CR = arma::eye<arma::mat>(R, R) - (1.0 / (double)R) * J_R;
            arma::mat CS = arma::eye<arma::mat>(S, S) - (1.0 / (double)S) * J_S;
            arma::mat A_R = arma::kron(CS, CR); // size (R*S) x (R*S)
            return A_R;
        }

        // Construct matrix Z (m x (1+R+S))
        arma::mat Z = arma::zeros<arma::mat>(m, 1 + R + S);
        // intercept
        Z.col(0).ones();

        // For each support position, find its (i,j)
        for (unsigned int t = 0; t < m; ++t) {
            arma::uword pos = idx_vec[t];
            // recover i,j: pos = j*R + i
            arma::uword j = pos / R;
            arma::uword i = pos % R;
            // set row indicator (1-based positions in R columns)
            Z(t, 1 + i) = 1.0;
            // set col indicator (after R columns)
            Z(t, 1 + R + j) = 1.0;
        }

        // Compute projection onto columns of Z: PZ = Z (Z'Z)^+ Z'
        arma::mat ZtZ = Z.t() * Z;
        arma::mat ZtZ_pinv = pinv_mat(ZtZ);
        arma::mat PZ = Z * ZtZ_pinv * Z.t();

        arma::mat A_support = arma::eye<arma::mat>(m, m) - PZ;

        // Embed into full RS x RS matrix
        arma::mat A_full = arma::zeros<arma::mat>(R * S, R * S);
        for (unsigned int a = 0; a < m; ++a) {
            for (unsigned int b = 0; b < m; ++b) {
                arma::uword row = idx_vec[a];
                arma::uword col = idx_vec[b];
                A_full(row, col) = A_support(a, b);
            }
        }
        return A_full;
    }
}

// -----------------------------
// Core: compute upsilon + variance + se from contingency H
// Returns list(upsilon, variance, se, structural_zeros)
// -----------------------------
// [[Rcpp::export]]
Rcpp::List upsilon_full_H_rcpp(const arma::mat& H_in,
    int ipfp_maxit = 2000,
    double ipfp_tol = 1e-12,
    double pinv_tol = 1e-12) {
    // H_in: contingency table (counts), numeric
    arma::mat H = H_in;
    unsigned int R = H.n_rows;
    unsigned int S = H.n_cols;
    double N = arma::accu(H);

    if (N <= 0.0) {
        Rcpp::stop("Empty contingency table (total count <= 0).");
    }

    // Empirical joint pmf
    arma::mat P = H / N;

    // Copula pmf via Armadillo IPFP (uniform rectangular targets)
    arma::mat Cop = copula_pmf_armadillo(H, ipfp_maxit, ipfp_tol);

    // Detect structural zeros
    bool has_structural = false;
    for (unsigned int i = 0; i < R && !has_structural; ++i)
        for (unsigned int j = 0; j < S; ++j)
            if (H(i, j) == 0.0) { has_structural = true; break; }

    // Compute upsilon (use indices 1..R and 1..S as in your R function)
    arma::vec r1 = arma::regspace<arma::vec>(1.0, (double)R);
    arma::vec r2 = arma::regspace<arma::vec>(1.0, (double)S);

    double Csum = 0.0;
    for (unsigned int i = 0; i < R; ++i)
        for (unsigned int j = 0; j < S; ++j)
            Csum += (double)(i + 1) * (double)(j + 1) * Cop(i, j);

    double kappa = 12.0 / std::sqrt((double)(R - 1) * (double)(R + 1) * (double)(S - 1) * (double)(S + 1));
    double rest = -3.0 * std::sqrt(((double)(R + 1) * (double)(S + 1)) / ((double)(R - 1) * (double)(S - 1)));
    double upsilon = kappa * Csum + rest;

    // Build A_Sp (RS x RS)
    arma::mat A_Sp = build_A_Sp_from_H(H);

    // Vectorize Cop and P (column-major)
    arma::vec g = arma::vectorise(Cop); // length R*S
    arma::vec p = arma::vectorise(P);

    // Build diag^-1 matrices (RS x RS)
    // Guard against zero entries in g or p by replacing zero with small positive to avoid Inf
    double eps = 1e-20;
    arma::vec g_safe = g;
    arma::vec p_safe = p;
    for (arma::uword i = 0; i < g_safe.n_elem; ++i) if (g_safe(i) <= 0.0) g_safe(i) = eps;
    for (arma::uword i = 0; i < p_safe.n_elem; ++i) if (p_safe(i) <= 0.0) p_safe(i) = eps;

    arma::mat Dg_inv = arma::diagmat(1.0 / g_safe);
    arma::mat Dp_inv = arma::diagmat(1.0 / p_safe);

    // Compute M1, M2
    arma::mat M1 = A_Sp.t() * Dg_inv * A_Sp;
    arma::mat M2 = A_Sp.t() * Dp_inv * A_Sp;

    // Moore-Penrose inverse
    arma::mat M1_inv = pinv_mat(M1, pinv_tol);

    // Sigma_gamma (eq. 4.11)
    arma::mat Sigma_gamma = A_Sp * M1_inv * M2 * M1_inv * A_Sp.t();

    // k_prod = kron(r2, r1) (vector of length R*S)
    arma::vec k_prod = kron_vec(r2, r1);

    // sigma^2 (eq. 5.5)
    double sigma2 = kappa * kappa * as_scalar(k_prod.t() * Sigma_gamma * k_prod);

    // standard error scaled by sqrt(N)
    double se = 0.0;
    if (sigma2 >= 0.0) se = std::sqrt(sigma2) / std::sqrt(N);
    else se = std::numeric_limits<double>::quiet_NaN();

    return Rcpp::List::create(
        Rcpp::Named("upsilon") = upsilon,
        Rcpp::Named("variance") = sigma2,
        Rcpp::Named("se") = se,
        Rcpp::Named("structural_zeros") = has_structural
    );
}

// -----------------------------
// Overload: compute from two integer-coded vectors (R factors -> as.integer)
// xi, yi are IntegerVector (R-side) - may contain arbitrary integer codes
// -----------------------------
// [[Rcpp::export]]
Rcpp::List upsilon_full_xy_rcpp(const IntegerVector& x_in,
    const IntegerVector& y_in,
    int ipfp_maxit = 2000,
    double ipfp_tol = 1e-12,
    double pinv_tol = 1e-12) {
    // convert to integer vectors
    int N = x_in.size();
    if ((int)y_in.size() != N) Rcpp::stop("x and y must have same length.");

    // Map unique sorted categories to 0..R-1 and 0..S-1
    IntegerVector ux = sort_unique(x_in);
    IntegerVector uy = sort_unique(y_in);
    int R = ux.size();
    int S = uy.size();

    std::unordered_map<int, int> mapx;
    for (int i = 0; i < R; ++i) mapx[ux[i]] = i;
    std::unordered_map<int, int> mapy;
    for (int j = 0; j < S; ++j) mapy[uy[j]] = j;

    arma::mat H = arma::zeros<arma::mat>((unsigned int)R, (unsigned int)S);
    for (int k = 0; k < N; ++k) {
        int xi = x_in[k];
        int yi = y_in[k];
        int i = mapx[xi];
        int j = mapy[yi];
        H(i, j) += 1.0;
    }

    return upsilon_full_H_rcpp(H, ipfp_maxit, ipfp_tol, pinv_tol);
}

// -----------------------------
// Pairwise matrix: accepts R data.frame (handles factor columns)
// Returns list(cor=matrix, se=matrix) similar to R version
// -----------------------------
// [[Rcpp::export]]
Rcpp::List yule_cor_full_rcpp(const Rcpp::DataFrame& X_df,
    int ipfp_maxit = 2000,
    double ipfp_tol = 1e-12,
    double pinv_tol = 1e-12) {
    int P = X_df.size();
    int N = X_df.nrows();

    arma::mat cor_mat = arma::eye<arma::mat>(P, P);
    arma::mat se_mat = arma::eye<arma::mat>(P, P);

    // Pre-convert each column to integer codes (as in R: as.integer(factor))
    std::vector< IntegerVector > cols(P);
    for (int a = 0; a < P; ++a) {
        SEXP col = X_df[a];
        // If factor or integer, as<IntegerVector> will work
        // If numeric, we convert to integer via as<IntegerVector> (may truncate)
        cols[a] = as<IntegerVector>(col);
    }

    for (int a = 0; a < P; ++a) {
        for (int b = a + 1; b < P; ++b) {
            Rcpp::List out = upsilon_full_xy_rcpp(cols[a], cols[b], ipfp_maxit, ipfp_tol, pinv_tol);
            double u = as<double>(out["upsilon"]);
            double se = as<double>(out["se"]);
            cor_mat(a, b) = u;
            cor_mat(b, a) = u;
            se_mat(a, b) = se;
            se_mat(b, a) = se;
        }
    }

    // Convert to R objects
    return Rcpp::List::create(
        Rcpp::Named("cor") = cor_mat,
        Rcpp::Named("se") = se_mat
    );
}

// -----------------------------
// Recover observed joint pmf from copula pmf + marginals
// -----------------------------
// [[Rcpp::export]]
arma::mat observed_pmf_from_copula_rcpp(const arma::mat& Cop,      // copula pmf (R x S)
    const arma::vec& pr,       // row marginals (length R)
    const arma::vec& pc,       // column marginals (length S)
    int max_iter = 2000,
    double tol = 1e-12) {
    unsigned int R = Cop.n_rows;
    unsigned int S = Cop.n_cols;

    if (pr.n_elem != R)
        Rcpp::stop("Length of pr must match number of rows of Cop.");
    if (pc.n_elem != S)
        Rcpp::stop("Length of pc must match number of columns of Cop.");

    if (std::abs(arma::accu(pr) - 1.0) > 1e-10)
        Rcpp::stop("Row marginals must sum to 1.");
    if (std::abs(arma::accu(pc) - 1.0) > 1e-10)
        Rcpp::stop("Column marginals must sum to 1.");

    arma::mat P = Cop;

    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        // --- Row scaling ---
        for (unsigned int i = 0; i < R; ++i) {
            double rs = arma::accu(P.row(i));
            if (rs > 0.0) {
                double f = pr(i) / rs;
                P.row(i) *= f;
                max_diff = std::max(max_diff, std::abs(rs - pr(i)));
            }
        }

        // --- Column scaling ---
        for (unsigned int j = 0; j < S; ++j) {
            double cs = arma::accu(P.col(j));
            if (cs > 0.0) {
                double f = pc(j) / cs;
                P.col(j) *= f;
                max_diff = std::max(max_diff, std::abs(cs - pc(j)));
            }
        }

        if (max_diff < tol)
            break;
    }

    // Numerical safeguard
    double total = arma::accu(P);
    if (total > 0.0)
        P /= total;

    return P;
}

// -----------------------------
// Simulate contingency table from copula + marginals
// Asymptotically: H / N -> observed_pmf_from_copula
// -----------------------------
// [[Rcpp::export]]
arma::mat simulate_table_from_copula_rcpp(const arma::mat& Cop,   // copula pmf (R x S)
    const arma::vec& pr,    // row marginals
    const arma::vec& pc,    // column marginals
    int N                   // sample size
) {
    if (N <= 0)
        Rcpp::stop("N must be positive.");

    // Step 1: construct target joint pmf
    arma::mat P =
        observed_pmf_from_copula_rcpp(Cop, pr, pc);

    unsigned int R = P.n_rows;
    unsigned int S = P.n_cols;
    unsigned int RS = R * S;

    // Step 2: vectorize probabilities (column-major)
    arma::vec prob = arma::vectorise(P);

    // Numerical safeguard
    double psum = arma::accu(prob);
    if (std::abs(psum - 1.0) > 1e-10)
        prob /= psum;

    // Step 3: multinomial draw
    arma::vec counts(RS, arma::fill::zeros);

    // Draw N categorical samples
    for (int n = 0; n < N; ++n) {
        double u = R::runif(0.0, 1.0);
        double csum = 0.0;

        for (arma::uword k = 0; k < RS; ++k) {
            csum += prob(k);
            if (u <= csum) {
                counts(k) += 1.0;
                break;
            }
        }
    }

    // Step 4: reshape back to R x S
    arma::mat H = arma::reshape(counts, R, S);

    return H;
}
