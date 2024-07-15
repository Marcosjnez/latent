#' @title
#' Generate random orthogonal matrices.
#' @description
#'
#' Generate random orthogonal matrices from a standard normal distribution. First, a matrix of random standard normal variables is simulated and then, the Q factor from the QR decomposition is returned.
#'
#' @usage
#'
#' random_orth(p, q)
#'
#' @param p Number of rows.
#' @param q Number of columns. Should not be greater than p.
#'
#' @return An orthogonal matrix with normally distributed data.
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
random_orth <- function(p, q) {
  .Call(`_latent_random_orth`, p, q)
}

#' @title
#' Generate random oblique matrices.
#' @description
#'
#' Generate random oblique matrices from a standard normal distribution.
#'
#' @usage
#'
#' random_oblq(p, q)
#'
#' @param p Number of rows.
#' @param q Number of columns. Should not be greater than p.
#'
#' @return An oblique matrix with normally distributed data.
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
random_oblq <- function(p, q) {
  .Call(`_latent_random_oblq`, p, q)
}

#' @title
#' Generate a random partially oblique matrix.
#' @description
#'
#' First, a matrix is simulated from a standard normal distribution. Second,the matrix is normalized and the Gram-Schmidt process is performed between the oblique blocks. Finally, the orthogonal blocks correspond to those columns of the Q matrix from the QR decomposition.
#'
#' @usage
#'
#' random_poblq(p, q, oblq_factors)
#'
#' @param p Number of rows.
#' @param q Number of columns. Should not be greater than p.
#' @param oblq_factors A vector with the number of factors for each oblique block. E.g.: c(2, 4) means that there are two blocks of oblique factors: one with 2 factors and another with 4 factors. Everything else is orthogonal.
#'
#' @return A partially oblique matrix.
#'
#' @examples
#'
#' random_poblq(p = 7, q = 7, oblq_factors = c(3, 2))
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
random_poblq <- function(p, q, oblq_factors) {
  .Call(`_latent_random_poblq`, p, q, oblq_factors)
}

#' @title
#' Retraction of a matrix onto the orthogonal manifold.
#' @description
#'
#' Transform a matrix into an orthogonal matrix.
#'
#' @usage
#'
#' retr_orth(X)
#'
#' @param X A matrix.
#'
#' @return An orthogonal matrix.
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
retr_orth <- function(X) {
  .Call(`_latent_retr_orth`, X)
}

#' @title
#' Retraction of a matrix onto the oblique manifold.
#' @description
#'
#' Transform a matrix into an oblique matrix.
#'
#' @usage
#'
#' retr_oblq(X)
#'
#' @param X A matrix.
#'
#' @return An oblique matrix.
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
retr_oblq <- function(X) {
  .Call(`_latent_retr_oblq`, X)
}

#' @title
#' Retraction of a matrix onto the partially oblique manifold.
#' @description
#'
#' Transform a matrix into a partially oblique matrix.
#'
#' @usage
#'
#' retr_poblq(X, oblq_factors, PhiTarget)
#'
#' @param X A matrix.
#' @param oblq_factors A vector with the number of factors for each oblique block. E.g.: c(2, 4) means that there are two blocks of oblique factors: one with 2 factors and another with 4 factors. Everything else is orthogonal.
#'
#' @return A partially oblique matrix.
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @examples
#'
#' X <- replicate(8, rnorm(8))
#' retr_poblq(X, oblq_factors = c(2, 3, 3))
#'
#' @export
retr_poblq <- function(X, oblq_factors = NULL, PhiTarget = NULL) {
  .Call(`_latent_retr_poblq`, X, oblq_factors, PhiTarget)
}

#' @title
#' Schmid-Leiman Transformation.
#' @description
#'
#' Schmid-Leiman transformation into a bi-factor pattern with one or multiple general factors.
#'
#' @usage
#'
#' sl(X, n_generals, n_groups, cor = "pearson", estimator = "uls", nobs = NULL,
#' first_efa = NULL, second_efa = NULL, cores = 1L)
#'
#' @param X Raw data matrix or correlation matrix.
#' @param n_generals Number of general factors.
#' @param n_groups Number of group factors.
#' @param cor Correlation method. Available correlations: c("pearson", "poly"). Defaults to "pearson".
#' @param estimator EFA fitting estimator: "ml" (maximum likelihood for multivariate normal variables), "uls" (minimum residuals), "pa" (principal axis) and "minrank" (minimum rank). Defaults to "uls".
#' @param missing The way to handle missing data. Options: c("impute.mean", "impute.median", "complete.cases", "pairwise.complete.cases"). Defaults to "pairwise.complete.cases".
#' @param nobs Sample size. Defaults to NULL.
#' @param first_efa Arguments to pass to \code{efast} in the first-order factor extraction. See \code{efast} for the default arguments.
#' @param second_efa Arguments to pass to \code{efast} in the second-order factor extraction. See \code{efast} for the default arguments.
#' @param cores Number of cores for the polychorics estimation.
#'
#' @details First, a hierarchical factor model is fitted using a second-order factor analysis on the factor correlation obtained from a first-order factor analysis. Then, the item loadings on the general factors are assumed to be the direct effects of the general factors according to such hierarchical model.
#' On the other hand, the item loadings on the group factors become the originally first-order loadings post-multiplied by the diagonal matrix containing the root of the item uniquenesses.
#'
#' Obviously, the first-order factor analysis should be oblique to perform a second exploratory factor analysis.
#'
#' If the second-order solution does not use an orthogonal projection, then the correlation matrix among the general factors for the Schmid-Leiman solution is simply that obtained from such second-order solution.
#'
#' @return
#'
#' \item{loadings}{Loading matrix of the Schmid-Leiman solution.}
#' \item{first_order_solution}{Object of class \code{efast} with the first-order solution.}
#' \item{second_order_solution}{Object of class \code{efast} with the second-order solution.}
#' \item{uniquenesses}{Vector of uniquenesses.}
#' \item{Rhat}{Correlation matrix predicted by the (hierarchical) model.}
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @examples
#'
#' \dontrun{
#' # Simulate data:
#' sim <- sim_factor(n_generals = 2, groups_per_general = 3, items_per_group = 5)
#' lambda <- sim$lambda
#' Target <- ifelse(lambda > 0, 1, 0)
#'
#' # Target rotation for the first-order efa and oblimin for the second-order efa:
#' first <- list(rotation = "target", projection = "oblq", Target = Target[, -c(1:2)])
#' second <- list(rotation = "oblimin", projection = "oblq", gamma = 0)
#'
#' SL <- <- sl(sim$R, n_generals = 2, n_groups = 6, nobs = 100,
#'          first_efa = first, second_efa = second)
#'}
#'
#' @export
sl <- function(X, n_generals, n_groups, cor = "pearson", estimator = "uls", missing = "pairwise.complete.cases", nobs = NULL, first_efa = NULL, second_efa = NULL, cores = 1L) {
  .Call(`_latent_sl`, X, n_generals, n_groups, cor, estimator, missing, nobs, first_efa, second_efa, cores)
}

#' @title
#' Fast rotation algorithm for factor analysis.
#' @description
#'
#' Riemannian Newton Trust-Region algorithm to quickly perform (parallel) rotations with different random starting values.
#'
#' @usage
#'
#' rotate(lambda, rotation = "oblimin", projection = "oblq",
#' gamma = 0, epsilon = 0.01, k = 0, w = 1,
#' Target = NULL, Weight = NULL, PhiTarget = NULL, PhiWeight = NULL,
#' blocks = NULL, block_weights = NULL, oblq_factors = NULL,
#' normalization = "none",
#' rot_control = NULL, random_starts = 1L, cores = 1L)
#'
#' @param lambda Unrotated loading matrix.
#' @param rotation Rotation criterion. Available rotations: "varimax", "cf" (Crawford-Ferguson), "oblimin", "geomin", "target", "xtarget" (extended target) and "none". Defaults to "oblimin".
#' @param projection Projection method. Available projections: "orth" (orthogonal), "oblq" (oblique), "poblq" (partially oblique). Defaults to "oblq".
#' @param gamma \eqn{\gamma} parameter for the oblimin criterion. Defaults to 0 (quartimin).
#' @param epsilon \eqn{\epsilon} parameter for the geomin criterion. Defaults to 0.01.
#' @param k \eqn{k} parameter for the Crawford-Ferguson family of rotation criteria. Defaults to 0.
#' @param w \eqn{w} parameter for the extended target criterion ("xtarget"). Defaults to 1.
#' @param Target Target matrix for the loadings. Defaults to NULL.
#' @param Weight Weight matrix for the loadings. Defaults to NULL.
#' @param PhiTarget Target matrix for the factor correlations. Defaults to NULL.
#' @param PhiWeight Weight matrix for the factor correlations. Defaults to NULL.
#' @param blocks Vector with the number of factors for which separately applying the rotation criterion. Defaults to NULL.
#' @param block_weights Vector of weights for each block of factors.
#' @param oblq_factors Vector with the number of factors for each oblique block. E.g.: c(2, 4) means that there are two blocks of oblique factors: one block with 2 factors and another block with 4 factors. Everything else is orthogonal. Defaults to NULL.
#' @param normalization Available normalizations: "kaiser". Defaults to "none".
#' @param rot_control List of control parameters for the rotation algorithm. Defaults to NULL.
#' @param random_starts Number of rotations with different random starting values. The rotation with the smallest cost function value is returned. Defaults to 1L.
#' @param cores Number of cores for parallel execution of random starts. Defaults to 1L.
#'
#' @details
#'
#' If \code{rot_control = NULL}, then \code{list(maxit = 1000, eps = 1e-05)} is passed to \code{rot_control}, where \code{eps} is the absolute tolerance. When the objective function does not make a larger improvement than \code{eps}, the algorithm is assumed to converge.
#' If \code{Target} is provided but not \code{Weight}, then \code{Weight = 1 - Target} by default, which means a partially specified target rotation is performed. The same applies for \code{PhiTarget} and \code{PhiWeight}.
#'
#' @return List of class \code{rotation} with the following components:
#' \item{loadings}{Rotated loading matrix.}
#' \item{Phi}{Correlation matrix among the factors.}
#' \item{T}{Rotation matrix.}
#' \item{f}{Objective value at the minimum.}
#' \item{iterations}{Number of iterations for the rotation algorithm to converge.}
#' \item{convergence}{TRUE if the algorithm converged and FALSE otherwise.}
#' \item{elapsed}{Total amount of time spent for execution (in nanoseconds).}
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' Zhang, G., Hattori, M., Trichtinger, L. A., & Wang, X. (2019). Target rotation with both factor loadings and factor correlations. Psychological Methods, 24(3), 390–402. https://doi.org/10.1037/met0000198
#'
#' @export
rotate <- function(lambda, rotation = "oblimin", projection = "oblq", gamma = 0, epsilon = 0.01, k = 0, w = 1, Target = NULL, Weight = NULL, PhiTarget = NULL, PhiWeight = NULL, blocks = NULL, block_weights = NULL, oblq_factors = NULL, normalization = "none", rot_control = NULL, random_starts = 1L, cores = 1L) {
  .Call(`_latent_rotate`, lambda, rotation, projection, gamma, epsilon, k, w, Target, Weight, PhiTarget, PhiWeight, blocks, block_weights, oblq_factors, normalization, rot_control, random_starts, cores)
}

#' @title
#' Fast exploratory factor analysis.
#' @description
#'
#' Fast exploratory factor analysis.
#'
#' @usage
#'
#' efast(X, nfactors, cor = "pearson", estimator = "uls",
#' rotation = "oblimin", projection = "oblq", nobs = NULL,
#' Target = NULL, Weight = NULL, PhiTarget = NULL, PhiWeight = NULL,
#' blocks = NULL, block_weights = NULL,
#' oblq_factors = NULL, gamma = 0,
#' epsilon = 1e-02, k = 0, w = 1,
#' random_starts = 1L, cores = 1L,
#' init = NULL, efa_control = NULL, rot_control = NULL)
#'
#' @param X Raw data matrix or correlation matrix.
#' @param nfactors Number of common factors to extract.
#' @param cor Correlation method. Available correlations: c("pearson", "poly"). Defaults to "pearson".
#' @param estimator EFA fitting estimator: "ml" (maximum likelihood for multivariate normal variables), "uls" (minimum residuals), "pa" (principal axis) and "minrank" (minimum rank). Defaults to "uls".
#' @param rotation Rotation criterion. Available rotations: "varimax", "cf" (Crawford-Ferguson), "oblimin", "geomin", "target", "xtarget" (extended target) and "none". Defaults to "oblimin".
#' @param projection Projection method. Available projections: "orth" (orthogonal), "oblq" (oblique), "poblq" (partially oblique). Defaults to "oblq".
#' @param missing The way to handle missing data. Options: c("impute.mean", "impute.median", "complete.cases", "pairwise.complete.cases"). Defaults to "pairwise.complete.cases".
#' @param nobs Sample size. Defaults to NULL.
#' @param Target Target matrix for the loadings. Defaults to NULL.
#' @param Weight Weight matrix for the loadings. Defaults to NULL.
#' @param PhiTarget Target matrix for the factor correlations. Defaults to NULL.
#' @param PhiWeight Weight matrix for the factor correlations. Defaults to NULL.
#' @param blocks List of vectors with the indexes of the factors for which separately applying the rotation criterion. Defaults to NULL.
#' @param block_weights Vector of weights for each block of factors.
#' @param oblq_factors Vector with the number of factors for each oblique block. E.g.: c(2, 4) means that there are two blocks of oblique factors: one block with 2 factors and another block with 4 factors. Everything else is orthogonal. Defaults to NULL.
#' @param gamma \eqn{\gamma} parameter for the oblimin criterion. Defaults to 0 (quartimin).
#' @param epsilon \eqn{\epsilon} parameter for the geomin criterion. Defaults to 0.01.
#' @param k \eqn{k} parameter for the Crawford-Ferguson family of rotation criteria. Defaults to 0.
#' @param w \eqn{w} parameter for the extended target criterion ("xtarget"). Defaults to 1L.
#' @param random_starts Number of rotations with different random starting values. The rotation with the smallest cost function value is returned. Defaults to 1L.
#' @param cores Number of cores for parallel execution of random starts. Defaults to 1L.
#' @param init Initial uniquenesses values for exploratory factor analysis estimation. Defaults to NULL.
#' @param efa_control List of control parameters for efa fitting. Defaults to NULL.
#' @param rot_control List of control parameters for the rotation algorithm. Defaults to NULL.
#'
#' @details
#'
#' If \code{efa.control = NULL}, then \code{list(maxit = 1e4)} is passed to \code{efa.control}. If \code{rot_control = NULL}, then \code{list(maxit = 1000, eps = 1e-05)} is passed to \code{rot_control}, where \code{eps} is the absolute tolerance. When the objective function does not make a larger improvement than \code{eps}, the algorithm is assumed to converge.
#'
#' If \code{Target} is provided but not \code{Weight}, then \code{Weight = 1 - Target} by default, which means a partially specified target rotation is performed. The same applies for \code{PhiTarget} and \code{PhiWeight}.
#'
#' If \code{init = NULL}, then the squared multiple correlations of each item with the remaining ones are used as initial values (These are known to be upper bounds).
#'
#' If a Heywood case is encountered, then \code{estimator =} "minrank" is automatically applied to ensure positive uniquenesses.
#'
#' @return List of class \code{efast} with the following components:
#' \item{efa}{List containing the following objects:}
#'
#' \itemize{
#' \item loadings - Unrotated loadings.
#' \item uniquenesses - Vector of uniquenesses.
#' \item Rhat - Correlation matrix predicted by the model.
#' \item residuals - Residual correlation matrix.
#' \item f - Objective value at the minimum.
#' \item Heywood - TRUE if any Heywood case is encountered and FALSE otherwise.
#' \item iterations - Number of iterations for the L-BFGS-B algorithm to converge.
#' \item convergence - TRUE if the L-BFGS-B algorithm converged and FALSE otherwise.
#' \item estimator - Method used to fit the exploratory factor analysis.
#' }
#'
#' \item{rotation}{List of class \code{rotation}. Only if the argument \code{rotation} is not "none". See \code{rotate} for the components.}
#' \item{elapsed}{Total amount spent for execution (in nanoseconds).}
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @examples
#'
#' \dontrun{
#' # Simulate data:
#' sim <- sim_factor(n_generals = 0, groups_per_general = 5, items_per_group = 6)
#' scores <- MASS::mvrnorm(1e3, rep(0, nrow(sim$R)), Sigma = sim$R)
#' s <- cor(scores)
#'
#' # Fit efa:
#' fit <- efast(s, nfactors = 5, estimator = "uls", rotation = "oblimin",
#' projection = "oblq", gamma = 0, random_starts = 10L, cores = 1L)
#'}
#'
#' @export
efast <- function(X, nfactors, cor = "pearson", estimator = "uls", rotation = "oblimin", projection = "oblq", missing = "pairwise.complete.cases", nobs = NULL, Target = NULL, Weight = NULL, PhiTarget = NULL, PhiWeight = NULL, blocks = NULL, block_weights = NULL, oblq_factors = NULL, gamma = 0, epsilon = 1e-02, k = 0, w = 1, random_starts = 1L, cores = 1L, init = NULL, efa_control = NULL, rot_control = NULL) {
  .Call(`_latent_efast`, X, nfactors, cor, estimator, rotation, projection, missing, nobs, Target, Weight, PhiTarget, PhiWeight, blocks, block_weights, oblq_factors, gamma, epsilon, k, w, random_starts, cores, init, efa_control, rot_control)
}

#' @title
#' Get a target from a loading matrix.
#' @description
#'
#' Get a target for the loading matrix using a custom or empirical cut-off.
#'
#' @param loadings A matrix of loadings.
#' @param Phi A correlation matrix among the factors. Defaults to NULL.
#' @param cutoff The cut-off used to create the target matrix. Defaults to 0.
#'
#' @details
#'
#' If \code{cutoff} is not 0, loadings smaller than such a cut-off are fixed to 0. When \code{cutoff = 0}, an empirical cut-off is used for each column of the loading matrix. They are the mean of the one-lagged differences of the sorted squared normalized loadings. Then, the target is determined by fixing to 0 the squared normalized loadings smaller than such cut-offs.
#'
#' @return A target matrix.
#'
#' @references
#'
#' Garcia-Garzon, E., Abad, F. J., & Garrido, L. E. (2019). Improving bi-factor exploratory modeling: Empirical target rotation based on loading differences. Methodology: European Journal of Research Methods for the Behavioral and Social Sciences, 15(2), 45–55. https://doi.org/10.1027/1614-2241/a000163
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @export
get_target <- function(loadings, Phi = NULL, cutoff = 0) {
  .Call(`_latent_get_target`, loadings, Phi, cutoff)
}

#' @title
#' Fit an exploratory bi-factor model with one or multiple general factors.
#' @usage
#'
#' bifactor(X, n_generals, n_groups, method = "GSLiD", cor = "pearson",
#' estimator = "uls", projection = "oblq", nobs = NULL, PhiTarget = NULL,
#' PhiWeight = NULL, blocks = NULL, block_weights = NULL,
#' oblq_factors = NULL, init_Target = NULL, maxit = 20L, cutoff = 0,
#' normalization = "none", w = 1, random_starts = 1L, cores = 1L,
#' init = NULL, efa_control = NULL, rot_control = NULL, first_efa = NULL,
#' second_efa = NULL, verbose = TRUE)
#'
#' @description
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @param X Raw data matrix or correlation matrix.
#' @param n_generals Number of general factors to extract.
#' @param n_groups Number of group factors to extract.
#' @param method "GSLiD", "SL" (Schmid-Leiman), and "botmin" (bifactor-oblimin-target minimization). Defaults to "GSLiD".
#' @param cor Correlation method. Available correlations: c("pearson", "poly"). Defaults to "pearson".
#' @param estimator EFA fitting estimator: "ml" (maximum likelihood for multivariate normal variable), "uls" (minimum residuals), "pa" (principal axis) or "minrank" (minimum rank). Defaults to "uls".
#' @param projection Projection method. Available projections: "orth" (orthogonal), "oblq" (oblique) and "poblq" (partially oblique). Defaults to "oblq".
#' @param missing The way to handle missing data. Options: c("impute.mean", "impute.median", "complete.cases", "pairwise.complete.cases"). Defaults to "pairwise.complete.cases".
#' @param nobs Sample size. Defaults to NULL.
#' @param PhiTarget Target matrix for the factor correlations. Defaults to NULL.
#' @param PhiWeight Weight matrix for the factor correlations. Defaults to NULL.
#' @param init_Target Initial target matrix for the loadings. Defaults to NULL.
#' @param maxit Maximum number of iterations for the GSLiD algorithm. Defaults to 20L.
#' @param cutoff Cut-off used to update the target matrix upon each iteration. Defaults to 0.
#' @param normalization Available normalizations: "kaiser". Defaults to "none".
#' @param w \eqn{w} parameter for the extended target criterion ("xtarget"). Defaults to 1L.
#' @param random_starts Number of rotations with different random starting values. The rotation with the smallest cost function value is returned. Defaults to 1L.
#' @param blocks Vector with the number of factors for which separately applying the rotation criterion. Defaults to NULL.
#' @param block_weights Vector of weights for each block of factors.
#' @param oblq_factors Vector with the number of factors for each oblique block. E.g.: c(2, 4) means that there are two blocks of oblique factors: one block with 2 factors and another block with 4 factors. Everything else is orthogonal. Defaults to NULL.
#' @param cores Number of cores for parallel execution of multiple rotations. Defaults to 1L.
#' @param init Initial uniquenesses values for exploratory factor analysis estimation. Defaults to NULL.
#' @param efa_control List of control parameters for efa fitting. Defaults to NULL.
#' @param rot_control List of control parameters for the rotation algorithm. Defaults to NULL.
#' @param first_efa List of arguments to pass to \code{efast} to perform the first-order solution for the Schmid-Leiman method. Defaults to NULL.
#' @param second_efa List of arguments to pass to \code{efast} to perform the second-order solution for the Schmid-Leiman method. Defaults to NULL.
#' @param verbose Print the convergence progress information. Defaults to TRUE.
#'
#' @details
#'
#' If \code{efa.control = NULL}, then \code{list(maxit = 1e4)} is passed to \code{efa.control}. If \code{rot_control = NULL}, then \code{list(maxit = 1000, eps = 1e-05)} is passed to \code{rot_control}, where \code{eps} is the absolute tolerance. When the objective function does not make a larger improvement than \code{eps}, the algorithm is assumed to converge.
#'
#' If \code{Target} is provided but not \code{Weight}, then \code{Weight = 1 - Target} by default, which means a partially specified target rotation is performed. The same applies for \code{PhiTarget} and \code{PhiWeight}.
#'
#' If \code{init = NULL}, then the squared multiple correlations of each item with the remaining ones are used as initial values (These are known to be upper bounds).
#'
#' If \code{init_Target} is provided, then an initial target by means of the Schmid-Leiman transformation is not necessary.
#'
#' If \code{cutoff} is not 0, loadings smaller than such a cut-off are fixed to 0. When \code{cutoff} = 0, an empirical cut-off is used for each column of the loading matrix. They are the mean of the one-lagged differences of the sorted squared normalized loadings. Then, the target is determined by fixing to 0 the squared normalized loadings smaller than such cut-offs.
#'
#' @return List of class \code{bifactor}.
#' \item{efa}{List containing objects related to the exploratory factor analysis estimation. See \code{efast}.}
#' \item{bifactor}{List with the following components:}
#' \itemize{
#' \item loadings - Rotated loading matrix.
#' \item Phi - Factor correlation matrix.
#' \item T - Transformation matrix.
#' \item f - Objective value at the minimum.
#' \item iterations - Number of iterations performed by the rotation algorithm.
#' \item convergence - Convergence of the rotation algorithm.
#' \item uniquenesses - Vector of uniquenesses.
#' \item Rhat - Correlation matrix predicted by the model.
#' \item Target - Updated target matrix.
#' \item Weight - Weight matrix. It is the complement of the updated target.
#' \item GSLiD_iterations - Number of iterations performed by the GSLiD algorithm.
#' \item GSLiD_convergence - Convergence of the GSLiD algorithm.
#' \item min_congruences - Vector containing, for each iteration, the minimum Tucker's congruence between
#'  the current loading matrix and the previous loading matrix.
#' \item max_abs_diffs - Vector containing, for each iteration, the maximum absolute difference between the
#' current loading matrix and the previous loading matrix.
#' }
#'
#' \item{elapsed}{Total amount of time spent for execution (in nanoseconds).}
#'
#' @references
#'
#' Jiménez, M., Abad, F. J., Garcia-Garzon, E., & Garrido, L. E. (2023). Exploratory Bi-factor Analysis with Multiple General Factors. Multivariate behavioral research, 1–18. Advance online publication. https://doi.org/10.1080/00273171.2023.2189571
#'
#' @examples
#'
#' \dontrun{# Simulate data:
#' sim <- sim_factor(n_generals = 3, groups_per_general = 5, items_per_group = 6,
#' generals_rho = 0.3)
#' scores <- MASS::mvrnorm(1e4, rep(0, nrow(sim$R)), Sigma = sim$R)
#' s <- cor(scores)
#'
#' # Fit an exploratory bi-factor model with GSLiD:
#' fit <- bifactor(s, n_generals = 3, n_groups = 15, method = "GSLiD",
#' estimator = "uls", projection = "poblq", nobs = NULL, oblq_factors = 3,
#' random_starts = 10, cores = 8, w = 1, maxit = 20, verbose = TRUE, normalization = "none")
#'}
#'
#' @export
bifactor <- function(X, n_generals, n_groups, method = "GSLiD", cor = "pearson", estimator = "uls", projection = "oblq", missing = "pairwise.complete.cases", nobs = NULL, PhiTarget = NULL, PhiWeight = NULL, blocks = NULL, block_weights = NULL, oblq_factors = NULL, init_Target = NULL, maxit = 20L, cutoff = 0, normalization = "none", w = 1, random_starts = 1L, cores = 1L, init = NULL, efa_control = NULL, rot_control = NULL, first_efa = NULL, second_efa = NULL, verbose = TRUE) {
  .Call(`_latent_bifactor`, X, n_generals, n_groups, method, cor, estimator, projection, missing, nobs, PhiTarget, PhiWeight, blocks, block_weights, oblq_factors, init_Target, maxit, cutoff, normalization, w, random_starts, cores, init, efa_control, rot_control, first_efa, second_efa, verbose)
}

#' @title
#' Asymptotic standard errors for correlation matrices.
#' @usage
#'
#' asymp_cov(R, X = NULL, eta = 1, type = "normal")
#'
#' @description
#'
#' Get the asymptotic standard errors of correlation matrices of normal or arbitrary random deviates.
#'
#' @param R Correlation matrix.
#' @param X Optional raw data matrix.
#' @param eta Skewness parameter for elliptical data distributions.
#' @param type Type of random deviates: "normal", "elliptical" or "general".
#'
#' @details
#'
#' If \code{type = "normal"}, the calculation assumes that the raw data follows a multivariate normal distribution.
#' If \code{type = "elliptical"}, the calculation assumes that the raw data follows an elliptical distribution with skewness parameter \code{eta}.
#' If \code{type = "general"}, no assumption is made but need to provide the raw data via the \code{X} argument.
#'
#' @return The asymptotic covariance matrix of \code{R}.
#'
#' @references
#'
#' M.W. Browne and A. Shapiro (1986). The asymptotic covariance matrix of sample correlation coefficients under general conditions. Linear Algebra and its Applications, 82, 169-176. https://doi.org/10.1016/0024-3795(86)90150-3
#'
#' @export
asymp_cov <- function(R, X = NULL, eta = 1, type = "normal") {
  .Call(`_latent_asymp_cov`, R, X, eta, type)
}

#' @title
#' Standard errors for rotated factor loadings, factor correlations and uniquenesses.
#' @usage
#'
#' se(fit = NULL, n = NULL, X = NULL, type = "normal", eta = 1)
#'
#' @description
#'
#' Compute the sandwich standard errors of factor loadings, factor correlations and uniquenesses.
#'
#' @param fit Optional \code{efast} model.
#' @param n Sample size.
#' @param X Raw data matrix.
#' @param type Type of random deviates: "normal", "elliptical" or "general".
#' @param eta Skewness parameter for elliptical data distributions.
#'
#' @details
#'
#' Currently, only available for \code{estimator = uls}.
#'
#' @return A list with the standard errors of the rotated factor loadings, factor correlations and uniquenesses.
#'
#' @references
#'
#' Zhang G, Preacher KJ, Hattori M, Jiang G, Trichtinger LA (2019). A sandwich standard error estimator for exploratory factor analysis with nonnormal data and imperfect models. Applied Psychological Measurement, 43, 360–373. https://doi.org/10.1177/0146621618798669
#'
#' @export
se <- function(fit = NULL, nobs = NULL, X = NULL, type = "normal", eta = 1) {
  .Call(`_latent_se`, fit, nobs, X, type, eta)
}

#' @title
#' Hierarchical parallel analysis using either principal components (PCA) or principal axis factoring (PAF).
#' @usage
#'
#' parallel(X, nboot = 100L, cor = "pearson", quant = NULL, mean = TRUE, replace = FALSE,
#' PA = NULL, hierarchical = FALSE, efa = NULL, cores = 1L)
#'
#' @description
#'
#' Perform hierarchical parallel analysis to detect dimensionality using either principal components or principal axis factoring.
#'
#' @param X Raw data matrix.
#' @param nboot Number of bootstrap samples.
#' @param type Type of correlations: "pearson" or "poly".
#' @param missing The way to handle missing data. Options: c("impute.mean", "impute.median", "complete.cases", "pairwise.complete.cases"). Defaults to "pairwise.complete.cases".
#' @param quant Vector of quantiles of the distribution of bootstrap eigenvalues to which the compare the sample eigenvalues.
#' @param mean Logical. Compare the sample eigenvalues to the mean of the bootstrap eigenvalues. Defaults to TRUE.
#' @param replace Logical indicating whether the columns of \code{X} should be permuted with replacement.
#' @param PA Parallel analysis method. It can be either principal components ("PCA"), principal axis ("PAF") or both ("PCA" and "PAF"). Defaults to NULL, which sets c("PCA", "PAF").
#' @param hierarchical Logical indicating whether a second parallel analysis should be performed from the factor scores obtained after a first factor analysis analysis.
#' @param efa A list of arguments to pass to \code{efast} when \code{hierarchical = TRUE}.
#' @param cores Number of cores to perform the parallel bootstrapping.
#'
#' @details
#'
#' Not yet.
#'
#' @return A list with the bootstrapped eigenvalues and the estimated dimensionality.
#'
#' @references
#'
#' Horn, J. L. (1965). A Rationale and Test For the Number of Factors in Factor Analysis, Psychometrika, 30, 179-85. https://doi.org/10.1007/BF02289447
#'
#' @export
parallel <- function(X, nboot = 100L, cor = "pearson", missing = "pairwise.complete.cases", quant = NULL, mean = TRUE, replace = FALSE, PA = NULL, hierarchical = FALSE, efa = NULL, cores = 1L) {
  .Call(`_latent_parallel`, X, nboot, cor, missing, quant, mean, replace, PA, hierarchical, efa, cores)
}

#' @title
#' Fast polychoric correlations.
#' @usage
#'
#' polyfast(data, acov = "none", smooth = "pd", min_eigval = 0.001, nboot = 1000L, fit = FALSE, cores = 1L)
#'
#' @description
#'
#' Compute huge polychoric correlation matrices very fast.
#'
#' @param X Matrix of categorical scores. The lowest score must start at 0.
#' @param missing The way to handle missing data. Options: c("impute.mean", "impute.median", "complete.cases", "pairwise.complete.cases"). Defaults to "pairwise.complete.cases".
#' @param acov Use acov = 'cov' to obtain the asymptotic covariance matrix and acov = 'var' to simply obtain the asymptotic variances. Use "bootstrap" for estimating the asymptotic covariance matrix by resampling. Defaults to "none".
#' @param smooth Smooth the matrix to be positive definite ("pd"), positive semi-definite ("psd"), or estimate the maximum likely polychoric correlation matrix under the positive semi-definite constraint ("analytical"). Defaults to "none".
#' @param min_eigval Minimum eigenvalue when smooth = "pd". Defaults to 0.001.
#' @param nboot Number of bootstrap samples to compute the standard errors. It only works if acov = "bootstrap". Defaults to 1000L.
#' @param fit Should the fit value be calculated? Defaults to FALSE.
#' @param cores Number of parallel cores to compute the polychoric correlations.
#'
#' @details
#'
#' None yet.
#'
#' @return A list with the polychoric correlations, the thresholds, and the elapsed time in nanoseconds.
#' @export
polyfast <- function(data, missing = "pairwise.complete.cases", acov = "none", smooth = "none", min_eigval = 0.001, nboot = 1000L, fit = FALSE, cores = 1L) {
  .Call(`_latent_polyfast`, data, missing, acov, smooth, min_eigval, nboot, fit, cores)
}
