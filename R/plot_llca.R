# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 02/02/2026
#'
#' @title
#' Plot regression coefficients
#' @description
#'
#' Plot the coefficients of latent models.
#'
#' @usage
#'
#' plot(fit)
#'
#' @param fit model fitted with lca.
#' @param confidence Coverage of the confidence interval.
#'
#' @details Obtain a plot with the coefficients of latent models.
#'
#' @return List with the following objects:
#' \item{vcov}{Variance-covariance matrix between the parameters.}
#' \item{se}{Standard errors.}
#' \item{SE}{Standard errors in the model list.}
#'
#' @references
#'
#' None yet.
#'
#' @method plot llca
#' @export
plot.llca <- function(fit,
                      type = "standard",
                      what = "OR",
                      effects = "coding",
                      confidence = 0.95,
                      show_est_ci = TRUE,
                      est_ci_header_cex = 0.5,
                      cex_y = 0.5,
                      mfrow = c(1, 1),
                      ...) {

  # Get standard errors:
  SE <- se(fit = fit, type = type, model = "model", digits = 9)
  select_betas <- match(fit@modelInfo$lca_trans$beta, colnames(SE$vcov))
  select_betas <- select_betas[!is.na(select_betas)]

  #### log(odds ratio) ####

  if(effects == "coding") {

    EF <- effects_coding(fit@transformed_pars$beta,
                         SE$vcov[select_betas, select_betas])
    betas <- EF$beta_new
    vcov <- EF$vcov_new
    colnames(betas) <- colnames(fit@transformed_pars$beta)

  } else if (effects == "dummy") {

    betas <- fit@transformed_pars$beta[, -1] # Remove the intercept
    vcov <- SE$vcov[select_betas, select_betas]

  } else {
    stop("Available effects: 'coding' and 'dummy'")
  }

  # Create new names for the parameters:
  allnames <- paste(rep(rownames(betas), times = ncol(betas)),
                    rep(colnames(betas), each = nrow(betas)),
                    sep = "|")
  parameters <- c(betas)
  names(parameters) <- allnames
  rownames(vcov) <- colnames(vcov) <- allnames

  group <- factor(rep(colnames(betas), each = nrow(betas)))
  coef_names <- rep(rownames(betas), times = ncol(betas))

  # Compute standard errors:
  se <- sqrt(diag(vcov))

  # Compute confidence intervals (asumming asymptotic normality):
  lower <- parameters - sqrt(qchisq(confidence, df = 1)) * se
  upper <- parameters + sqrt(qchisq(confidence, df = 1)) * se

  if(what == "log") {

    xlab <- "log(odds ratio) (95% CI)"
    est_ci_header <- "log(odds ratio) (95% CI)"
    refline <- 0

  } else if(what == "OR") {

    xlab <- "odds ratio (95% CI)"
    est_ci_header <- "odds ratio (95% CI)"
    refline <- 1

    parameters <- exp(parameters)

    # Compute standard errors:
    J <- diag(parameters)
    vcov <- J %*% vcov %*% t(J)
    se <- sqrt(diag(vcov))

    # Compute confidence intervals (asumming asymptotic normality):
    lower <- exp(lower)
    upper <- exp(upper)

  } else {
    stop("Available what: 'log' and 'OR'")
  }

  # Returned list
  result <- list(se = se,
                 vcov = vcov,
                 upper = upper,
                 lower = lower,
                 parameters = parameters)

  par(mfrow = mfrow)

  # Plot
  forestplot(parm = parameters,
             lower = lower, upper = upper,
             group = group,
             labels = coef_names,
             refline = refline,
             show_est_ci = show_est_ci,
             est_ci_header_cex = est_ci_header_cex,
             cex_y = cex_y,
             ...)

  invisible(result)

}

