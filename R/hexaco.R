#' Example dataset: hexaco
#'
#' A dataset from https://osf.io/72zp3/.
#'
#' @format A data frame with 5 different samples (28143 rows) and 100 HEXACO items:
#' \describe{
#'   \item{x}{An integer from 1 to 5}
#'   \item{y}{A lowercase letter}
#' }
#' @source Generated example
#'
#' @usage
#'
#' "hexaco"
#'
#' @examples
#'
#' \dontrun{
#' # Check samples and sample sizes
#' samples <- unique(hexaco$sample) # industry mooc fire student dutch
#' Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
#' names(Ns) <- samples
#'
#' # Subset the items pertaining to the HEXACO-100
#' selection <- 5:104
#' full <- hexaco[, selection]
#'
#' mooc <- full[hexaco$sample == samples[2], ]
#' dim(mooc)
#'
#' model.EM <- "FEA ~= hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
#'              ANX ~= hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
#'              DEP ~= hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
#'              SEN ~= hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"
#'
#' fit <- lcfa(model = model.EM, data = mooc,
#'             ordered = TRUE, estimator = "ml", do.fit = TRUE,
#'             control = NULL)
#' fit@loglik # -0.283407
#' fit@penalized_loglik # -0.283407
#' fit@loss # 0.1574787
#' fit@Optim$opt$iterations
#' fit@Optim$opt$convergence
#' fit@timing
#' }
"hexaco"
