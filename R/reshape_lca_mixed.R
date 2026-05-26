
#' Split and reshape mixed (multinomial + Gaussian) LCA item output
#'
#' @param item_output Result of \code{latInspect(fit, "item")}.
#' @param item_types Character vector of indicator types (e.g., \code{fit@dataList$item}).
#' @return A list with components \code{multinomial} and \code{gaussian}, each a data frame.
#' @export
reshape_lca_mixed <- function(item_output, item_types) {
  stopifnot(length(item_output) == length(item_types))

  is_multinom <- item_types == "multinomial"
  is_gauss    <- item_types == "gaussian"

  mult_items <- item_output[is_multinom]
  gauss_items<- item_output[is_gauss]

  mult_df <- if (length(mult_items) > 0) reshape_lca_multinomial(mult_items) else data.frame()
  gauss_df<- if (length(gauss_items) > 0) reshape_lca_continuous(gauss_items)  else data.frame()

  list(multinomial = mult_df, gaussian = gauss_df)
}


reshape_lca_multinomial <- function(object) {
  do.call(rbind, lapply(names(object), function(varname) {
    mat <- object[[varname]]
    # Create all combinations of response options and latent classes
    df <- expand.grid(
      response = rownames(mat),
      class    = colnames(mat),
      stringsAsFactors = FALSE
    )
    df$probability <- as.vector(mat)
    df$variable <- varname
    # Return with desired column order
    df[, c("variable", "response", "class", "probability")]
  }))
}

reshape_lca_continuous <- function(object) {
  do.call(rbind, lapply(names(object), function(varname) {
    mat <- object[[varname]]
    data.frame(
      variable = varname,
      class    = colnames(mat),
      mean     = mat["mean", ],
      sd       = mat["stdv", ],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}
