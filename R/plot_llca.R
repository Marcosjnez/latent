#' Plot method for latent class analysis (llca) objects
#'
#' @param x An object of class \code{"llca"}, typically from \code{latent}.
#' @param ... Additional arguments passed to \code{plot_lca_mixed}.
#' @return A list of ggplot objects (as from \code{plot_lca_mixed}).
#' @method plot llca
#' @export
plot.llca <- function(x, ...) {

  item_output <- latInspect(x, what = "item")
  item_types  <- x@dataList$item

  item_resp   <- reshape_lca_mixed(item_output, item_types)

  plot_lca_mixed(item_resp, ...)

}
