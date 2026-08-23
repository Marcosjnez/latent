# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026
#'
#' Exploratory Factor Analysis Class
#'
#' The \code{"lefa"} class stores a rotated exploratory factor analysis and
#' inherits the deterministic multistep covariance machinery. The unrotated
#' \code{lcfa} model is retained in the \code{extra} slot.
#'
#' @name lefa-class
#' @rdname lefa-class
setClass("lefa",
         contains = "multistep")
