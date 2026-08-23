# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 23/08/2026

#### Retrieve the unrotated CFA model ####

lefa_source_lcfa <- function(fit) {

  if(!inherits(fit, "lefa")) {
    stop("fit must inherit from class 'lefa'.")
  }

  source <- fit@extra$efa

  if(is.null(source)) {

    source_models <- fit@extra[
      vapply(fit@extra,
             FUN = \(x) inherits(x, "lcfa"),
             FUN.VALUE = logical(1L))
    ]

    if(length(source_models) != 1L) {
      stop("The lefa object does not contain one unrotated lcfa model.")
    }

    source <- source_models[[1L]]

  }

  if(!inherits(source, "lcfa")) {
    stop("The unrotated model stored in the lefa object is not an lcfa object.")
  }

  #### Result ####

  return(source)

}

#### Rotated parameter blocks ####

lefa_rotated_block_names <- function(fit) {

  data_param <- fit@modelInfo$data_param

  result <- unique(c(data_param$lambda_group,
                     data_param$psi_group,
                     data_param$alpha_group))
  result <- intersect(result,
                      names(fit@modelInfo$trans))

  #### Result ####

  return(result)

}

lefa_rotation_block_names <- function(fit) {

  data_param <- fit@modelInfo$data_param

  result <- intersect(data_param$X_group,
                      names(fit@modelInfo$trans))

  #### Result ####

  return(result)

}

lefa_rotated_parameters <- function(fit) {

  blocks <- lefa_rotated_block_names(fit)
  result <- fit@modelInfo$trans[blocks]

  #### Result ####

  return(result)

}

#### Estimator label ####

lefa_estimator_label <- function(fit) {

  source <- lefa_source_lcfa(fit)

  result <- switch(
    source@dataList$estimator,
    ml = "Maximum likelihood",
    fml = "Maximum likelihood",
    uls = "Unweighted least squares",
    dwls = "Diagonally weighted least squares",
    source@dataList$estimator
  )

  #### Result ####

  return(result)

}

#### Rotated parameter table ####

lefa_parameter_table <- function(fit) {

  parameters <- lefa_rotated_parameters(fit)
  labels <- unique(unname(c(unlist(parameters))))

  if(length(labels) == 0L) {

    result <- data.frame(
      label = character(0L),
      estimate = numeric(0L),
      se = numeric(0L),
      z = numeric(0L),
      pvalue = numeric(0L),
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    #### Result ####

    return(result)

  }

  trans_labels <- fit@modelInfo$transparameters_labels
  transparameters <- fit@Optim$transparameters

  if(is.null(names(transparameters))) {
    names(transparameters) <- trans_labels
  }

  estimates <- transparameters[labels]

  standard_errors <- tryCatch(
    se(fit,
       parameters = parameters,
       digits = NULL)$se,
    error = function(e) {
      warning("Standard errors for the rotated parameters could not be computed: ",
              conditionMessage(e))
      setNames(rep(NA_real_, length(labels)), labels)
    }
  )

  if(is.null(names(standard_errors))) {
    names(standard_errors) <- labels
  }

  standard_errors <- standard_errors[labels]
  z <- estimates/standard_errors
  z[!is.finite(z) | standard_errors <= 0] <- NA_real_
  pvalue <- 2*stats::pnorm(-abs(z))

  result <- data.frame(
    label = labels,
    estimate = as.numeric(estimates),
    se = as.numeric(standard_errors),
    z = as.numeric(z),
    pvalue = as.numeric(pvalue),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  #### Result ####

  return(result)

}
