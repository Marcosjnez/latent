# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/08/2026

#### Display labels for substantive groups ####

lcfa_group_display_labels <- function(fit) {

  labels <- as.character(fit@dataList$group_label)

  if(length(labels) == 1L && (is.na(labels) || labels == "")) {
    labels <- "group1"
  }

  empty <- is.na(labels) | labels == ""
  labels[empty] <- paste0("group", which(empty))

  #### Result ####

  return(labels)

}

#### Extract estimator-level fit components ####

lcfa_fit_components <- function(fit) {

  controls <- fit@modelInfo$control_estimator
  outputs <- tryCatch(
    fit@Optim$outputs$estimators$doubles,
    error = function(e) NULL
  )

  if(is.null(outputs) || length(outputs) != length(controls)) {
    stop("Estimator outputs are unavailable or do not match the fitted model.")
  }

  result <- vector("list", length(controls))

  infer_group <- function(control) {

    group_index <- control$group_index

    if(is.numeric(group_index) &&
       length(group_index) == 1L &&
       is.finite(group_index) &&
       group_index >= 1L &&
       group_index <= fit@dataList$ngroups) {
      return(as.integer(group_index))
    }

    if(fit@dataList$ngroups == 1L) {
      return(1L)
    }

    covariance_name <- control$covariance_name

    if(!is.null(covariance_name)) {

      suffixes <- paste0(".", fit@dataList$group_label)
      matched <- which(vapply(
        suffixes,
        FUN = \(suffix) grepl(suffix, covariance_name, fixed = TRUE),
        FUN.VALUE = logical(1L)
      ))

      if(length(matched) == 1L) {
        return(matched)
      }

    }

    return(NA_integer_)

  }

  for(i in seq_along(controls)) {

    control <- controls[[i]]
    values <- as.numeric(outputs[[i]])
    role <- control$role

    if(is.null(role)) {
      role <- if(startsWith(control$estimator, "cfa")) {
        "cfa"
      } else {
        "penalty"
      }
    }

    group_index <- infer_group(control)

    if(is.na(group_index)) {
      stop("The substantive group of estimator ", i,
           " could not be identified. Add group_index metadata when the ",
           "estimator is constructed.")
    }

    component <- list(
      estimator = control$estimator,
      role = role,
      group_index = group_index,
      group_label = fit@dataList$group_label[group_index],
      pattern_index = control$pattern_index,
      covariance_name = control$covariance_name,
      means_name = control$means_name,
      values = values
    )

    if(identical(role, "cfa")) {

      if(length(values) < 7L) {
        stop("CFA estimator ", i,
             " did not return the expected seven fit components.")
      }

      component$fit <- setNames(
        values[seq_len(7L)],
        c("loss", "loss_base", "loss_sat",
          "loglik", "loglik_base", "loglik_sat",
          "penalty")
      )

    } else {

      component$fit <- c(
        loss = 0,
        loss_base = 0,
        loss_sat = 0,
        loglik = 0,
        loglik_base = 0,
        loglik_sat = 0,
        penalty = if(length(values) > 0L) values[1L] else 0
      )

    }

    result[[i]] <- component

  }

  #### Result ####

  return(result)

}

#### Saturated incomplete-data model for direct FIML fit statistics ####

lcfa_direct_h1 <- function(fit) {

  if(!isTRUE(fit@dataList$direct_fiml)) {

    #### Result ####

    return(NULL)

  }

  control <- list(
    rstarts = 1L,
    cores = 1L,
    print = FALSE
  )

  group <- if(fit@dataList$ngroups > 1L) {
    fit@dataList$group
  } else {
    NULL
  }

  result <- lmvnorm(
    data = fit@dataList$data,
    group = group,
    variables = fit@dataList$item_label,
    se = FALSE,
    do.fit = TRUE,
    message = FALSE,
    control = control
  )

  if(!isTRUE(result@Optim$convergence) ||
     !is.finite(result@Optim$f)) {
    stop("The saturated incomplete-data H1 model did not converge.")
  }

  #### Result ####

  return(result)

}

lmvnorm_loglik_by_group <- function(fit) {

  outputs <- fit@Optim$outputs$estimators$doubles
  pattern_counts <- lengths(fit@dataList$patterns)
  result <- numeric(length(pattern_counts))

  offset <- 0L

  for(i in seq_along(pattern_counts)) {

    count <- pattern_counts[i]
    index <- offset+seq_len(count)

    result[i] <- sum(vapply(
      outputs[index],
      FUN = \(x) as.numeric(x[[2L]]),
      FUN.VALUE = numeric(1L)
    ))

    offset <- offset+count

  }

  names(result) <- if(length(result) > 1L) {
    fit@dataList$group_label
  } else {
    "group1"
  }

  #### Result ####

  return(result)

}

#### Independence-model log-likelihood ####

lcfa_independence_loglik <- function(fit) {

  ngroups <- fit@dataList$ngroups
  result <- numeric(ngroups)

  for(i in seq_len(ngroups)) {

    X <- fit@dataList$data_per_group[[i]]

    if(!is.null(X) && nrow(X) > 0L) {

      loglik <- 0

      for(j in seq_len(ncol(X))) {

        x <- X[[j]]
        x <- x[is.finite(x)]

        if(length(x) < 2L) {
          loglik <- NA_real_
          break
        }

        mu <- mean(x)
        variance <- mean((x-mu)^2)

        if(!is.finite(variance) || variance <= 0) {
          loglik <- NA_real_
          break
        }

        loglik <- loglik-0.5*sum(
          log(2*pi)+log(variance)+(x-mu)^2/variance
        )

      }

      result[i] <- loglik

    } else {

      S <- fit@dataList$sample.cov[[i]]
      n <- fit@dataList$nobs[[i]]
      variances <- diag(S)

      if(any(!is.finite(variances)) ||
         any(variances <= 0)) {
        result[i] <- NA_real_
      } else {
        result[i] <- -0.5*n*sum(
          log(2*pi)+log(variances)+1
        )
      }

    }

  }

  names(result) <- lcfa_group_display_labels(fit)

  #### Result ####

  return(result)

}

#### Aggregate fit components by substantive group ####

lcfa_fit_matrix <- function(fit, compute_h1 = TRUE) {

  components <- lcfa_fit_components(fit)
  group_labels <- lcfa_group_display_labels(fit)
  rows <- c("loss", "loss_base", "loss_sat",
            "loglik", "loglik_base", "loglik_sat",
            "penalty")

  result <- matrix(
    0,
    nrow = length(rows),
    ncol = fit@dataList$ngroups,
    dimnames = list(rows, group_labels)
  )

  missing_value <- matrix(
    FALSE,
    nrow = length(rows),
    ncol = fit@dataList$ngroups,
    dimnames = dimnames(result)
  )

  for(component in components) {

    group_index <- component$group_index

    for(row in rows) {

      value <- component$fit[[row]]

      if(!is.finite(value)) {
        missing_value[row, group_index] <- TRUE
      } else {
        result[row, group_index] <-
          result[row, group_index]+value
      }

    }

  }

  result[missing_value] <- NA_real_

  likelihood_model <- fit@dataList$estimator %in%
    c("ml", "fml")

  normal_likelihood <- is.null(fit@dataList$likelihood) ||
    identical(tolower(fit@dataList$likelihood), "normal")

  if(likelihood_model && normal_likelihood) {

    result["loglik_base", ] <-
      lcfa_independence_loglik(fit)

    if(isTRUE(fit@dataList$direct_fiml) &&
       compute_h1) {

      h1 <- tryCatch(
        lcfa_direct_h1(fit),
        error = function(e) {
          warning("The saturated incomplete-data model could not be fitted: ",
                  conditionMessage(e))
          NULL
        }
      )

      if(!is.null(h1)) {
        result["loglik_sat", ] <-
          lmvnorm_loglik_by_group(h1)
      }

    }

  } else if(likelihood_model && !normal_likelihood) {

    warning("Likelihood-based lcfa fit indices are currently available only ",
            "for likelihood = 'normal'.")
    result[c("loglik", "loglik_base", "loglik_sat"), ] <- NA_real_

  }

  result <- rbind(
    result,
    penalized_loss =
      result["loss", ]+result["penalty", ],
    penalized_loglik =
      result["loglik", ]-result["penalty", ]
  )

  overall <- vapply(
    seq_len(nrow(result)),
    FUN = \(i) {
      x <- result[i, ]
      if(any(!is.finite(x))) NA_real_ else sum(x)
    },
    FUN.VALUE = numeric(1L)
  )

  result <- cbind(result, overall = overall)

  #### Result ####

  return(result)

}

#### Model residuals ####

lcfa_residuals <- function(fit) {

  data_param <- fit@dataList$data_param
  group_labels <- lcfa_group_display_labels(fit)
  result <- vector("list", fit@dataList$ngroups)
  names(result) <- group_labels

  for(i in seq_len(fit@dataList$ngroups)) {

    model_name <- data_param$model_group[i]
    meanshat_name <- data_param$meanshat_group[i]

    model_cov <- fit@transformed_pars[[model_name]]
    sample_cov <- fit@dataList$sample.cov[[i]]

    covariance <- sample_cov-model_cov
    dimnames(covariance) <- dimnames(model_cov)

    means <- NULL
    thresholds <- NULL

    if(fit@dataList$cor == "pearson" &&
       isTRUE(fit@dataList$meanstructure)) {

      sample_mean <- fit@dataList$sample.mean[[i]]
      model_mean <- c(fit@transformed_pars[[meanshat_name]])
      names(model_mean) <- rownames(fit@transformed_pars[[meanshat_name]])
      means <- sample_mean[names(model_mean)]-model_mean

    }

    if(fit@dataList$cor == "poly") {

      tauhat_name <- data_param$tauhat_group[i]
      model_thresholds <- c(fit@transformed_pars[[tauhat_name]])
      names(model_thresholds) <- rownames(
        fit@transformed_pars[[tauhat_name]]
      )

      sample_thresholds <- lcfa_order_thresholds(
        threshold_blocks = fit@transformed_pars[data_param$taus_group[[i]]],
        model_tau = fit@dataList$model[[i]]$tau
      )
      names(sample_thresholds) <- names(model_thresholds)
      thresholds <- sample_thresholds-model_thresholds

    }

    result[[i]] <- list(
      covariance = covariance,
      means = means,
      thresholds = thresholds
    )

  }

  #### Result ####

  return(result)

}

#### Parameter table ####

lcfa_parameter_table <- function(fit) {

  labels <- fit@modelInfo$parameters_labels
  estimates <- fit@Optim$parameters
  names(estimates) <- labels

  standard_errors <- tryCatch(
    fit@Optim$SE$se,
    error = function(e) NULL
  )

  if(is.null(standard_errors) ||
     length(standard_errors) != length(labels)) {

    standard_errors <- tryCatch(
      se(fit, parameters = labels, digits = NULL)$se,
      error = function(e) {
        warning("Standard errors could not be computed: ",
                conditionMessage(e))
        setNames(rep(NA_real_, length(labels)), labels)
      }
    )

  }

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
