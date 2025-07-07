

setMethod("show", "llca", function(object) {

  # Print header with model name and version
  cat(sprintf("%s %s ended normally after %d iterations\n\n",
              "latent", as.character( packageVersion('latent') ),
              object@Optim$iterations))

  # Print Estimator, Optimization, and Parameters section
  cat(sprintf("  %-45s %s\n", "Estimator", "Multinomial"))
  cat(sprintf("  %-45s %s\n", "Optimization method", object@Optim$control$opt))
  cat(sprintf("  %-45s %d\n\n", "Number of model parameters", object@modelInfo$nparam))

  # Print Number of Observations
  cat(sprintf("  %-45s %d\n\n", "Number of observations", object@modelInfo$nobs))
  cat(sprintf("  %-45s %d\n\n", "Number of response patterns", object@modelInfo$npatterns))

  # Print Model Test Section
  chisq <- -2*object@loglik
  pval <- 1 - pchisq(chisq, df = object@modelInfo$df)
  cat("Model Test User Model:\n")
  cat("  ", paste(rep("-", 54), collapse = ""), "\n\n", sep = "")
  cat(sprintf("  %-45s %.3f\n", "Test statistic", chisq))
  cat(sprintf("  %-45s %d\n", "Degrees of freedom", object@modelInfo$df))
  cat(sprintf("  %-45s %.3f\n", "P-value (Chi-square)", pval))

  invisible(object)

})


inspect.llca <- function(object, ...) {
  lInspect.llca(object, ...)
}

setMethod("inspect", "llca", lInspect.llca)

## use lInspect everywhere we can:
## add more things as we go

lInspect.llca <- function(object,
                          what = "free",
                          digits = 3){


  # object must inherit from class llca
  stopifnot(inherits(object, "llca"))



  # be case insensitive
  what <- tolower(what)


  ### result objects
  if (what == "classconditional" ||
      what == "profile") {

    temp <- object@ClassConditional

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(temp)


  } else if (what == "classes"){

    temp <- round(object@transformed_pars$classes, digits)
    names(temp) <- paste0("Class", 1:length(temp))

    return(temp)

  } else if (what == "respconditional") {

    temp <- object@RespConditional

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

    return(temp)

  } else if (what == "convergence") {

    object@Optim$convergence

  } else if (what == "data"){

    object@Optim$data_list$dt

  } else if (what == "posterior"){

    object@posterior

  } else if (what == "state"){

    object@state

  }else if (what == "loglik_case" ||
             what == "casell"){

    object@loglik_case

  } else if (what == "loglik" ||
             what == "ll"){

    object@loglik

  } else if (what == "probcat"){

    object@probCat

  } else if (what == "timing"){

    object@timing

    }


}


