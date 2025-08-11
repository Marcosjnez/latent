## use latInspect everywhere we can:
## add more things as we go

latInspect <- function(object,
                       what = "classconditional",
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

    round(object@posterior, digits)

  } else if (what == "state"){

    object@state

  }else if (what == "loglik_case" ||
            what == "casell"){

    object@loglik_case

  } else if (what == "loglik" ||
             what == "ll" ||
             what == "LL"){

    object@loglik

  } else if (what == "probcat"){

    temp <- object@probCat

    for(j in 1:length(temp)){
      temp[[j]] <- round(temp[[j]], digits)
    }

  } else if (what == "timing"){

    object@timing

  }


}

