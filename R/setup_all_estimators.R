setup_all_estimators <- function(estimators, arguments) {

  result <- list()

  for(i in 1:length(estimators)) {

    estimator <- estimators[i]
    args <- arguments[[i]]
    result[[i]] <- setup_estimator(estimator, args)

  }

  return(result)

}
