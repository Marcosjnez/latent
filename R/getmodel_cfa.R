# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 22/05/2025

getmodel_cfa <- function(obj) {

  ng <- obj@Data@ngroups
  tech <- lavaan::lavTech(obj)
  par_mat <- lavaan::lavMatrixRepresentation(lavaan::partable(obj),
                                             representation = "LISREL")
  # mat_names <- subset(unique(par_mat$mat), unique(par_mat$mat) != "")
  mat_names <- c("lambda", "psi", "theta")

  matrices <- vector("list", length = ng)

  for(j in 1:ng) { ## per group

    matrices[[j]] <- tech[mat_names]
    par_group <- subset(par_mat, group == j)

    for(k in 1:length(mat_names)){ ## per matrix
      # matrices[[j]][[ mat_names[k] ]]

      par_mat2 <- subset(par_group, mat == mat_names[k])

      for(e in 1:nrow(par_mat2)){## per element of the matrix
        if(par_mat2[e,"label"] == ""){
          matrices[[j]][[ mat_names[k] ]][par_mat2[e,"row"], par_mat2[e,"col"]] <- par_mat2[e,"plabel"]
        }else{
          matrices[[j]][[ mat_names[k] ]][par_mat2[e,"row"], par_mat2[e,"col"]] <- par_mat2[e,"label"]
        }

        if(par_mat2[e,"free"] == 0 & !is.na(par_mat2[e,"ustart"])){
          matrices[[j]][[ mat_names[k] ]][par_mat2[e,"row"], par_mat2[e,"col"]] <- par_mat2[e,"ustart"]
        }


      }
    }

    ## fill in the lower diagonal of symetric matrices
    matrices[[j]]$psi[lower.tri(matrices[[j]]$psi)] <- t(matrices[[j]]$psi)[lower.tri(matrices[[j]]$psi)]
    matrices[[j]]$theta[lower.tri(matrices[[j]]$theta)] <- t(matrices[[j]]$theta)[lower.tri(matrices[[j]]$theta)]

  }

  # V <- list(group1 = matrices[[1]], group2 = matrices[[2]])

  if(ng > 1) {
    return(matrices)
  } else {
    return(matrices[[1]])
  }

}
