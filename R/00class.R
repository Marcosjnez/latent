# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/08/2026

# class definitions
#

setClass("latent",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           extra            = "list"
         )
)

setClass("llca",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           extra            = "list"
         ),
         contains = "latent"
)

setClass("lcfa",
         slots = c(
           version          = "character",
           call             = "call",
           timing           = "numeric",
           dataList         = "list",
           modelInfo        = "list",
           Optim            = "list",
           parameters       = "list",
           transformed_pars = "list",
           loglik           = "numeric",
           penalized_loglik = "numeric",
           loss             = "numeric",
           penalized_loss   = "numeric",
           extra            = "list"
         )
)
