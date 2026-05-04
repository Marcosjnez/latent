# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 04/05/2026

# class definitions
#

setClass("llca",
         slots = c(
           version            = "character",
           call               = "call",
           timing             = "numeric",
           dataList           = "list",
           modelInfo          = "list",
           Optim              = "list",
           parameters         = "list",
           transformed_pars   = "list",
           loglik             = "numeric",
           penalized_loglik   = "numeric",
           loss               = "numeric",
           penalized_loss     = "numeric",
           extra              = "list"
         )
)

setClass("lcfa",
         slots = c(
           version            = "character",
           call               = "call",
           timing             = "numeric",
           dataList           = "list",
           modelInfo          = "list",
           Optim              = "list",
           parameters         = "list",
           transformed_pars   = "list",
           loglik             = "numeric",
           penalized_loglik   = "numeric",
           loss               = "numeric",
           penalized_loss     = "numeric",
           extra              = "list"
         )
)

setClass("latent",
         slots = c(
           version            = "character",
           call               = "call",
           timing             = "numeric",
           dataList           = "list",
           modelInfo          = "list",
           Optim              = "list",
           parameters         = "list",
           transformed_pars   = "list",
           loglik             = "numeric",
           penalized_loglik   = "numeric",
           loss               = "numeric",
           penalized_loss     = "numeric",
           extra              = "list"
         )
)

