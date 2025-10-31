# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 06/10/2025

# class definitions
#

setClass("llca",
         slots = c(
           version            = "character", # version
           call               = "call", # matched call
           timing             = "numeric", # timing information
           modelInfo          = "list", # modelInfo
           Optim              = "list", # opt
           user_model         = "list",
           parameters         = "list",
           transformed_pars   = "list",
           posterior          = "matrix",
           state              = "vector",
           loglik             = "numeric", # loglik values
           penalized_loglik   = "numeric", # penalized loglik
           loglik_case        = "numeric",
           summary_table      = "data.frame",
           ClassConditional   = "list",
           RespConditional    = "list",
           probCat            = "list"
         )
)

setClass("lcfa",
         slots = c(
           version            = "character", # version
           call               = "call", # matched call
           timing             = "numeric", # timing information
           modelInfo          = "list", # modelInfo
           Optim              = "list", # opt
           parameters         = "list",
           transformed_pars   = "list",
           loglik             = "numeric", # loglik values
           penalized_loglik   = "numeric", # penalized loglik
           loss               = "numeric",
           penalized_loss     = "numeric"
         )
)

#setClass("lcfa",
#         contains = "lavaan"
#)


