# Author: Mauricio Garnier-Villarreal
# Modified by: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/04/2026

# class definitions
#

setClass("llca",
         slots = c(
           version            = "character",
           call               = "call",
           timing             = "numeric",
           data_list          = "list",
           modelInfo          = "list",
           Optim              = "list",
           user_model         = "list",
           parameters         = "list",
           transformed_pars   = "list",
           posterior          = "matrix",
           state              = "vector",
           loglik             = "numeric",
           penalized_loglik   = "numeric",
           loglik_case        = "numeric",
           summary_table      = "data.frame",
           ClassConditional   = "list",
           RespConditional    = "list",
           probCat            = "list"
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
           penalized_loss     = "numeric"
         )
)

#setClass("lavaan",
#         contains = "lavaan"
#)

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

