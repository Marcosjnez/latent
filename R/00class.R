

# class definitions
#

setClass("llca",
         slots = c(
           version            = "character", # version
           call               = "call", # matched call
           timing             = "numeric", # timing information
           modelInfo          = "list", # modelInfo
           Optim              = "list", # opt
           parameters         = "list",
           transformed_pars   = "list",
           posterior          = "matrix",
           state              = "vector",
           loglik             = "numeric", # loglik values
           penalized_loglik   = "numeric", # penalized loglik
           loglik_case        = "numeric",
           summary_table      = "list",
           ClassConditional   = "list",
           RespConditional    = "list",
           probCat            = "list"
         )
)



