

# class definitions
#

setClass("llca",
         slots = c(
           version            = "character", # version
           call               = "call", # matched call
           timing             = "list", # timing information
           modelInfo          = "list", # modelInfo
           Optim              = "list", # opt
           parameters         = "list",
           transformed_pars   = "list",
           posterior          = "list",
           state              = "list",
           loglik             = "list", # loglik values and info
           loglik_case        = "list",
           summary_table      = "list",
           ClassConditional   = "list",
           RespConditional    = "list",
           probCat            = "list"
         )
)



