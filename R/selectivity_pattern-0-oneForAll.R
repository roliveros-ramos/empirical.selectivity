
fit_selectivity_0 = function(object, pattern, ...) {
  # create a list of model parameters

  out  = list()
  fit  = list()
  npar = list()

  for(patt in pattern) {
    nm = as.character(patt)
    message("--- Selectivity pattern ", nm, "\n")
    FUN = match.fun(sprintf("fit_selectivity_%s", nm))
    outx  = FUN(object, ...)
    out[[nm]] = outx
    fit[[nm]] = outx$fit
    npar[[nm]] = outx$npar
  }

  x = as.numeric(colnames(object))

  fit = as.data.frame(fit, check.names=FALSE)
  npar = as.data.frame(npar, check.names=FALSE)

  best_model = apply(fit, 1, which.min)
  names(best_model) = rownames(object)

  xo = out[[1]]$selectivity
  output = out[[1]]$models
  best_pattern = out[[1]]$pattern
  names(best_pattern) = rownames(object)

  if(length(pattern)>1) {
    message("Best models:\n")
    for(i in seq_along(best_model)) {
      this = out[[best_model[i]]]
      xo[i, ] = as.numeric(this$selectivity[i, ])
      output[[i]] = this$models[[i]]
      best_pattern[i] = this$pattern[i]
      message(names(best_model)[i], ": ", this$pattern[i])
    }
  }

  output = list(selectivity=xo, models=output, y=object, x=x, pattern=best_pattern,
                fit=fit, npar=npar)

  return(output)

}


# Internal functions ------------------------------------------------------


# .compare_fit = function(object) {
#   fit = object$selectivity
#   y = object$y
#   out = rowSums((log(fit) - log(y))^2)
#   return(out)
# }

.compare_fit = function(fit, y) {
  out = rowSums((log(fit) - log(y))^2)
  return(out)
}




