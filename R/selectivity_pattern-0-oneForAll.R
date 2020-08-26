
fit_selectivity_0 = function(object, pattern, ...) {
  # create a list of model parameters

  out = list()
  fit = list()

  for(patt in pattern) {
    nm = as.character(patt)
    message("--- Selectivity pattern ", nm, "\n")
    FUN = match.fun(sprintf("fit_selectivity_%s", nm))
    outx  = FUN(object, ...)
    out[[nm]] = outx
    fit[[nm]] = .compare_fit(outx)
  }

  fit = as.data.frame(fit)

  best_model = apply(fit, 1, which.min)
  names(best_model) = rownames(object)

  xo = out[[1]]$selectivity
  output = out[[1]]$models
  best_pattern = numeric(length(best_model))

  message("Best models:\n")

  for(i in seq_along(best_model)) {
    this = out[[best_model[i]]]
    xo[i, ] = as.numeric(this$selectivity[i, ])
    output[[i]] = this$models[[i]]
    best_pattern[i] = this$pattern
    message(names(best_model)[i], ": ", this$pattern)
  }

  names(best_pattern) = names(best_model)

  output = list(selectivity=xo, models=output, y=object, x=x, pattern=best_pattern)

  return(output)

}


# Internal functions ------------------------------------------------------


.compare_fit = function(object) {
  fit = object$selectivity
  y = object$y
  out = rowSums((log(fit) - log(y))^2)
  return(out)
}

