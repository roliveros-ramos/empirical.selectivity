
#' @export
fit_selectivity_27 = function(object, k=7, thr=1e-3, span=3, ...) {

  if(any(k<3)) {
    k = pmax(k, 3)
    warning("k must be greater or equal to 3, setting to 3.")
  }

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_27 = function(x, y, k, thr, span) {
    # knots, values, derivatives at extremes
    # create a list of model parameters
    x = as.numeric(x)
    y = as.numeric(y)
    if(all(y==0)) {
      output = list(fitted=y, model=NULL)
      return(output)
    }
    ind = .nonNullPoints(y, thr=thr, span=span)
    x0 = x[ind]
    y0 = y[ind]
    mod = fks(x0, log(y0), k = k-2, degree=3, prec=0)
    pred = predict(mod, newdata = data.frame(x=x), type="response")
    pred = exp(pred - max(pred, na.rm=TRUE))
    output = list(fitted=pred, x=x, y=y, model=mod)
    return(output)
  }

  x = as.numeric(colnames(object))

  xo = object*0
  out = vector("list", nrow(object))

  for(i in seq_len(nrow(object))) {
    cat("year =",i, "\n")
    tmp = .fit_selectivity_27(x=x, y=object[i, ], k=k, thr=thr, span=span)
    out[[i]] = tmp$model
    xo[i, ]  = tmp$fitted
  }

  output = list(selectivity=xo, models=out, y=object, x=x)

  return(output)
}
