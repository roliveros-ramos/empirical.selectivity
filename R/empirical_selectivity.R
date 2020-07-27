
# Class empirical selectivity ---------------------------------------------

empirical_selectivity = function(object, ...) {
  UseMethod("empirical_selectivity")
}

empirical_selectivity.default = function(object, ...) {
  # this function probably is gonna be used within SS_output
  # get the empirical selectivity for all fleets and all times,

  lbins = object$lbins
  years = seq(from=object$startyr, to=object$endyr, by=1)

  nlend = object$lendbase
  natlen = object$natlen
  natlen = natlen[natlen$'Beg/Mid'=="B" & natlen$Era=="TIME", ]
  natlen = natlen[order(natlen$Yr), ]

  nbase = expand.grid(Yr=years, Sex=c(1,2))
  nbase = merge(nbase, natlen, all=TRUE)
  nbase = nbase[order(nbase$Yr), ]

  yr = sort(unique(natlen$Yr))
  natlen_1 = natlen[natlen$Sex==1, -(1:12)]
  natlen_2 = natlen[natlen$Sex==2, -(1:12)]

  natlen = as.matrix(natlen_1) + as.matrix(natlen_2)
  natlen = natlen[, as.character(lbins)]

  xx = split(nlend, f = list(nlend$Fleet, nlend$Sex))
  output = lapply(xx, FUN = .mta, bins=lbins, years=years, scale=natlen)

  class(output) = "SS_empirical_selectivity"

  return(output)

}

empirical_selectivity.SS_empirical_selectivity = function(object, fleet=NULL, sex=1) {

  if(is.null(fleet)) return(object)

  xcode = paste(fleet, sex, sep=".")
  check = xcode %in% names(object)
  msg = sprintf("Empirical selectivity for fleet %s and sex %d is not available",
                fleet, sex)
  if(!check) stop(msg)

  return(object[[xcode]])

}

# aggregation is done outside.

weighted.mean.empirical_selectivity = function(x, w, ...) {

  if(any(is.na(w))) stop("NA in weights are not allowed.")

  N = nrow(x)
  if(length(w)==1) w = rep(w, N)
  if(length(w)!=N) stop("One weight for each x row must be indicated.")

  rs = rowSums(x)
  w[rs==0] = 0 # rows with zero values get 0 weight
  w = w/sum(w) # weigths normalized to sum 1.

  output = colSums(x*w)
  output = output/max(output) # max to 1.
  output = matrix(output, nrow=1)
  colnames(output) = colnames(x)
  rownames(output) = "weighted.mean"

  class(output) = c("empirical_selectivity", "SS_timexsize", "matrix")
  return(output)

}

weighted.mean.SS_empirical_selectivity = function(x, w, ...) {

  output = lapply(x, FUN=weighted.mean, w=w)
  class(output) = "SS_empirical_selectivity"
  return(output)

}


# Methods

'[.empirical_selectivity' = function(x, i, j, drop = FALSE) {

  attList = attributes(x)
  class(x) =  "matrix"
  out = '['(x, i, j, drop=drop)
  class(out) = c("empirical_selectivity", "SS_timexsize", "matrix")
  return(out)
}

plot.empirical_selectivity = function(x, type=1, col="blue", p=1.5, w=1, ...) {

  if(nrow(x)==1) type = 1

  switch (type,
    "1" = plot_TimexSize_type1(x=x, col=col, p=p, ...),
    "2" = plot_TimexSize_type2(x=x, col=col, w=w, ...))


  return(invisible())
}

plot_TimexSize_type1 = function(x, col, p, ...) {

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(nrow(z)==1) {
    msg = if(is.na(x)) rownames(z) else sprintf("year = %d", x)
    plot(y, z, type="l", xlab="Time", ylab="Size (cm)", main=msg,
         las=1, col=col,
         ...)
    return(invisible())
  }

  col = directionalPalette(col=col, p=p)
  image.plot(x, y, z, xlab="Time", ylab="Size (cm)", las=1, col=col,
             ...)
  return(invisible())

}

plot_TimexSize_type2 = function(x, col, w, alpha, ...) {

  if(missing(alpha)) alpha = min(3/nrow(x), 0.9)
  zz = weighted.mean(x, w=w)
  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  plot.new()
  plot.window(xlim=range(y), ylim=c(0,1), xlab="Time", ylab="Size (cm)")
  for(i in seq_len(nrow(z))) {
    lines(y, as.numeric(z[i,]), col=.makeTransparent(alpha, col))
  }
  lines(y, zz, lwd=2, col=col)
  axis(1)
  axis(2, las=1)
  box()

  return(invisible())

}

# model uses the SS code
fit_selectivity = function(object, pattern = 27, ...) {
  # this function creates all the parameters (e.g. knots, values)
  # use a switch for every pattern
  if(!inherits(object, "empirical_selectivity"))
    stop("Object to fit must be of class 'empirical_selectivity'.")

  output = switch(as.character(pattern),
    '24' = fit_selectivity_24(object, ...),
    '27' = fit_selectivity_27(object, ...),
    stop(sprintf("Selectivity pattern %s not implemented.", pattern)))

  # output includes a function to predict.
  class(output) = "selectivity_model"
  return(output)
}

predict.selectivity_model = function(object, x, ...) {
  # use the internal function(size/age)
}

plot.selectivity_model = function(object, x=NULL, ...) {
  # standard plot function common to all methods
  plot(object$selectivity, ...)
  if(nrow(object$selectivity)==1) points(object$x,object$y, pch=19, cex=0.5)
  return(invisible())
}


# For each pattern

fit_selectivity_27 = function(object, k=7, thr=0.05, span=3, ...) {

  if(any(k<3)) {
    k = pmax(k, 3)
    warning("k must be greater or equal to 3, setting to 3.")
  }

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_27 = function(x, y, k, thr, span, ...) {
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
    mod = fks(x0, log(y0), k = k-2, degree=3)
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

fit_selectivity_24 = function(object, ...) {
  # create a list of model parameters
  return(output)
}

SS_writeselec = function(object, ...) {
  # write a selectivity_model object into lines for a ctl file.
  # can be printed in the console or to a file.

  return(invisible(output))
}

SS_run = function(...) {

}

