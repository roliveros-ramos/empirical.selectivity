#' Fitting models to empirical selectivity data
#'
#' @param object An empirical selectivity object, extracted with
#' the \code{empirical_selectivity} function.
#' @param pattern Integer, the code for the selectivity model, as in SS.
#' @param blocks Integer, with the maximum number of blocks to be fitted,
#' or a vector with the year breaks for each block.
#' @param w A vector of weights (one per year) or a character
#' ('yield', 'Neff', 'Nsamp' or 'equal').
#' @param method Method for creating the blocks, current available methods are
#' 'equal' (equidistant placement) and 'cluster'.
#' @param control A list with additional options for every method. See details.
#' @param FUN Error function used for data fitting.
#' @param ... Additional arguments, not used currently.
#'
#' @export
fit_selectivity = function(object, pattern=0, blocks=NULL,
                           w=NULL, method=NULL, control=list(),
                           FUN=NULL, ...) {

  if(!inherits(object, "empirical_selectivity"))
    stop("Object to fit must be of class 'empirical_selectivity'.")

  if(is.null(method)) method = "equal"

  if(is.null(blocks)) {
    breaks = NULL
  } else {
    if(length(blocks)>1) {
      breaks = blocks
      blocks = length(breaks) - 1
    } else breaks = NULL
    if(blocks==0) blocks = NULL
  }

  if(is.null(breaks) & is.null(blocks) & !is.null(w))
    object = weighted.mean(object, w=w)

  if(!is.null(blocks)) {
    if(is.null(breaks)) {

      breaks = .getBlockBreaks(object, n=blocks, method=method,
                               control=control)
    }
    if(method=="optim") breaks = .optimBlocks(object, breaks, w)

    if(!is.null(attr(breaks, "n", exact=TRUE)))
      blocks = attr(breaks, "n", exact=TRUE) # temporal solution

    if(length(breaks)!=(blocks+1))
      stop(sprintf("Incorrect number of breaks, %d were expected.", blocks+1))

    # create new object for fitting using blocks
    object = split(object, breaks=breaks)
    if(is.null(w)) stop("Weighting 'w' must be specified")
    object = weighted.mean(object, w=w)
    nm = names(object)
    xmods = lapply(object, attr, which="model")
    ys = do.call(rbind, lapply(xmods, FUN="[[", i="y"))
    att = attributes(object[[1]])
    att$model$y = ys
    att$model$years = nm
    object = do.call(rbind, object)
    att$dim = dim(object)
    att$dimnames[[1]] = nm
    attributes(object) = att

  }

  if(identical(as.integer(pattern), 0L)) pattern = c(1, 24, 27)

  # error function
  if(is.null(FUN)) FUN = function(x, y) (log(x) - log(y))^2

  FUN = match.fun(FUN)
  if(!is.function(FUN)) stop("FUN must be a function.")

  output = fit_selectivity_0(object, pattern=pattern,
                             control=control, FUN=FUN, ...)

  attr(output, "fleet") = attr(object, "fleet")
  attr(output, "blocks") = ifelse(is.null(blocks), 0, blocks)
  attr(output, "breaks") = breaks
  # output includes a function to predict.
  class(output) = "selectivity_model"
  return(output)
}


#' Predict the selectivity for a new length or age vector.
#' @param object A 'selectivity_model' or 'empirical_selectivity_model' object.
#' @param x The length or age vector.
#' @param ... Additional arguments, currently not used.
#'
#' @export
predict.selectivity_model = function(object, x, ...) {
  out = c(x=list(x), lapply(object$models, FUN=predict, x=x))
  out = as.data.frame(out, optional = TRUE)
  class(out) = c("predict_empirical_selectivity", class(out))
  return(out)
}

#' @rdname predict.selectivity_model
#' @export
predict.empirical_selectivity_model = function(object, x, ...) {
  out = object$predict(x)
  names(out) = x
  return(out)
}



#' Title
#'
#' @param object
#' @param file
#' @param t
#' @param phase
#' @param ...
#'
#' @export
SS_writeselec = function(object, file=NULL, t=1, phase=2, ...) {

  FUN = match.fun(sprintf(".SS_writeselec_%d", object$pattern[t]))
  out = FUN(object, file=file, phase=phase, t=t, ...)
  return(invisible(out))

}


# Methods -----------------------------------------------------------------



#' @export
plot.selectivity_model = function(object, col="blue", cex=1, ylim=c(0,1.2), ...) {
  # standard plot function common to all methods
  if(attr(object, "blocks")>0) {
    nb = attr(object, "blocks")
    m1 = floor(sqrt(nb))
    m2 = ceiling(nb/m1)
    mfrow = sort(c(m1, m2), decreasing = TRUE)
    opar = par(no.readonly = TRUE)
    par(mfrow=mfrow, mar=c(3,3,2,1), oma=c(1,1,1,1))
    on.exit(par(opar))
    for(i in seq_len(nb)) {
      main = sprintf("yr = %s", rownames(object$selectivity)[i])
      plot(object$selectivity[i, ], main=main, col=col,
           ylim=ylim, cex=0.5*cex, ...)
      mod = attr(object$selectivity[i, ], "model")
      lines(mod, lwd=2, col="black", lty=1)
      points(object$x, object$y[i, ]/object$scale[i], pch=19, cex=0.75, col=col)
      if(object$pattern[i]==27)
        abline(v=object$models[[i]]$knots$knots, lty=3, col="red")
      legend("topleft",
             c("model", sprintf("fitted (pattern=%s)", object$pattern[i])),
             col=c("black", col), lty=1, lwd=c(3, 1), bty="n")
    }

    return(invisible())

  }

  plot(object$selectivity, col=col, ylim=ylim, ...)

  if(nrow(object$selectivity)==1) {
    points(object$x,object$y/object$scale, pch=19, cex=0.5*cex)
    if(object$pattern[1]==27)
      abline(v=object$models[[1]]$knots$knots, lty=3, col="red")
  }

  return(invisible())
}

#' @export
lines.selectivity_model = function(object, ...) {
  # standard plot function common to all methods
  if(nrow(object$selectivity)==1)
    lines(object$x, as.numeric(object$selectivity), ...)
  return(invisible())
}

#' @export
summary.selectivity_model = function(object, ...) {

  nblock  = nrow(object$fit)
  nmodels = ncol(object$fit)

  .tab = function(i, object) rbind(fit=object$fit[i,  ,drop=FALSE],
                                   npar=object$npar[i, , drop=FALSE])

  out = lapply(seq_len(nblock), FUN = .tab, object=object)

  names(out) = rownames(object$fit)

  class(out) = "summary.empirical_selectivity"
  return(out)

}

#' @export
coef.selectivity_model = function(object, ...) {

  nblock  = nrow(object$fit)
  nmodels = ncol(object$fit)

  .tab = function(i, object) object$models[[i]]$par

  out = lapply(seq_len(nblock), FUN = .tab, object=object)

  names(out) = rownames(object$fit)

  return(out)



}


#' @export
print.summary.empirical_selectivity = function(x, ...) {

  for(i in seq_along(x)) {
    cat("-- Block", names(x)[i], "\n")
    print(x[[i]])
  }

  return(invisible())

}

#' @export
plot.predict_empirical_selectivity = function(object, xlim=NULL, ylim=NULL, lwd=3, col=NULL, ...) {

  n = ncol(object)-1
  if(is.null(xlim)) xlim = range(object$x)
  if(is.null(ylim)) ylim = c(0,1.12)
  if(is.null(col)) col = seq_len(n)
  col = rep_len(col, n)
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  for(i in seq_len(n)) {
    mtext(colnames(object)[i+1], 3, line=-i-1, adj=0.05, col=col[i])
    lines(x=object$x, y=object[,i+1], col=col[i], lwd=lwd)
  }
  axis(1)
  axis(2, las=2)
  box()
  return(invisible(NULL))

}
