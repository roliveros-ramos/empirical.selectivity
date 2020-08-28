#' Title
#'
#' @param object
#'
#' @param pattern
#' @param blocks
#' @param w
#' @param method
#' @param control
#' @param ...
#'
#' @export
fit_selectivity = function(object, pattern=0, blocks=NULL,
                           w=NULL, method=NULL, control=list(), ...) {
  # this function creates all the parameters (e.g. knots, values)
  # use a switch for every pattern
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
      # equal is hardcoded temporarily
      breaks = .getBlockBreaks(object, n=blocks, method="equal")
    }
    if(method=="optim") breaks = .optimBlocks(object, breaks, w)

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

  output = fit_selectivity_0(object, pattern=pattern, control=control, ...)

    # switch(as.character(pattern),
    #               '0' = fit_selectivity_0(object, pattern=patts, control=control, ...),
    #               '1' = fit_selectivity_1(object, control=control, ...),
    #               '24' = fit_selectivity_24(object, control=control, ...),
    #               '27' = fit_selectivity_27(object, control=control, ...),
    #               stop(sprintf("Selectivity pattern %s not implemented.", pattern)))

  attr(output, "fleet") = attr(object, "fleet")
  attr(output, "blocks") = ifelse(is.null(blocks), 0, blocks)
  # output includes a function to predict.
  class(output) = "selectivity_model"
  return(output)
}


#' Title
#' @param object
#'
#' @param x
#' @param ...
#'
#' @export
predict.selectivity_model = function(object, x, ...) {
  # use the internal function(size/age)
}

#' Title
#' @param object
#'
#' @param file
#' @param phase
#' @param t
#' @param ...
#'
#' @export
SS_writeselec = function(object, file=NULL, phase=2, t=1, ...) {

  FUN = match.fun(sprintf(".SS_writeselec_%d", object$pattern[t]))
  out = FUN(object, file=file, phase=phase, t=t, ...)
  return(invisible(out))

}


#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
SS_run = function(...) {

}



# Methods -----------------------------------------------------------------



#' @export
plot.selectivity_model = function(object, col="blue", ...) {
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
      plot(object$selectivity[i, ], main=main, col=col, ...)
      mod = attr(object$selectivity[i, ], "model")
      lines(mod, lwd=2, col="black", lty=1)
      points(object$x, object$y[i, ], pch=19, cex=0.75, col=col)
      if(object$pattern[i]==27)
        abline(v=object$models[[i]]$knots$knots, lty=3, col="red")
      legend("topleft",
             c("model", sprintf("fitted (pattern=%s)", object$pattern[i])),
             col=c("black", col), lty=1, lwd=c(3, 1), bty="n")
    }

    return(invisible())

  }

  plot(object$selectivity, ...)

  if(nrow(object$selectivity)==1) {
    points(object$x,object$y, pch=19, cex=0.5)
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

  .tab = function(i, object) rbind(fit=object$fit[i, ], npar=object$npar[i, ])

  out = lapply(seq_len(nblock), FUN = .tab, object=object)

  names(out) = rownames(object$fit)

  class(out) = "summary.empirical_selectivity"
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
