
# Class empirical selectivity ---------------------------------------------

#' @export
empirical_selectivity = function(object, ...) {
  UseMethod("empirical_selectivity")
}

#' @export
empirical_selectivity.default = function(object, thr=1e-5, ...) {
  # this function probably is gonna be used within SS_output
  # get the empirical selectivity for all fleets and all times,

  # if already exists, use it
  if(!is.null(object$empirical_selectivity) &
     inherits(object$empirical_selectivity, "SS_empirical_selectivity"))
    return(empirical_selectivity(object$empirical_selectivity, ...))

  # if does not exist, create it (for older versions)
  output = list(length = .getES(object, by="length", use="B", thr=thr),
                age    = .getES(object, by="age", use="B", thr=thr))

  output = unlist(output, recursive = FALSE)

  class(output) = "SS_empirical_selectivity"

  return(empirical_selectivity(output, ...))

}

#' @export
empirical_selectivity.SS_empirical_selectivity = function(object, fleet=NULL, sex=1, by="length") {

  by = match.arg(by, c("length", "age"))

  if(is.null(fleet)) return(object)

  xcode = paste(by, fleet, sex, sep=".")
  check = xcode %in% names(object)
  msg = sprintf("Empirical selectivity by %s for fleet %s and sex %d is not available",
                by, fleet, sex)
  if(!check) stop(msg)

  return(object[[xcode]])

}


#' @export
empirical_selectivity.SS_output = function(object, fleet=NULL, sex=1, by="length") {

  return(empirical_selectivity(object$empirical_selectivity, fleet=fleet, sex=sex, by=by))

}


# Methods for empirical_selectivity ---------------------------------------

#' @export
weighted.mean.empirical_selectivity = function(x, w, ...) {

  if(is.character(w)) {
    w = match.arg(w, c("catch", "Neff", "Nsamp", "equal"))
    lab = sprintf("%s (%s weighted)", attr(x, "fleet"), w)
    w = if(w=="equal") 1 else attr(x, "weights")[, w]
  } else lab = attr(x, "fleet")

  if(any(is.na(w))) stop("NA in weights are not allowed.")

  N = nrow(x)
  if(length(w)==1) w = rep(w, N)
  if(length(w)!=N) stop("One weight for each x row must be indicated.")

  rs = rowSums(x)
  w[rs==0] = 0 # rows with zero values get 0 weight
  w = w/sum(w) # weigths normalized to sum 1.

  output = colSums(x*w)
  output = output/max(output, na.rm=TRUE) # max to 1.
  output = matrix(output, nrow=1)
  colnames(output) = colnames(x)
  rownames(output) = lab

  attr(output, "fleet") = attr(x, "fleet")
  attr(output, "by") = attr(x, "by")
  attr(output, "model") = attr(x, "model")

  class(output) = c("empirical_selectivity", "SS_timexsize", "matrix")
  return(output)

}

#' @export
weighted.mean.SS_empirical_selectivity = function(x, w, ...) {

  output = lapply(x, FUN=weighted.mean, w=w)
  class(output) = "SS_empirical_selectivity"
  return(output)

}


# Methods

#' @export
'[.empirical_selectivity' = function(x, i, j, yr, age, length, drop = FALSE) {

  attList = attributes(x)
  class(x) =  "matrix"

  if(!missing(yr)) {
    years = suppressWarnings(as.numeric(rownames(x)))
    if(all(is.na(years))) stop("No year data is available for indexing.")
    i = na.omit(match(yr, years))
  }

  if(!missing(age)) {
    if(attr(x, "by")!="age") stop("No age data is available for indexing.")
    ages = suppressWarnings(as.numeric(colnames(x)))
    # extract ages withing the range
    j = which(ages >= min(age, na.rm=TRUE) & ages <= max(age, na.rm=TRUE))
  }

  if(!missing(length)) {
    if(attr(x, "by")!="length") stop("No length data is available for indexing.")
    lengths = suppressWarnings(as.numeric(colnames(x)))
    # extract lengths withing the range
    j = which(lengths >= min(length, na.rm=TRUE) & lengths <= max(length, na.rm=TRUE))
  }

  out = '['(x, i, j, drop=drop)
  class(out) = c("empirical_selectivity", "SS_timexsize", "matrix")

  attr(out, "by") = attList$by
  attr(out, "fleet") = attList$fleet
  attr(out, "model") = attList$model
  if(!is.null(attList$weights)) attr(out, "weights") = attList$weights[i,]

  return(out)
}

#' @export
plot.empirical_selectivity = function(x, type=1, col="blue", p=1.5,
                                      xlim=NULL, ylim=NULL, ...) {

  if(nrow(x)==1) type = 2

  units = if(attr(x, "by")=="length") "cm" else "years"
  lab = sprintf(if(attr(x, "by")=="length") "Size (%s)" else "Age (%s)", units)

  switch (type,
    "1" = plot_TimexSize_type1(x=x, col=col, lab=lab, xlim=xlim, ylim=ylim, ...),
    "2" = plot_TimexSize_type2(x=x, col=col, p=p, lab=lab, xlim=xlim, ylim=ylim, ...),
    "3" = plot_TimexSize_type3(x=x, col=col, p=p, lab=lab, xlim=xlim, ylim=ylim, ...),
    stop("Invalid plot type."))

  return(invisible())
}


#' @export
plot_TimexSize_type1 = function(x, col, alpha, lab, xlim, ylim, ...) {

  if(missing(alpha)) alpha = max(min(3/nrow(x), 0.9), 0.1)

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(is.null(xlim)) xlim = range(y, na.rm=TRUE)
  if(is.null(ylim)) ylim = c(0,1)

  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  title(xlab=lab, main=attr(z, "fleet"))
  for(i in seq_len(nrow(z))) {
    lines(y, as.numeric(z[i,]), col=.makeTransparent(alpha, col))
  }
  lines(y, weighted.mean(z, w="catch"), lwd=2, col="black", lty=3)
  lines(y, weighted.mean(z, w="Nsamp"), lwd=2, col="red", lty=1)
  lines(y, weighted.mean(z, w="Neff"), lwd=2, col="red", lty=3)
  lines(y, weighted.mean(z, w="equal"), lwd=2, col=col, lty=1)
  lines(attr(z, "model"), lwd=3, col="black", lty=1)
  axis(1)
  axis(2, las=1)
  box()

  legend("topleft", c("catch", "N (sample)", "N (effective)", "equal", "model"),
         col=c(1,2,2,col,1), lty=c(3,1,3,1,1), lwd=c(2,2,2,2,3), bty="n")

  return(invisible())

}


#' @export
plot_TimexSize_type2 = function(x, col, p, lab, xlim, ylim, ...) {

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(nrow(z)==1) {
    msg = if(is.na(x)) rownames(z) else sprintf("year = %d", x)
    plot(y, z, type="l", xlab=lab, ylab="", main=msg,
         las=1, col=col, xlim=xlim, ylim=ylim, ...)
    return(invisible())
  }

  col = directionalPalette(col=col, p=p)
  if(is.null(xlim)) xlim = range(x, na.rm=TRUE)
  if(is.null(ylim)) ylim = range(y, na.rm=TRUE)

  image.plot(x, y, as.matrix(z), xlab="Time", ylab=lab, las=1, col=col,
             xlim=xlim, ylim=ylim, ...)
  return(invisible())

}

#' @export
plot_TimexSize_type3 = function(x, col, p, lab, xlim, ylim, ...) {

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(nrow(z)==1) {
    msg = if(is.na(x)) rownames(z) else sprintf("year = %d", x)
    plot(y, z, type="l", xlab=lab, ylab="", main=msg,
         las=1, col=col, xlim=xlim, ylim=ylim, ...)
    return(invisible())
  }

  col = directionalPalette(col=col, p=p)
  if(is.null(xlim)) xlim = range(x, na.rm=TRUE)
  if(is.null(ylim)) ylim = range(y, na.rm=TRUE)

  z = as.matrix(z)
  z[is.na(z)] = 0

  mountains(xvec=y, yvec=x, zmat=z, xlab=lab, ylab="Time", las=1, col=col,
            ...)
  return(invisible())

}

# model uses the SS code
#' @export
fit_selectivity = function(object, pattern = 27, ...) {
  # this function creates all the parameters (e.g. knots, values)
  # use a switch for every pattern
  if(!inherits(object, "empirical_selectivity"))
    stop("Object to fit must be of class 'empirical_selectivity'.")

  output = switch(as.character(pattern),
     '1' = fit_selectivity_1(object, ...),
    '24' = fit_selectivity_24(object, ...),
    '27' = fit_selectivity_27(object, ...),
    stop(sprintf("Selectivity pattern %s not implemented.", pattern)))

  attr(output, "fleet") = attr(object, "fleet")
  # output includes a function to predict.
  class(output) = "selectivity_model"
  return(output)
}

#' @export
predict.selectivity_model = function(object, x, ...) {
  # use the internal function(size/age)
}

#' @export
plot.selectivity_model = function(object, ...) {
  # standard plot function common to all methods
  plot(object$selectivity, ...)
  if(nrow(object$selectivity)==1) {
    points(object$x,object$y, pch=19, cex=0.5)
    abline(v=object$models[[1]]$knots$knots, lty=3, col="red")
  }
  return(invisible())
}

#' @export
lines.selectivity_model = function(object, ...) {
  # standard plot function common to all methods
  if(nrow(object$selectivity)==1)
    lines(object$x, as.numeric(object$selectivity),
          ...)
  return(invisible())
}

#' @export
SS_writeselec = function(object, file=NULL, phase=2, t=1, ...) {

  FUN = match.fun(sprintf(".SS_writeselec_%d", object$pattern))
  out = FUN(object, file=file, phase=phase, t=t, ...)
  return(invisible(out))

  }


SS_run = function(...) {

}

