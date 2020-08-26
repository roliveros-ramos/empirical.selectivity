
# Class empirical selectivity ---------------------------------------------

#' Title
#' @param object
#' @param fleet
#' @param sex
#' @param by
#'
#' @param ...
#'
#' @export
empirical_selectivity = function(object, ...) {
  UseMethod("empirical_selectivity")
}

#' @describeIn empirical_selectivity
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

#' @describeIn empirical_selectivity
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

#' @describeIn empirical_selectivity
#' @export
empirical_selectivity.SS_output = function(object, fleet=NULL, sex=1, by="length") {

  return(empirical_selectivity(object$empirical_selectivity, fleet=fleet, sex=sex, by=by))

}


# Methods for empirical_selectivity ---------------------------------------

#' Title
#' @param x
#'
#' @param w
#' @param ...
#'
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

#' @export
weighted.mean.list = function(x, w, ...) {

  output = lapply(x, FUN=weighted.mean, w=w)
  return(output)

}


# Methods -----------------------------------------------------------------

#' Title
#' @param x
#'
#' @param i
#' @param j
#' @param yr
#' @param age
#' @param length
#' @param drop
#'
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


#' Title
#' @param x
#'
#' @param f
#' @param breaks
#' @param drop
#' @param sep
#' @param lex.order
#' @param ...
#'
#' @export
split.empirical_selectivity = function (x, f, breaks, drop = FALSE, sep = ".",
                                        lex.order = FALSE, ...) {

  if(!missing(breaks)) {
    if(is.null(attr(x, "weights")))
      stop("No time dimension available for computing blocks.")
    yr = attr(x, "weights")$year
    f = cut(yr, breaks = breaks, right=FALSE, include.lowest = TRUE)
  }

  if (is.list(f))
    f = interaction(f, drop = drop, sep = sep, lex.order = lex.order)
  else if(!is.factor(f))
    f = as.factor(f)
  else if(drop)
    f = factor(f)
  storage.mode(f) <- "integer"
  ind = .Internal(split(seq_len(nrow(x)), f))
  out = lapply(ind, function(i) x[i])
  return(out)

}


#' Title
#' @param object
#'
#' @param pattern
#' @param blocks
#' @param breaks
#' @param w
#' @param ...
#'
#' @export
fit_selectivity = function(object, pattern = 27, blocks=NULL, breaks=NULL,
                           w=NULL, method=NULL, ...) {
  # this function creates all the parameters (e.g. knots, values)
  # use a switch for every pattern
  if(!inherits(object, "empirical_selectivity"))
    stop("Object to fit must be of class 'empirical_selectivity'.")

  if(is.null(method)) method = "equal"

  if(is.null(breaks) & is.null(blocks) & !is.null(w))
    object = weighted.mean(object, w=w)

  if(!is.null(breaks) & is.null(blocks)) blocks = length(breaks) - 1

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
    att = attributes(object[[1]])
    nm = names(object)
    object = do.call(rbind, object)
    att$dim = dim(object)
    att$dimnames[[1]] = nm
    attributes(object) = att

  }

  if(identical(as.integer(pattern), 0L)) patts = c(1, 24, 27)

  if(length(pattern)>1) {
    patts = pattern
    pattern = 0
  }

  output = switch(as.character(pattern),
     '0' = fit_selectivity_0(object, pattern=patts, ...),
     '1' = fit_selectivity_1(object, ...),
    '24' = fit_selectivity_24(object, ...),
    '27' = fit_selectivity_27(object, ...),
    stop(sprintf("Selectivity pattern %s not implemented.", pattern)))

  attr(output, "fleet") = attr(object, "fleet")
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

