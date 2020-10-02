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

  if(any(is.na(w)))
    stop("NA in weights are not allowed.")

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
  xmod = attr(x, "model")
  if(!is.null(xmod)) {
    xmod$y = matrix(colSums(w*xmod$y), nrow=1)
    attr(output, "model") = xmod
  }

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

  attList$model$y = attList$model$y[i, ]
  attList$model$years = attList$model$years[i]

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
  storage.mode(f) = "integer"
  ind = split.default(seq_len(nrow(x)), f)
  out = lapply(ind, function(i) x[i])
  return(out)

}


