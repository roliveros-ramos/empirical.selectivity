

.getBlockBreaks = function(object, n, method="equal") {

  if(is.null(attr(object, "weights")))
    stop("No time dimension available for computing blocks.")

  yr = attr(object, "weights")$year

  if(method=="equal") {
    breaks = ceiling(quantile(yr, prob=seq(0, 1, len=n+1)))
    return(breaks)
  }

  # method using clustering

}

.optimBlocks = function(object, breaks, w) {


  return(breaks)

}

