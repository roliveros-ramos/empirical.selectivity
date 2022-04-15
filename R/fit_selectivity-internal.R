

.getBlockBreaks = function(object, n, method="equal", control=list()) {

  if(is.null(attr(object, "weights")))
    stop("No time dimension available for computing blocks.")

  yr = attr(object, "weights")$year

  if(method=="equal") {
    breaks = ceiling(quantile(yr, prob=seq(0, 1, len=n+1)))
    return(breaks)
  }

  min_block_size = if(!is.null(control$min_block_size)) control$min_block_size else 10
  breaks = .clusterBlocks(object, n, min_block_size)

  return(breaks)

  # method using clustering

}

.optimBlocks = function(object, breaks, w) {


  return(breaks)

}


