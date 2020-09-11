
directionalPalette = function(n=64, col="green3", p=0.8, x=NULL, ...) {

  alpha = seq(from=0, to=1, length=n)^p

  cols = sapply(alpha, FUN=.makeTransparent, col=col)

  if(!is.null(x)) {
    ind = cut(x, breaks = n, labels = FALSE)
    return(cols[ind])
  }

  return(unlist(cols))

}


# Internal ----------------------------------------------------------------

.makeTransparent = function(alpha, col) {
  col = col2rgb(col)
  alpha = floor(255*alpha)
  rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
}
