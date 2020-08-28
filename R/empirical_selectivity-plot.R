
#' @export
plot.empirical_selectivity = function(x, type=1, col="blue", p=1.5,
                                      xlim=NULL, ylim=NULL, main=NULL, ...) {

  if(nrow(x)==1) type = 2

  units = if(attr(x, "by")=="length") "cm" else "years"
  lab = sprintf(if(attr(x, "by")=="length") "Size (%s)" else "Age (%s)", units)

  switch (type,
          "1" = plot_TimexSize_type1(x=x, col=col, lab=lab, xlim=xlim, ylim=ylim, main=main, ...),
          "2" = plot_TimexSize_type2(x=x, col=col, p=p, lab=lab, xlim=xlim, ylim=ylim, main=main, ...),
          "3" = plot_TimexSize_type3(x=x, col=col, p=p, lab=lab, xlim=xlim, ylim=ylim, main=main, ...),
          stop("Invalid plot type."))

  return(invisible())
}

# Internal ----------------------------------------------------------------

plot_TimexSize_type1 = function(x, col, alpha, lab, xlim, ylim, main, ...) {

  if(missing(alpha)) alpha = max(min(3/nrow(x), 0.9), 0.1)

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(is.null(xlim)) xlim = range(y, na.rm=TRUE)
  if(is.null(ylim)) ylim = c(0,1)

  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  title(xlab=lab, main=ifelse(is.null(main), attr(z, "fleet"), main))

  if(!is.null(attr(z, "weights"))) {
    cols = rep(.makeTransparent(alpha, col), nrow(z))
  } else {
    cols = seq_len(nrow(z)) + 1
  }

  for(i in seq_len(nrow(z))) {
    lines(y, as.numeric(z[i,]), col=cols[i])
  }
  if(!is.null(attr(z, "weights"))) {
    lines(y, weighted.mean(z, w="catch"), lwd=2, col="black", lty=3)
    lines(y, weighted.mean(z, w="Nsamp"), lwd=2, col="red", lty=1)
    lines(y, weighted.mean(z, w="Neff"), lwd=2, col="red", lty=3)
    lines(y, weighted.mean(z, w="equal"), lwd=2, col=col, lty=1)

    legend("topleft", c("catch", "N (sample)", "N (effective)", "equal", "model"),
           col=c(1,2,2,col,1), lty=c(3,1,3,1,1), lwd=c(2,2,2,2,3), bty="n")

  } else {

    legend("topleft", c("model",rownames(z)),
           col=c(1, cols), lty=1, lwd=c(3, rep(1, nrow(z))), bty="n")

  }
  mod = attr(z, "model")
  mod$y = colMeans(mod$y)
  lines(mod, lwd=3, col="black", lty=1)
  axis(1)
  axis(2, las=1)
  box()

  return(invisible())

}


plot_TimexSize_type2 = function(x, col, p, lab, xlim, ylim, main, ...) {

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(nrow(z)==1) {
    if(is.null(main))
      main = if(is.na(x)) rownames(z) else sprintf("year = %d", x)
    plot(y, z, type="l", xlab=lab, ylab="", main=main,
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

plot_TimexSize_type3 = function(x, col, p, lab, xlim, ylim, main, ...) {

  z = x
  x = suppressWarnings(as.numeric(rownames(z)))
  y = suppressWarnings(as.numeric(colnames(z)))

  if(nrow(z)==1) {
    if(is.null(main))
      main = if(is.na(x)) rownames(z) else sprintf("year = %d", x)
    plot(y, z, type="l", xlab=lab, ylab="", main=main,
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
