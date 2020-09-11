
.clusterBlocks = function(object, n, thr) {

  clus = kmeans(object, centers=n, nstart = 100)

  kc = as.numeric(factor(clus$cluster,
                         levels=unique(clus$cluster),
                         labels = 1:n))

  kc2 = .findBlocks(kc, n=thr)

  yr = attr(object, "weights")$year

  breaks = c(sort(tapply(yr, INDEX = kc2, FUN=min)), max(yr))

  attr(breaks, "cluster") = data.frame(year=yr, cluster=kc, block=kc2)
  attr(breaks, "n") = length(table(kc2))

  return(breaks)

}



# Internal functions ------------------------------------------------------

.findBlocks = function(x, n) {

  S = max(x)
  spp = seq_len(S)

  thr = n/length(x)

  out = matrix(0, ncol=S, nrow=length(x))

  for(iSp in spp) {
    out[, iSp] = 0 + (x == iSp)
  }
  colnames(out) = seq_len(S)

  out0 = apply(out, 2, .getRates, n=n)
  out1 = apply(out, 2, .getBlocks, n=n)
  out2 = rle(rowSums(out1))
  ones = out2$values == 1
  out2$values[ones] =  seq_len(sum(ones))
  out2$values[!ones] = NA
  out2 = inverse.rle(out2)
  if(sum(!is.na(out2))<n) {
    # warning("No blocks identified, returning NA")
    return(out2)
  }
  out3 = round(approx(x=out2, xout = seq_along(out2), rule=2)$y, 0)

  spx = table(out3, x)
  sp0 = setNames(apply(spx, 1, which.max), nm=NULL)
  gg0 = setNames(seq_along(unique(sp0)), nm = unique(sp0))
  out4 = setNames(gg0[as.character(sp0[out3])], nm=NULL)
  out4 = .removeBumps(out4, thr=thr)

  mainSp = apply(table(out4, factor(x, levels=spp)), 1, which.max)

  out = setNames(mainSp[as.character(out4)], nm=NULL)

  return(out)

}


.getBlocks = function(x, n) {
  x1 = 0 + (.getRates(x=x, n=n) > 0.5)
  x1 = .completeTails(x1)
  return(x1)
}

.getRates = function(x, n) {
  x0 = as.numeric(filter(x, filter=rep(1,n)/n))
  return(x0)
}

.removeBumps = function(x, thr) {

  x = rle(x)

  px = x$lengths/sum(x$lengths)
  x$values[px<thr] = NA
  validGroups = unique(na.omit(x$values))
  validGroups = setNames(seq_along(validGroups), nm=validGroups)
  validGroups = validGroups[as.character(x$values)]
  names(validGroups) = NULL
  x$values = validGroups
  x = inverse.rle(x)
  x = round(approx(x=x, xout = seq_along(x), rule=2)$y, 0)

  x = rle(x)
  if(length(x$values)<3) {
    x$values = seq_along(x$values)
    x = inverse.rle(x)
    return(x)
  }

  while(any(table(x$values)>1)) {
    bumps = numeric(length(x$values))
    for(i in seq_len(length(x$values)-2)) {
      # if(x$values[i]==x$values[i+2]) x$values[i+1] = x$values[i]
      if(x$values[i]==x$values[i+2]) bumps[i+1] = 1
    }
    thisFirst = which.min(x$lengths/bumps)
    x$values[thisFirst] = x$values[thisFirst+1]
    x = inverse.rle(x)
    x = rle(x)
  }
  # x = inverse.rle(x)
  # x = rle(x)
  x$values = seq_along(x$values)
  x = inverse.rle(x)

  return(x)
}


.completeTails = function(x) {
  nn = length(x)
  nona = which(!is.na(x))
  h = head(nona, 1)
  t = tail(nona, 1)
  x[seq_len(h-1)] = x[h]
  x[t+seq_len(nn-t)] = x[t]
  return(x)
}


