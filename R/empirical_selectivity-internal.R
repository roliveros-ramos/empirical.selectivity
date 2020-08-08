

.mta = function(x, bins, years, var="Obs", FUN=mean, scale=1, by="length") {

  FUN = match.fun(FUN)
  out = matrix(0, nrow=length(years), ncol=length(bins))
  rownames(out) = years
  colnames(out) = bins

  xout  = tapply(X = x[, var], INDEX=list(x$Yr, x$Bin), FUN=FUN)

  obins = as.numeric(colnames(xout))
  oyear = as.numeric(rownames(xout))

  rrows = match(oyear, years)
  rcols = match(obins, bins)

  out[rrows, rcols] = xout

  out = out/scale
  xmax = apply(out, 1, max, na.rm=TRUE)
  xmax[xmax==0] = 1
  out = out/xmax

  out[is.nan(out)] = 0

  names(dimnames(out)) = c("year", by)

  attr(out, "by") = by

  class(out) = c("empirical_selectivity", "SS_timexsize", "matrix")

  return(out)

}


.nonNullPoints0 = function(y) {
  # copy y
  yx = y
  yx[yx<thr] = thr # round to thr

  # remove isolated 0s
  y0 = rle(yx==thr)
  y0$values[y0$lengths==1 & y0$values] = FALSE
  yx = inverse.rle(y0)
  # remove isolated 1s
  y0 = rle(yx)
  y0$values[y0$lengths==1 & !y0$values] = TRUE
  yx = inverse.rle(y0)
  y0 = rle(yx)
  y1 = y0
  y1$lengths[y1$values] = -1
  ind = which.max(y1$lengths)
  ind = cumsum(y0$lengths[c(ind-1, ind)])
  ind = seq(from=ind[1]-1, to=ind[2]+1, by=1)
  ind = pmin(pmax(ind, 1), length(y))
  return(ind)
}

.nonNullPoints = function(y, thr, span) {
  # copy y

  if(length(span)==1) span = c(span, span)
  if(any(span<0)) stop("span must be positive.")

  yx = cumsum(y)/sum(y)

  ind0 = which.min(yx<thr/2)
  ind1 = which.max(yx>=(1-thr/2))

  ind = c(ind0, ind1)
  ind = seq(from=ind[1]-span[1], to=ind[2]+span[2], by=1)
  ind = pmin(pmax(ind, 1), length(y))
  ind = sort(unique(ind))
  return(ind)
}

.getES = function(object, by="length", use="B") {

  db  = sprintf("%sdbase", substr(by, 1, 3))
  nat = sprintf("nat%s", substr(by, 1, 3))
  years = seq(from=object$startyr, to=object$endyr, by=1)

  naged  = object[[db]]
  natage = object[[nat]]

  natage = natage[natage$'Beg/Mid'==use & natage$Era=="TIME", ]
  natage = natage[order(natage$Yr), ]

  nbase = expand.grid(Yr=years, Sex=c(1,2))
  nbase = merge(nbase, natage, all=TRUE)
  nbase = nbase[order(nbase$Yr), ]

  yr = sort(unique(natage$Yr))
  natage_1 = natage[natage$Sex==1, -(1:12)]
  natage_2 = natage[natage$Sex==2, -(1:12)]

  natage = as.matrix(natage_1) + as.matrix(natage_2)

  abins = sort(as.numeric(colnames(natage))) # use pop bins

  natage = natage[, as.character(abins)]

  output = lapply(split(naged, f = list(naged$Fleet, naged$Sex)),
                  FUN = .mta, bins=abins, years=years, scale=natage, by=by)

  for(i in seq_along(output))
    attr(output[[i]], "fleet") = object$FleetNames[i]

  return(output)

}


