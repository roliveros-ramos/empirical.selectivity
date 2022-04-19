
# Class empirical selectivity ---------------------------------------------

#' @title Empirical selectivity diagnostic
#' @description Function to compute the empirical selectivity diagnostic for SS
#' outputs and methods to extract them for a particular fleet and time.
#' @param object Typically, the result of reading Stock Synthesis model outputs
#' with \code{SS_output}.
#' @param fleet Integer, fleet number to compute the empirical selectivity.
#' @param sex If the fleet targets a specific sex, defaults to 1.
#' @param by Either 'length' or 'age', defaults to 'length'.
#' @param ... Additional arguments for specific methods, see Details.
#' @details This function is called on run time by \code{SS_output}, so the
#' empirical selectivities are already computed (for length and age) in
#' your 'rep' object. An additional call to the function can be used to extract
#' the empirical selectivities for a particular fleet or sex.
#' @examples \dontrun{
#'
#' rep = BigEye_2020
#' es = empirical_selectivity(rep, fleet=2, by="length")
#' }
#' @export
empirical_selectivity = function(object, ...) {
  UseMethod("empirical_selectivity")
}

#' @describeIn empirical_selectivity General constructor
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
                age    = .getES(object, by="age", use="B", thr=thr),
                gsize  = .getES(object, by="gsize", use="B", thr=thr))

  output = unlist(output, recursive = FALSE)

  class(output) = "SS_empirical_selectivity"

  return(empirical_selectivity(output, ...))

}

#' @describeIn empirical_selectivity Extract from empirical_selectivity objects
#' @export
empirical_selectivity.SS_empirical_selectivity = function(object, fleet=NULL, sex=1, by="length", ...) {

  by = match.arg(by, c("length", "age", "gsize"))

  if(is.null(fleet)) return(object)

  xcode = paste(by, fleet, sex, sep=".")
  check = xcode %in% names(object)
  msg = sprintf("Empirical selectivity by %s for fleet %s and sex %d is not available",
                by, fleet, sex)
  if(!check) stop(msg)

  return(object[[xcode]])

}

#' @describeIn empirical_selectivity Extracts from SS_output object
#' @export
empirical_selectivity.SS_output = function(object, fleet=NULL, sex=1, by="length", ...) {

  return(empirical_selectivity(object$empirical_selectivity, fleet=fleet, sex=sex, by=by, ...))

}

#' @describeIn empirical_selectivity Transform a matrix in a empirical_selectivity object.
#' @export
empirical_selectivity.matrix = function(object, fleet=NULL, sex=1, by="length",
                                        years=NULL, bins=NULL, weights=NULL,
                                        model=NULL,  ...) {

  if(is.null(fleet)) stop("Fleet name ('character') must be provided")
  if(is.null(years)) stop("Years (rows) must be provided")
  if(is.null(bins)) stop("Bin marks (columns) must be provided")

  out = object
  colnames(out) = bins
  rownames(out) = years

  attr(out, "fleet") = fleet
  attr(out, "by") = by
  attr(out, "model") = model
  attr(out, "weights") = weights

  class(out) = c("empirical_selectivity", "SS_timexsize", "matrix")

  return(out)

}

# Internal functions ------------------------------------------------------


.mta = function(x, bins, years, var="Obs", FUN=mean, scale=1, by="length",
                mid = NULL, normalize=TRUE, bin.name = "Bin", thr=1e-5) {

  FUN = match.fun(FUN)
  out = matrix(0, nrow=length(years), ncol=length(bins))
  rownames(out) = years
  colnames(out) = bins

  if(is.null(x)) return(NA*out)

  msg = sprintf("Variable %s not found.", var)
  if(!(var %in% names(x))) stop(msg)

  xout = tapply(X = x[, var], INDEX=list(x$Yr, x[, bin.name], x$Seas), FUN=FUN)
  xout = apply(xout, 1:2, FUN = FUN, na.rm=TRUE)

  obins = as.numeric(colnames(xout))
  oyear = as.numeric(rownames(xout))

  rrows = match(oyear, years)
  rcols = match(obins, bins)

  out[na.omit(rrows), na.omit(rcols)] = xout[!is.na(rrows), !is.na(rcols)]
  out[is.na(out)] = 0

  if(!is.null(mid)) colnames(out) = mid

  if(normalize) {

    if(is.matrix(scale)) {
      # convert scale to proportion
      scale = scale/rowSums(scale, na.rm=TRUE)
      ind = which(colSums(scale>thr) == 0)
      scale[, ind] = NA
    }

    out = out/scale
    xmax = apply(out, 1, max, na.rm=TRUE)
    xmax[xmax==0] = 1
    out = out/xmax

    out[is.nan(out)] = 0
    out[, ind] = NA


  } else {
    out = out/scale
    out[is.nan(out)] = 0
    out[is.na(out)]  = 0
    names(dimnames(out)) = c("year", by)
    return(out)
  }

  names(dimnames(out)) = c("year", by)
  attr(out, "by") = by

  if(by %in% c("length", "age"))
    class(out) = c("empirical_selectivity", "SS_timexsize", "matrix")

  return(out)

}

.getES = function(object, by="length", use="B", thr=1e-5) {

  by = match.arg(by, choices = c("length", "age", "gsize"))

  if(by=="gsize") return(NULL)

  code = switch (by,
    length = "len",
    age    = "age",
    gsize  = "size")

  # sizedbase
  db  = sprintf("%sdbase", code)
  nat = sprintf("nat%s", code)
  years = seq(from=object$startyr, to=object$endyr, by=1)

  naged  = object[[db]]
  natage = object[[nat]]

  nsamp = c("Nsamp_adj", "Nsamp_in", "Nsamp", "N")
  nsamp = nsamp[ min(which(nsamp %in% names(naged)))]

  # fill and zeros
  w_catch = .mta(object$catch, bins=seq_along(object$FleetNames),
                 years=years, var="Obs", FUN=sum, by="fleet", normalize=FALSE, bin.name = "Fleet")
  w_nsamp = .mta(naged, bins=seq_along(object$FleetNames),
                 years=years, var=nsamp, FUN=mean, by="fleet", normalize=FALSE, bin.name = "Fleet")
  w_effN  = .mta(naged, bins=seq_along(object$FleetNames),
                 years=years, var="effN", FUN=mean, by="fleet", normalize=FALSE, bin.name = "Fleet")

  thisFleets = unique(naged$Fleet)

  popbin = as.numeric(colnames(natage)[-(1:12)])
  fshbin = if(by=="age") object$agebins else object$lbins

  if(all(is.na(fshbin))) fshbin = popbin

  bins = .checkBins(fshbin = fshbin, popbin = popbin)

  natage = natage[natage$'Beg/Mid'==use & natage$Era=="TIME", ]
  natage = natage[order(natage$Yr), ]

  nbase = expand.grid(Yr=years, Seas=seq_len(object$nseasons), Sex=c(1,2))
  nbase = merge(nbase, natage, all=TRUE)
  nbase = nbase[order(nbase$Yr, nbase$Seas, nbase$Sex), ]

  yr = sort(unique(natage$Yr))
  natage_1 = nbase[nbase$Sex==1, -(1:12)]
  natage_2 = nbase[nbase$Sex==2, -(1:12)]

  natage_1[is.na(natage_1)] = 0
  natage_2[is.na(natage_2)] = 0

  natage = as.matrix(natage_1) + as.matrix(natage_2)
  natage = rowsum(natage, group = nbase[nbase$Sex==1, ]$Yr)
  natage = natage %*% bins$conv # reshape to new bins

  output = lapply(split(naged, f = list(naged$Fleet, naged$Sex)),
                  FUN = .mta, bins=bins$bin, years=years, scale=natage, by=by,
                  mid = bins$mid, thr=thr)

  for(i in seq_along(output)) {
    # thisFleets has the "fleet code"
    attr(output[[i]], "fleet") = object$FleetNames[thisFleets[i]]
    attr(output[[i]], "weights") = data.frame(year  = years,
                                              catch = w_catch[, thisFleets[i]],
                                              Nsamp = w_nsamp[, thisFleets[i]],
                                              Neff  = w_effN[, thisFleets[i]])
    sex = as.numeric(strsplit(names(output)[i], "\\.")[[1]][2])
    attr(output[[i]], "model") = .getSelex(object, fleet = thisFleets[i], sex=sex , by=by)

  }

  return(output)

}

.checkBins = function(fshbin, popbin, k=1000) {

  n = max(20, 2*((max(popbin) - max(fshbin))/min(diff(fshbin))))
  find0 = seq_len(length(fshbin)+2*n) - n
  fmod = lm(bin ~ ind, data=data.frame(bin=fshbin, ind=seq_along(fshbin)))

  pind0 = seq_len(length(popbin)+1)
  pmod = lm(bin ~ ind, data=data.frame(bin=popbin, ind=seq_along(popbin)))

  # nbin extends fshbins
  nbin0 = round(as.numeric(predict(fmod, newdata = data.frame(ind=find0))),3)
  mid0 = as.numeric(stats::filter(nbin0, filter=c(0.5, 0.5), sides=1))
  ind1 = which(nbin0 >= min(popbin) & nbin0 <= max(popbin))
  nbin = nbin0[ind1]
  mid = mid0[ind1+1]
  top = nbin0[ind1+1]
  pbin = popbin[popbin >= min(nbin) & popbin <= max(nbin)]

  cond = identical(nbin, pbin)

  popbin2 = round(as.numeric(predict(pmod, newdata = data.frame(ind=pind0))),3)

  # if cond is FALSE, regrid popbin

  # this will deal with irregular bins!
  xout = .rebinning(x=popbin2, y=na.omit(c(nbin, max(top))), k=100)

  # mini = head(approx(x=popbin2, n=k*length(popbin2)-(k-1))$y, -1)
  # out = cut(mini, breaks=c(nbin, max(top)), labels=FALSE, right = FALSE)
  # levels = unique(na.omit(out))
  # out = matrix(out, ncol=k, byrow = TRUE)
  # .mytable = function(x, levels) table(factor(x, levels=levels))
  # xout = t(apply(out, 1, .mytable, levels=levels))/k

  output = list(bin=nbin, mid=mid, conv = xout)
  return(output)
}


.getSelex = function(object, fleet, sex=1, by="length") {

  obj = if(by=="length") object$sizeselex else object$ageselex
  Asel = ifelse("Asel" %in% names(table(object$ageselex$Factor)), "Asel", "Asel2")
  var = if(by=="length") "Lsel" else Asel
  n = if(by=="length") 5 else 7

  tt = obj[obj$Factor %in% var & obj$Fleet %in% fleet & obj$Sex %in% sex, ]
  syr = tt$Yr
  years = seq(from=object$startyr, to=object$endyr, by=1)

  tt = cbind(Yrf=cut(syr, breaks=syr, right = FALSE, include.lowest = TRUE), tt[, -(1:n)])

  db = data.frame(Yr=years, Yrf=cut(years, breaks=syr, right = FALSE, include.lowest = TRUE))

  db = merge(db, tt, sort=FALSE, all.y=TRUE)
  db = db[!is.na(db$Yr), , drop=FALSE]
  db = db[order(db$Yr), -c(1,2), drop=FALSE]

  out = list(x=as.numeric(colnames(db)), y=as.matrix(db), years=years)

  return(out)

}

.rebinning = function(x, y, k=100) {
  .mini = function(x, k=100) head(approx(x=x, n=k*length(x)-(k-1))$y, -1)
  .mytable = function(x, levels) table(factor(x, levels=levels))
  xm = cbind(head(x, -1), tail(x, -1))
  mini = apply(xm, 1, .mini, k=k)
  out = cut(mini, breaks=y, right = FALSE)
  levels = levels(out)
  out = matrix(out, ncol=k, byrow = TRUE)
  xout = t(apply(out, 1, .mytable, levels=levels))/k
  return(xout)
}

