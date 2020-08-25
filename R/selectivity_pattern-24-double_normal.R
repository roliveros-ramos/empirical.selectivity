
#' @export
fit_selectivity_24 = function(object, ...) {
  # create a list of model parameters

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_24 = function(x, y, ...) {

    x = as.numeric(x)
    y = as.numeric(y)
    par = .doubleNorm24_guess(x, y)
    mod = optim(par, fn=.doubleNorm24_fit, x=x, y=y, method="L-BFGS-B")
    pred = .doubleNorm24(x, mod$par)

    output = list(fitted=pred, x=x, y=y, model=mod)
    return(output)
  }

  x = as.numeric(colnames(object))

  xo = object*0
  out = vector("list", nrow(object))

  for(i in seq_len(nrow(object))) {
    cat("year =",i, "\n")
    tmp = .fit_selectivity_24(x=x, y=object[i, ], ...)
    out[[i]] = tmp$model
    xo[i, ]  = tmp$fitted
  }

  output = list(selectivity=xo, models=out, y=object, x=x)

  return(output)
}


# Internal functions ------------------------------------------------------

.doubleNorm24 = function(x, par) {

  p1 = par[1]
  p2 = par[2]
  p3 = par[3]
  p4 = par[4]
  p5 = par[5]
  p6 = par[6]

  .inv_logit = function(x) 1/(1 + exp(-x))
  .rescaling = function(x, y, ind) y + (1 - y)*(x - x[ind])/(1 - x[ind])

  startbin = ifelse(as.integer(p5) < -1000, floor(-1001 - p5), 1)
  lastbin  = ifelse(as.integer(p6) < -1000, floor(-1000 - p6), length(x))

  binwidth2 = mean(diff(x), na.rm=TRUE) # we don't have this one locally

  peak1     = p1
  peak2     = peak1 + binwidth2 + (0.99*x[lastbin] - peak1 - binwidth2)*.inv_logit(p2)
  upselex   = exp(p3) # ensure positivity of variance
  downselex = exp(p4) # ensure positivity of variance
  point1    = .inv_logit(p5)
  point2    = .inv_logit(p6)

  join1 = 1/(1 + exp(-(20/(1 + abs(x - peak1))) * (x - peak1)))
  join2 = 1/(1 + exp(-(20/(1 + abs(x - peak2))) * (x - peak2)))

  asc = exp(-(x - peak1)^2/upselex)
  if(p5 > -999) asc = .rescaling(x=asc, y=point1, ind=max(startbin, 1))

  dsc = exp(-(x - peak2)^2/downselex)
  if(p6 > -999) dsc = .rescaling(x=dsc, y=point2, ind=lastbin)

  sel = (1 - join1)*asc + join1*(1 - join2 + dsc*join2)

  if(as.integer(p5) < -1000) sel[seq_along(x) <= startbin] = 1e-6
  if(as.integer(p6) < -1000) sel[seq_along(x) > lastbin] = sel[lastbin]

  return(sel)
}

.doubleNorm24_guess = function(x, y) {

  logit = function(x) log(x/(1-x))
  logS = function(x, y, ind, peak) log(mean(1/(-log(y[ind])/(x[ind]-peak)^2)))

  binwidth2 = mean(diff(x))
  plateau = x[which(y>0.98*max(y))]
  peak1 = min(plateau)
  peak2 = max(plateau)
  if(peak1==peak2) peak2 = peak1 + binwidth2
  p1 = peak1
  p2 = logit((peak2 - peak1)/(0.99*max(x) - peak1 - binwidth2))
  p3 = logS(x, y, ind=which(x<peak1 & y > 1e-3), peak=peak1)
  p4 = logS(x, y, ind=which(x>peak2 & y > 1e-3), peak=peak2)
  p5 = logit(y[1])
  p6 = logit(y[length(x)])

  out = c(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6)

  return(out)

}


.doubleNorm24_fit = function(par, x, y) {
  fit = .doubleNorm24(x, par)
  out = sum((log(fit) - log(y))^2)
  return(out)
}


