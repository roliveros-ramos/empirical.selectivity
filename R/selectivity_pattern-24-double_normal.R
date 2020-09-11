
fit_selectivity_24 = function(object, FUN, ...) {
  # create a list of model parameters

  FUN = match.fun(FUN)

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_24 = function(x, y, ...) {

    x = as.numeric(x)
    y = as.numeric(y)
    par = .doubleNorm24_guess(x, y)
    mod = optim(par, fn=.doubleNorm24_fit, x=x, y=y, FUN=FUN, method="L-BFGS-B")
    mod$guess = par
    pred = .doubleNorm24(x, mod$par)

    output = list(fitted=pred, x=x, y=y, model=mod, npar=6)
    return(output)
  }

  x = as.numeric(colnames(object))

  xo = object*0
  out = vector("list", nrow(object))
  npar = vector("integer", nrow(object))

  labs = sprintf("year = %s", rownames(object))
  if(nrow(object)==1) labs = rownames(object)

  for(i in seq_len(nrow(object))) {
    message("Fitting ", labs[i], "\n")
    tmp = .fit_selectivity_24(x=x, y=object[i, ], ...)
    out[[i]] = tmp$model
    xo[i, ]  = tmp$fitted
    npar[i]  = tmp$npar
  }

  fit = .compare_fit(xo, object, FUN)

  output = list(selectivity=xo, models=out, y=object, x=x,
                pattern=rep(24, nrow(xo)), fit=fit, npar=npar)

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


.doubleNorm24_fit = function(par, x, y, FUN) {
  FUN = match.fun(FUN)
  fit = .doubleNorm24(x, par)
  out = sum(FUN(fit, y), na.rm=TRUE)
  return(out)
}

#' @export
.SS_writeselec_24 = function(object, file=NULL, phase=2, t=1, ...) {
  # write a selectivity_model object into lines for a ctl file.
  # can be printed in the console or to a file.

  if(nrow(object$selectivity)>1)
    message("Using parameters of year", rownames(object$selectivity)[t])

  n = 6
  m = min(object$x)
  M = max(object$x)

  fleet_nm = attr(object, "fleet")

  names = c("#_", "LO", "HI", "INIT", "PRIOR", "PR_SD", "PR_type", "PHASE", "env-var",
            "use_dev", "dev_mnyr", "dev_mxyr", "dev_PH", "Block", "Blk_Fxn", "# parm_name")

  lo   = c(m, rep(-15, n-1))
  hi   = c(M, rep(+15, n-1))
  init = object$models[[t]]$par
  init = pmax(pmin(init, hi), lo)
  init = round(init, 3)
  prior = c(round(0.5*(m+M)), rep(0, n-1))
  prior_sd = round(c(sqrt(((M-m)^2)/12), rep(0, n-1)))
  prior_ty = rep(0, n)

  PHASE = c(phase-1, rep(phase+1, 3), -rep(phase-1, 2))

  env_var  = use_dev = dev_mnyr = dev_mxyr = 0
  dev_PH = 0.5
  Block = Blk_fxn = 0

  nm = c("#  Size_DblN_peak_%s(1)", "#  Size_DblN_top_logit_%s(1)",
         "#  Size_DblN_ascend_se_%s(1)", "#  Size_DblN_descend_se_%s(1)",
         "#  Size_DblN_start_logit_%s(1)", "#  Size_DblN_end_logit_%s(1)")
  nm = sprintf(nm, fleet_nm)

  out = cbind("", lo, hi, init, prior, prior_sd, prior_ty, PHASE,
              env_var, use_dev, dev_mnyr, dev_mxyr, dev_PH, Block, Blk_fxn, nm)

  ofile = sprintf("%s_doubleNormal.ctl", fleet_nm)

  cat(sprintf("# %s Length Selex\n", fleet_nm), file=ofile)
  write(t(out), ncolumns = ncol(out), file=ofile, sep="\t", append = TRUE)

  return(invisible(out))

}
