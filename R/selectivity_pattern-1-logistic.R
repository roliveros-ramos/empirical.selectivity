
fit_selectivity_1 = function(object, prob=0.95, FUN, ...) {
  # create a list of model parameters

  FUN = match.fun(FUN)

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_1 = function(x, y, prob, ...) {

    x = as.numeric(x)
    y = as.numeric(y)
    par = .logistic1_guess(x, y, prob=prob)
    mod = list()
    mod$optim = optim(par, fn=.logistic1_fit, x=x, y=y, prob=prob, FUN=FUN, method="L-BFGS-B",
                lower=c(0.1, 0.1, 0.5), upper=c(Inf, Inf, 2))
    mod$par   = mod$optim$par
    mod$guess = par
    mod$predict = function(x) .logistic1(x, mod$par, prob=prob)
    mod$pattern = 1
    class(mod) = "empirical_selectivity_model"

    pred = mod$predict(x)
    scale = mod$par[3]

    output = list(fitted=pred, x=x, y=y, model=mod, npar=2, scale=scale)
    return(output)
  }

  x = as.numeric(colnames(object))

  xo = object*0
  yo = object*0
  out = vector("list", nrow(object))
  npar = vector("integer", nrow(object))
  scales = vector("double", nrow(object))

  labs = sprintf("year = %s", rownames(object))
  if(nrow(object)==1) labs = rownames(object)

  for(i in seq_len(nrow(object))) {
    message("Fitting ", labs[i], "\n")
    tmp = .fit_selectivity_1(x=x, y=object[i, ], prob=prob, ...)
    out[[i]]  = tmp$model
    xo[i, ]   = tmp$fitted
    yo[i, ]   = tmp$y
    npar[i]   = tmp$npar
    scales[i] = tmp$scale
  }

  fit = .compare_fit(xo, object, FUN)

  output = list(selectivity=xo, models=out, y=yo, x=x, pattern=rep(1, nrow(xo)),
                fit=fit, npar=npar, scale=scales)

  return(output)
}


# Internal functions ------------------------------------------------------

.logistic1 = function(x, par, prob=0.95) {

  p1 = par[1]
  p2 = par[2]
  p3 = par[3]

  const = 1/(1 - prob) - 1

  sel = 1/(1 + exp(-1*log(const)*(x-p1)/p2))

  return(sel)
}


.logistic1_guess = function(x, y, prob=0.95) {

  p1 = min(x[which(y >= 0.5)])
  p2 = min(x[which(y >= prob)]) - p1
  top = y[which(y >= 0.95)]
  p3 = mean(top, na.rm=TRUE)/max(top, na.rm=TRUE)

  if(is.nan(p3)) p3 = 1
  if(!is.finite(p3)) p3 = 1

  out = c(p1=p1, p2=p2, p3=p3)

  return(out)

}

.logistic1_fit = function(par, x, y, FUN, prob=0.95) {
  FUN = match.fun(FUN)
  removeZero = !is.finite(FUN(1,0))
  if(removeZero) {
    y = .tinyExp(y)
  }
  fit = par[3]*.logistic1(x, par, prob=prob)
  out = sum(FUN(fit, y), na.rm=TRUE)
  return(out)
}

#' @export
.SS_writeselec_1 = function(object, file=NULL, phase=2, t=1, ...) {
  # write a selectivity_model object into lines for a ctl file.
  # can be printed in the console or to a file.

  if(nrow(object$selectivity)>1)
    message("Using parameters of year", rownames(object$selectivity)[t])

  n = 2
  m = min(object$x)
  M = max(object$x)

  fleet_nm = attr(object, "fleet")

  names = c("#_", "LO", "HI", "INIT", "PRIOR", "PR_SD", "PR_type", "PHASE", "env-var",
            "use_dev", "dev_mnyr", "dev_mxyr", "dev_PH", "Block", "Blk_Fxn", "# parm_name")

  lo   = c(m, -(M-m))
  hi   = c(M, M-m)
  init = object$models[[t]]$par[1:2]
  init = pmax(pmin(init, hi), lo)
  init = round(init, 3)
  prior = c(round(0.5*(m+M)), 0)
  prior_sd = round(c(sqrt(((M-m)^2)/12), 0))
  prior_ty = rep(0, n)

  PHASE = rep(phase, n)

  env_var  = use_dev = dev_mnyr = dev_mxyr = 0
  dev_PH = 0.5
  Block = Blk_fxn = 0

  nm = c("#  Size_Logistic_inflection_%s(1)", "#  Size_Logistic_width95_%s(1)")
  nm = sprintf(nm, fleet_nm)

  out = cbind("", lo, hi, init, prior, prior_sd, prior_ty, PHASE,
              env_var, use_dev, dev_mnyr, dev_mxyr, dev_PH, Block, Blk_fxn, nm)

  ofile = sprintf("%s_logistic.ctl", fleet_nm)

  cat(sprintf("# %s Length Selex\n", fleet_nm), file=ofile)
  write(t(out), ncolumns = ncol(out), file=ofile, sep="\t", append = TRUE)

  return(invisible(out))

}
