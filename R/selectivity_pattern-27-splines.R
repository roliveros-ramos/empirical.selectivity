
fit_selectivity_27 = function(object, k=7, thr=1e-3, span=5, ...) {

  if(any(k<3)) {
    k = pmax(k, 3)
    warning("k must be greater or equal to 3, setting to 3.")
  }

  # main function to be applied to 'empirical_selectivity' object
  .fit_selectivity_27 = function(x, y, k, thr, span, tiny=1e-4) {
    # knots, values, derivatives at extremes
    # create a list of model parameters
    x = as.numeric(x)
    y = as.numeric(y)
    if(all(na.omit(y)==0)) {
      output = list(fitted=y, x=x, y=y, model=NULL)
      return(output)
    }
    ind = .nonNullPoints(y, thr=thr, span=span)
    x0 = x[ind]
    y0 = y[ind]
    mod = suppressMessages(fks(x0, log(y0), k = k-2, degree=3, prec=0))
    pred = predict(mod, newdata = data.frame(x=x), type="response")
    pred = exp(pred - max(pred, na.rm=TRUE))
    output = list(fitted=pred, x=x, y=y, model=mod)
    return(output)
  }

  x = as.numeric(colnames(object))

  xo = object*0
  out = vector("list", nrow(object))

  labs = sprintf("year = %s", rownames(object))
  if(nrow(object)==1) labs = rownames(object)

  for(i in seq_len(nrow(object))) {
    message("Fitting ", labs[i], "\n")
    tmp = .fit_selectivity_27(x=x, y=object[i, ], k=k, thr=thr, span=span, ...)
    out[[i]] = tmp$model
    xo[i, ]  = tmp$fitted
  }

  output = list(selectivity=xo, models=out, y=object, x=x, pattern=rep(27, nrow(xo)))

  return(output)
}


# Internal functions ------------------------------------------------------

.nonNullPoints = function(y, thr, span) {
  # copy y

  if(length(span)==1) span = c(span, span)
  if(any(span<0)) stop("span must be positive.")

  y[is.na(y)] = 0
  yx = cumsum(y)/sum(y)

  ind0 = which.min(yx<thr/2)
  ind1 = which.max(yx>=(1-thr/2))

  ind = c(ind0, ind1)
  ind = seq(from=ind[1]-span[1], to=ind[2]+span[2], by=1)
  ind = pmin(pmax(ind, 1), length(y))
  ind = sort(unique(ind))
  return(ind)
}


.SS_writeselec_27 = function(object, file=NULL, phase=2, fix_bounds=TRUE, t=1, ...) {
  # write a selectivity_model object into lines for a ctl file.
  # can be printed in the console or to a file.

  if(nrow(object$selectivity)>1)
    message("Using parameters of year", rownames(object$selectivity)[t])

  knots = object$models[[t]]$knots
  n = nrow(knots)
  m = min(object$x)
  M = max(object$x)

  fleet_nm = attr(object, "fleet")

  fixed = median(seq_len(n))
  if(fix_bounds) fixed = sort(unique(c(1, fixed, n)))
  fixed = fixed + 3 + n
  names = c("#_", "LO", "HI", "INIT", "PRIOR", "PR_SD", "PR_type", "PHASE", "env-var",
            "use_dev", "dev_mnyr", "dev_mxyr", "dev_PH", "Block", "Blk_Fxn", "# parm_name")


  lo   = c(0, -0.01, -1.00, rep(m, n), rep(-9, n))
  hi   = c(0, +1.00, +0.01, rep(M, n), rep(+7, n))
  init = c(0, knots$deriv[1], tail(knots$deriv,1),
           knots$knots, knots$value)
  init = pmax(pmin(init, hi), lo)
  init = round(init, 3)
  prior = c(0, 0, 0, knots$knots, rep(0, n))
  prior_sd = c(0, 1e-3, 1e-3, rep(0, n), rep(1, n))
  prior_ty = c(0, 1, 1, rep(0, n), rep(1, n))

  PHASE = c(-99, rep(phase+1, 2), rep(-99, n), rep(phase, n))

  PHASE[fixed] = -PHASE[fixed]

  env_var  = use_dev = dev_mnyr = dev_mxyr = 0
  dev_PH = 0.5
  Block = Blk_fxn = 0

  knot_nm = sprintf("# SizeSpline_Knot_%d_%%s(1)", seq_len(n))
  valu_nm = sprintf("# SizeSpline_Val_%d_%%s(1)", seq_len(n))

  nm = sprintf(c("# SizeSpline_Code_%s(1)", "# SizeSpline_GradLo_%s(1)",
                 "# SizeSpline_GradHi_%s(1)", knot_nm, valu_nm), fleet_nm)


  out = cbind("", lo, hi, init, prior, prior_sd, prior_ty, PHASE,
              env_var, use_dev, dev_mnyr, dev_mxyr, dev_PH, Block, Blk_fxn, nm)

  ofile = sprintf("%s_spline.ctl", fleet_nm)

  cat(sprintf("# %s Length Selex\n", fleet_nm), file=ofile)
  write(t(out), ncolumns = ncol(out), file=ofile, sep="\t", append = TRUE)

  return(invisible(out))

}
