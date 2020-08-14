
#' @export
SS_output2 = function(..., thr=1e-5) {

  returndat = r4ss::SS_output(...)
  returndat$empirical_selectivity = empirical_selectivity(returndat, thr=thr)
  class(returndat) = c("SS_output", "list")
  return(returndat)

}
