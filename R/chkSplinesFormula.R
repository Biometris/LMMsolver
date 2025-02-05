chkSplinesFormula <- function(spline) {
  splErr <- paste("spline should be a formula of form \"~ spl1D() + ... + ",
                  "spl1D()\", \"~ spl2D()\" or \"~spl3D()\"\n")
  if (!is.null(spline)) {
    if (!inherits(spline, "formula")) stop(splErr)
    spline <- formula(paste((gsub(pattern = "LMMsolver::",
                                  replacement = "",
                                  as.character(spline))), collapse = ""))
    splTrms <- terms(spline, specials = c("spl1D", "spl2D", "spl3D","torus"))
    splSpec <- attr(splTrms, "specials")
    if (length(terms(splTrms)) != 2 ||
        ## Spline formula should consist of splxD() terms and nothing else.
        length(unlist(splSpec)) != length(labels(terms(spline))) ||
        length(splSpec$spl2D) > 1 ||
        length(splSpec$spl3D) > 1) {
      stop(splErr)
    }
  }
}
