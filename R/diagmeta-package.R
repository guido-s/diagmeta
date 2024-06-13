#' diagmeta: Brief overview of methods and general hints
#'
#' @description
#' R package \bold{diagmeta} implements the method by Steinhauser et
#' al. (2016) for the meta-analysis of diagnostic test accuracy studies
#' with multiple cutoffs.
#'
#' @name diagmeta-package
#'
#' @details
#' Main function of R package \bold{diagmeta} is the eponymous
#' \code{\link{diagmeta}}. Corresponding functions for printing and
#' plotting are available: \code{\link{print.diagmeta}},
#' \code{\link{plot.diagmeta}}
#'
#' Furthermore, a summary function and corresponding print function
#' are available to provide a briefer output:
#' \code{\link{summary.diagmeta}},
#' \code{\link{print.summary.diagmeta}}
#'
#' Additional functions provided in \bold{diagmeta} are
#' \code{\link{diagstats}} to calculate additional statistical
#' measures for the diagnostic test accuracy meta-analysis and
#' \code{\link{ipd2diag}} to transform individual participant data to
#' the data format required by \code{\link{diagmeta}}.
#' 
#' Type \code{help(package = "diagmeta")} for a listing of R functions and
#' datasets available in \bold{diagmeta}.
#' 
#' Type \code{citation("diagmeta")} for the preferred citation of 
#' R package \bold{diagmeta}.
#' 
#' To report problems and bugs
#' \itemize{
#' \item type \code{bug.report(package = "diagmeta")} if you do not
#'   use RStudio,
#' \item send an email to Guido Schwarzer
#'     \email{guido.schwarzer@@uniklinik-freiburg.de} if you use RStudio.
#' }
#'
#' The development version of \bold{diagmeta} is available on GitHub
#' https://github.com/guido-s/diagmeta.
#' 
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de},
#' Susanne Steinhauser \email{susanne.steinhauser@@uni-koeln.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de},
#' Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Steinhauser S, Schumacher M, Rücker G (2016):
#' Modelling multiple thresholds in meta-analysis of diagnostic test
#' accuracy studies.
#' \emph{BMC Medical Research Methodology},
#' \bold{16}, 97
#'
#' @keywords package
#'

"_PACKAGE"


NULL
