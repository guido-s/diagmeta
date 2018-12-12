#' Print method for diagmeta objects
#' 
#' Print method for objects of class \code{diagmeta}.
#' 
#' @param x An object of class \code{diagmeta}.
#' @param digits Number of significant digits for printing.
#' @param \dots Additional arguments.
#'
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Susanne
#'   Steinhauser \email{susanne.steinhauser@@uni-koeln.de}, Srinath
#'   Kolampally \email{kolampal@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}} \code{\link{summary.diagmeta}}
#' 
#' @examples
#' 
#' # FENO dataset
#' #
#' data(Schneider2017)
#' 
#' diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
#'                   studlab = paste(author, year, group),
#'                   data = Schneider2017,
#'                   model = "DIDS", log.cutoff = TRUE)
#' 
#' diag1
#'
#' @export print.diagmeta
#' @export
#' 
#' @importFrom stats quantile


print.diagmeta <- function(x, digits = 3, ...) {
  
  meta:::chkclass(x, "diagmeta")
  ##
  meta:::chknumeric(digits, min = 0, single = TRUE)
  
  cat("\nList and distribution of cutoffs:", "\nCutoffs")
  prmatrix(table(x$cutoff), collab = c("Frequency"))
  
  cat("\nNumber of cutoffs per study:\n")
  prmatrix(table(x$studlab), collab = c(" "))
  
  cat("\nQuantiles of the number of cutoffs in a study:\n")
  print(quantile(table(list(x$studlab))))
  
  cat("\nNumber of studies by number of cutoffs:")
  print(table(table(x$studlab)))
  
  cat("\nResults:\n")
  print(x$result.lmer)
  
  print(summary(x))
  
  invisible(NULL)
}
