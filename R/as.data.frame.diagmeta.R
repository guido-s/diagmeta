#' Extract data frame from diagmeta objects
#' 
#' @description
#' Extract data frame from objects of class \code{diagmeta}.
#' 
#' @param x An object of class \code{diagmeta}.
#' @param row.names Argument of R function \code{\link{as.data.frame}}
#'   (ignored).
#' @param optional Argument of R function \code{\link{as.data.frame}}
#'   (ignored).
#' @param \dots Other arguments.
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}} \code{\link{summary.diagmeta}}
#'
#' @return
#' A data frame is returned by the function \code{as.data.frame}.
#' 
#' @examples
#' # FENO dataset
#' #
#' data(Schneider2017)
#' 
#' diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
#'                   studlab = paste(author, year, group),
#'                   data = Schneider2017,
#'                   log.cutoff = TRUE)
#' 
#' as.data.frame(diag1)
#'
#' @method as.data.frame diagmeta
#' @export
#' @export as.data.frame.diagmeta


as.data.frame.diagmeta <- function(x, row.names=NULL, optional=FALSE, ...) {
  
  chkclass(x, "diagmeta")
  
  ## Remove element 'call' from object of class meta to get rid
  ## of an error message in meta-analyses with six studies:
  ## 'Error: evaluation nested too deeply: infinite recursion ...'
  ##
  ## NB: Element 'call' which is of length six contains information
  ##     on the function call.
  ##
  x$call <- NULL
  
  sel <- as.vector(lapply(x, length) == length(x$TP))
  
  res <- as.data.frame(x[names(x)[sel]], ...)
  
  attr(res, "version") <- packageDescription("diagmeta")$Version
  
  res
}
