#' Summary method for diagmeta
#' 
#' @description
#' Summary method for objects of class \code{diagmeta}.
#' 
#' @param object An object of class \code{diagmeta}.
#' @param \dots Additional arguments.
#'
#' @return
#' A list with classes 'summary.diagmeta' and 'diagmeta' is
#' returned. The list elements are identical to a
#' \code{\link{diagmeta}} object.
#' 
#' @author
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de},
#' Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}}
#' 
#' @examples
#' # FENO dataset
#' #
#' data(Schneider2017)
#' 
#' diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
#'                   studlab = paste(author, year, group),
#'                   data = Schneider2017,
#'                   model = "DIDS", log.cutoff = TRUE)
#' 
#' summary(diag1)
#'
#' @method summary diagmeta
#' @export
#' @export summary.diagmeta


summary.diagmeta <- function(object, ...) {
  
  meta:::chkclass(object, "diagmeta")
  
  res <- object
  ##
  class(res) <- c("summary.diagmeta", "diagmeta")
  
  res
}
