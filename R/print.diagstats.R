#' Print method for diagstats objects
#' 
#' @description
#' Print method for objects of class \code{diagstats}.
#' 
#' @param x An object of class \code{diagmeta}.
#' @param sensspec A logical indicating whether sensitivities and
#'   specificies should be printed.
#' @param predicted A logical indicating whether predicted values
#'   should be printed.
#' @param density A logical indicating whether values of the
#'   model-based density functions should be printed.
#' @param digits Number of significant digits for printing.
#' @param \dots Additional arguments.
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagstats}} \code{\link{diagmeta}}
#' 
#' @keywords print
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
#' # Values for prevalence 10% at cutoffs 25 and 50
#' #
#' ds1 <- diagstats(diag1, c(25, 50), 0.10)
#' ds1
#' print(ds1, predicted = FALSE)
#'
#' @method print diagstats
#' @export
#' @export print.diagstats


print.diagstats <- function(x,
                            sensspec = TRUE,
                            predicted = TRUE,
                            density = FALSE,
                            digits = 3, ...) {
  
  
  meta:::chkclass(x, "diagstats")
  ##
  meta:::chknumeric(digits)
  
  
  drop.names <- c()
  ##
  if (!sensspec)
    drop.names <- c("Sens", "seSens", "lower.Sens", "upper.Sens",
                    "Spec", "seSpec", "lower.Spec", "upper.Spec")
  ##
  if (all(is.na(x$prevalence)) | !predicted)
    drop.names <- c(drop.names, c("prevalence", "PPV", "NPV", "PD"))
  ##
  if (!density)
    drop.names <- c(drop.names, c("dens.nondiseased", "dens.diseased"))
  ##  
  x <- x[, !(names(x) %in% drop.names), drop = FALSE]
  
  
  if (ncol(x) != 0) {
    x <- round(x, digits)
    prmatrix(x, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x)))
  }
  else
    cat("No variables to print.\n")
  
  
  invisible(NULL)
}
