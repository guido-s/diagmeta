#' Print method for summary of diagmeta objects
#' 
#' @description
#' Print method for objects of class \code{summary.diagmeta}.
#' 
#' @param x An object of class \code{summary.diagmeta}.
#' @param digits Number of significant digits for printing.
#' @param \dots Additional arguments.
#'
#' @author
#' Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de},
#' Susanne Steinhauser \email{susanne.steinhauser@@uni-koeln.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de},
#' Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}} \code{\link{summary.diagmeta}}
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
#'                   log.cutoff = TRUE)
#' 
#' summary(diag1)
#' print(summary(diag1), digits = 2)
#'
#' @method print summary.diagmeta
#' @export
#' @export print.summary.diagmeta


print.summary.diagmeta <- function(x, digits = 3, ...) {
  
  
  meta:::chkclass(x, "summary.diagmeta")
  ##
  meta:::chknumeric(digits, min = 0, length = 1)


  formatN <- meta:::formatN
  formatCI <- meta:::formatCI
  
  
  cat("\n*** Results of multiple cutoffs model ***\n")
  
  Tstudlab <- sum(table(table(x$studlab)))
  ##
  cat(paste0("\n", "Total number of studies: ", Tstudlab, "\n", sep = ""))
  
  Tcutoffs <- length(x$studlab)
  ##
  cat(paste0( "Total number of cutoffs: ", Tcutoffs, "\n", sep = ""))
  
  Nunicutoffs <- length(unique(x$cutoff))
  ##
  cat(paste0("Number of different cutoffs: ", Nunicutoffs, "\n\n", sep = ""))
  ##
  cat(paste0("Model: ", x$model, "\n\n", sep = ""))
  ##
  cat(paste0("Type of distribution: ", x$distr, "\n\n", sep = ""))
  ##
  cat(paste0("Cutoffs log transformed: ", x$log.cutoff, "\n\n", sep = ""))
  ##
  cat(paste0("The optimal cutoff value: ",
             formatN(round(x$optcut, digits), digits)))
  if (!is.na(x$lower.optcut))
    cat(paste0(" ",
               formatCI(formatN(round(x$lower.optcut, digits), digits),
                        formatN(round(x$upper.optcut, digits), digits))))
  cat("\n\n")
  ##
  cat("Sensitivity and specificity at optimal cutoff:\n")
  ##
  cat(paste0("\tSens: ",
             formatN(round(x$Sens.optcut, digits), digits),
             " ",
             formatCI(formatN(round(x$lower.Sens.optcut, digits), digits),
                      formatN(round(x$upper.Sens.optcut, digits), digits)),
             "\n", sep = ""))
  ##
  cat(paste0("\t", "Spec: ",
             formatN(round(x$Spec.optcut, digits), digits),
             " ",
             formatCI(formatN(round(x$lower.Spec.optcut, digits), digits),
                      formatN(round(x$upper.Spec.optcut, digits), digits)),
             "\n", sep = ""))
  
  
  cat("\nArea under the curve (AUC): \n")
  ##
  cat(paste0(" ",
             formatN(round(x$AUC, digits), digits),
             " ",
             formatCI(formatN(round(x$AUCSens.lower, digits), digits),
                      formatN(round(x$AUCSens.upper, digits), digits)),
             " - confidence region for sensitivity given specificity\n",
             sep = ""))
  ##
  cat(paste0(" ",
             formatN(round(x$AUC, digits), digits),
             " ",
             formatCI(formatN(round(x$AUCSpec.lower, digits), digits),
                      formatN(round(x$AUCSpec.upper, digits), digits)),
             " - confidence region for specificity given sensitivity\n",
             sep = ""))
  
  
  invisible(NULL)
}
