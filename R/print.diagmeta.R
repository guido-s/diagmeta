#' Print method for diagmeta objects
#' 
#' @description
#' Print method for objects of class \code{diagmeta}.
#' 
#' @param x An object of class \code{diagmeta}.
#' @param digits Number of significant digits for printing of optimal
#'   cutoff.
#' @param digits.prop Number of significant digits for proportions,
#'   e.g., sensitivities and specificities.
#' @param \dots Additional arguments.
#'
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de},
#' Susanne Steinhauser \email{susanne.steinhauser@@uni-koeln.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de}
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
#' diag1
#'
#' @method print diagmeta
#' @export
#' @export print.diagmeta
#' 
#' @importFrom stats quantile


print.diagmeta <- function(x,
                           digits = 3,
                           digits.prop = gs("digits.prop"),
                           ...) {
  
  
  chkclass(x, "diagmeta")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  #
  direction <- replaceNULL(x$direction, "increasing")
  #
  min.cutoff <- replaceNULL(x$min.cutoff, min(x$data.lmer$Cutoff, na.rm = TRUE))
  max.cutoff <- replaceNULL(x$max.cutoff, max(x$data.lmer$Cutoff, na.rm = TRUE))
  
  cat("*** Results of multiple cutoffs model ***\n")
  
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
  cat(paste0("Direction: ", direction, "\n\n", sep = ""))
  ##
  cat(paste0("Model: ", x$model, "\n\n", sep = ""))
  ##
  cat(paste0("Type of distribution: ", x$distr, "\n\n", sep = ""))
  ##
  cat(paste0("Cutoffs log transformed: ", x$log.cutoff, "\n\n", sep = ""))
  ##
  if (is.na(x$optcut))
    cat("The optimal cutoff iteration didn't converge.\n")
  else {
    cat(paste0("The optimal cutoff value: ",
               formatN(invert(x$optcut, direction, min.cutoff, max.cutoff),
                       digits)))
    if (!is.na(x$lower.optcut))
      cat(paste0(" ",
                 formatCI(
                   formatN(invert(x$lower.optcut, direction,
                                  min.cutoff, max.cutoff),
                           digits),
                   formatN(invert(x$upper.optcut, direction,
                                  min.cutoff, max.cutoff), digits))))
    cat("\n\n")
    ##
    cat("Sensitivity and specificity at optimal cutoff:\n")
    ##
    cat(paste0(" Sens: ",
               formatN(x$Sens.optcut, digits.prop),
               " ",
               formatCI(formatN(x$lower.Sens.optcut, digits.prop),
                        formatN(x$upper.Sens.optcut, digits.prop)),
               "\n", sep = ""))
    ##
    cat(paste0(" Spec: ",
               formatN(x$Spec.optcut, digits.prop),
               " ",
               formatCI(formatN(x$lower.Spec.optcut, digits.prop),
                        formatN(x$upper.Spec.optcut, digits.prop)),
               "\n", sep = ""))
  }
  
  cat("\nArea under the curve (AUC): \n")
  ##
  cat(paste0(" ",
             formatN(x$AUCSens, digits.prop),
             " ",
             formatCI(formatN(x$AUCSens.lower, digits.prop),
                      formatN(x$AUCSens.upper, digits.prop)),
             " - confidence region for sensitivity given specificity\n",
             sep = ""))
  ##
  cat(paste0(" ",
             formatN(x$AUCSpec, digits.prop),
             " ",
             formatCI(formatN(x$AUCSpec.lower, digits.prop),
                      formatN(x$AUCSpec.upper, digits.prop)),
             " - confidence region for specificity given sensitivity\n",
             sep = ""))
  
  
  invisible(NULL)
}
