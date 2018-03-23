#' Calculate statistical measures of test performance for objects of
#' class \code{diagmeta}
#' 
#' The user can provide cutoffs and prevalences to calculate
#' sensitivities and specificities (with confidence intervals),
#' positive predictive values (PPV), negative predictive values (NPV),
#' and probabilities of disease (PD) for a \code{diagmeta} object.
#' 
#' @param x An object of class \code{diagmeta}
#' @param cutoff A numeric or vector with cutoff value(s)
#' @param prevalence A numeric or vector with the prevalence(s)
#' @param level The level used to calculate confidence intervals
#' 
#' @return A data frame of class "diagstats" with the following
#'   variables:
#' 
#' \item{cutoff}{As defined above.}
#' \item{Sens}{Model-based estimate of the sensitivity for given
#'   cutoff}
#' \item{seSens}{Standard error of sensitivity}
#' \item{lower.Sens, upper.Sens}{Lower and upper confidence limits of
#'   the sensitivity}
#' \item{Spec}{Model-based estimate of the specificity for given
#'   cutoff}
#' \item{seSpec}{Standard error of sensitivity}
#' \item{lower.Spec, upper.Spec}{Lower and upper confidence limits of
#'   the specificity}
#' \item{prevalence}{As defined above.}
#' \item{PPV}{Positive predictive value, given the cutoff}
#' \item{NPV}{Negative predictive value, given the cutoff}
#' \item{PD}{Probability of disease if the given "cutoff" was observed
#'   as the measurement for an individual}
#' \item{dens.nondiseased}{Value of the model-based density function at the
#'   given cutoff for non-diseased individuals}
#' \item{dens.diseased}{Value of the model-based density function at the
#'   given cutoff for diseased individuals}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Srinath
#'   Kolampally \email{kolampal@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}} \code{\link{print.diagstats}}
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
#' # Values at the optimal cutoff
#' #
#' diagstats(diag1)
#' 
#' # Values for prevalence 10% at cutoffs 25 and 50
#' #
#' diagstats(diag1, c(25, 50), 0.10)
#' 
#' @export
#'
#' @importFrom meta ci


diagstats <- function(x,
                      cutoff = x$optcut,
                      prevalence = NA,
                      level = 0.95) {
  
  
  meta:::chkclass(x, "diagmeta")
  ##
  meta:::chknumeric(cutoff)
  ##
  if (!missing(prevalence))
    meta:::chklevel(prevalence, single = FALSE)
  ##
  meta:::chklevel(level)
  
  
  if (length(cutoff) != length(prevalence)) {
    if (length(cutoff) > length(prevalence))
      bad <- length(cutoff) %% length(prevalence) != 0
    else
      bad <- length(prevalence) %% length(cutoff) != 0
    if (bad)
      stop("Lengths of arguments 'cutoff' and 'prevalence' do not match.")
  }

  
  regr <- x$regr
  ##
  alpha0 <- regr$alpha0
  alpha1 <- regr$alpha1
  beta0 <- regr$beta0
  beta1 <- regr$beta1
  ##
  distr <- x$distr
  
  
  if (x$log.cutoff)
    cutoff <- log(cutoff)
  
  
  Sens <- 1 - pdiag(beta1 * cutoff + alpha1, distr)
  Spec <- pdiag(beta0 * cutoff + alpha0, distr)
  ##
  seSens <- sqrt(regr$var.alpha1 +
                 cutoff^2 * regr$var.beta1 +
                 2 * cutoff * regr$cov.alpha1.beta1 +
                 x$var.diseased)
  ##
  seSpec <- sqrt(regr$var.alpha0 +
                 cutoff^2 * regr$var.beta0 +
                 2 * cutoff * regr$cov.alpha0.beta0 +
                 x$var.nondiseased)
  ##
  ci1 <- ci(beta1 * cutoff + alpha1, seSens, level = level)
  ci0 <- ci(beta0 * cutoff + alpha0, seSpec, level = level)
  ##
  lower.Sens <- 1 - pdiag(ci1$upper, distr)
  upper.Sens <- 1 - pdiag(ci1$lower, distr)
  ##
  lower.Spec <- pdiag(ci0$lower, distr)
  upper.Spec <- pdiag(ci0$upper, distr)
  
  
  ## PPV, NPV and probability of disease (depending on cutoff, because
  ## Sens and Spec depend on cutoff)
  ##
  dens.nondiseased <- beta0 * ddiag(beta0 * cutoff + alpha0, distr)
  dens.diseased    <- beta1 * ddiag(beta1 * cutoff + alpha1, distr)
  ##
  PPV <- Sens * prevalence /
    (prevalence * Sens + (1 - prevalence) * (1 - Spec))
  ##
  NPV <- Spec * (1 - prevalence) /
    (prevalence * (1 - Sens) + (1 - prevalence) * Spec)
  ##
  PD <- prevalence * dens.diseased /
    (prevalence * dens.diseased + (1 - prevalence) * dens.nondiseased)
  
  
  if (x$log.cutoff)
    cutoff <- exp(cutoff)
  
  
  res <- data.frame(cutoff = cutoff,
                    Sens = Sens, seSens = seSens,
                    lower.Sens = lower.Sens, upper.Sens = upper.Sens,
                    Spec = Spec, seSpec = seSpec,
                    lower.Spec = lower.Spec, upper.Spec = upper.Spec,
                    prevalence = prevalence,
                    PPV = PPV, NPV = NPV, PD = PD,
                    dens.nondiseased = dens.nondiseased,
                    dens.diseased = dens.diseased)
  ##
  class(res) <- c("diagstats", "data.frame")
  ##
  res
}
