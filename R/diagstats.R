#' Calculate statistical measures of test performance for objects of
#' class \code{diagmeta}
#' 
#' @description
#' The user can provide cutoffs, sensitivities, and / or specificities
#' to calculate the respective quantities (with confidence
#' intervals). Furthermore, positive predictive values (PPV), negative
#' predictive values (NPV), and probabilities of disease (PD) are
#' calculated if the prevalence is provided.
#' 
#' @param x An object of class \code{diagmeta}.
#' @param cutoff A numeric or vector with cutoff value(s).
#' @param sens A numeric or vector with sensitivity value(s).
#' @param spec A numeric or vector with specificity value(s).
#' @param prevalence A numeric or vector with the prevalence(s).
#' @param level The level used to calculate confidence intervals.
#' 
#' @return
#' A data frame of class "diagstats" with the following variables:
#' \item{cutoff}{Cutoffs provided in argument \code{cutoff} and / or
#'   model-based cutoff values for given sensitivities /
#'   specificities.}
#' \item{Sens}{Sensitivities provided in argument \code{sens} and / or
#'   model-based estimates of the sensitivity for given cutoffs /
#'   specificities.}
#' \item{lower.Sens, upper.Sens}{Lower and upper confidence limits of
#'   the sensitivities.}
#' \item{Spec}{Specificities provided in argument \code{spec} and / or
#'   model-based estimates of the specificity for given cutoffs /
#'   sensitivities.}
#' \item{lower.Spec, upper.Spec}{Lower and upper confidence limits of
#'   the specificities.}
#' \item{prevalence}{As defined above.}
#' \item{PPV}{Positive predictive value (based on the prevalence).}
#' \item{lower.PPV, upper.PPV}{Lower and upper confidence limits of
#'   positive predictive values.}
#' \item{NPV}{Negative predictive value (based on the prevalence)}
#' \item{lower.NPV, upper.NPV}{Lower and upper confidence limits of
#'   negative predictive values.}
#' \item{PD}{Probability of disease if the given cutoff value was
#'   observed as the measurement for an individual.}
#' \item{lower.PD, upper.PD}{Lower and upper confidence limits of
#'   probabilities of disease.}
#' \item{dens.nondiseased}{Value of the model-based density function at the
#'    cutoff(s) for non-diseased individuals.}
#' \item{dens.diseased}{Value of the model-based density function at the
#'    cutoff(s) for diseased individuals.}
#' 
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de},
#' Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}} \code{\link{print.diagstats}}
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
#' # Results at the optimal cutoff
#' #
#' diagstats(diag1)
#' 
#' # Results for cutoffs 25 and 50 (and a prevalence of 10%)
#' #
#' diagstats(diag1, c(25, 50), prevalence = 0.10)
#' 
#' # Results for sensitivity and specificity of 0.95
#' #
#' diagstats(diag1, sens = 0.95, spec = 0.95)
#' 
#' @export


diagstats <- function(x,
                      cutoff = x$optcut,
                      sens, spec,
                      prevalence,
                      level = 0.95) {
  
  
  chkclass(x, "diagmeta")
  ##
  cutoff.given <- !missing(cutoff)
  sens.given <- !missing(sens)
  spec.given <- !missing(spec)
  ##
  chknumeric(cutoff)
  ##
  if (!missing(prevalence))
    chklevel(prevalence, length = 0)
  else
    prevalence <- NA
  ##
  if (sens.given)
    chklevel(sens, length = 0)
  ##
  if (spec.given)
    chklevel(spec, length = 0)
  ##
  chklevel(level)
  #
  direction <- replaceNULL(x$direction, "increasing")
  #
  min.cutoff <- replaceNULL(x$min.cutoff, min(x$data.lmer$Cutoff, na.rm = TRUE))
  max.cutoff <- replaceNULL(x$max.cutoff, max(x$data.lmer$Cutoff, na.rm = TRUE))
  #
  log.cutoff <- x$log.cutoff
  
  regr <- x$regr
  ##
  alpha0 <- regr$alpha0
  alpha1 <- regr$alpha1
  beta0 <- regr$beta0
  beta1 <- regr$beta1
  #
  var.alpha0 <- regr$var.alpha0
  var.alpha1 <- regr$var.alpha1
  var.beta0 <- regr$var.beta0
  var.beta1 <- regr$var.beta1
  cov.alpha0.alpha1 <- regr$cov.alpha0.alpha1
  cov.alpha0.beta0 <- regr$cov.alpha0.beta0
  cov.alpha0.beta1 <- regr$cov.alpha0.beta1
  cov.alpha1.beta0 <- regr$cov.alpha1.beta0 
  cov.alpha1.beta1 <- regr$cov.alpha1.beta1
  cov.beta0.beta1 <- regr$cov.beta0.beta1  
  #
  distr <- x$distr
  
  
  if (sens.given & spec.given) {
    cutoff1 <- (qdiag(1 - sens, distr) - alpha1) / beta1
    cutoff2 <- (qdiag(spec, distr) - alpha0) / beta0
    #
    cutoff1 <-
      backtransf(cutoff1, direction, log.cutoff, min.cutoff, max.cutoff)
    cutoff2 <-
      backtransf(cutoff2, direction, log.cutoff, min.cutoff, max.cutoff)
    #
    if (cutoff.given)
      cutoff <- c(cutoff, cutoff1, cutoff2)
    else
      cutoff <- c(cutoff1, cutoff2)
  }
  else if (sens.given) {
    cutoff1 <- (qdiag(1 - sens, distr) - alpha1) / beta1
    #
    cutoff1 <-
      backtransf(cutoff1, direction, log.cutoff, min.cutoff, max.cutoff)
    #
    if (cutoff.given)
      cutoff <- c(cutoff, cutoff1)
    else
      cutoff <- cutoff1
  }
  else if (spec.given) {
    cutoff2 <- (qdiag(spec, distr) - alpha0) / beta0
    #
    cutoff2 <-
      backtransf(cutoff2, direction, log.cutoff, min.cutoff, max.cutoff)
    #
    if (cutoff.given)
      cutoff <- c(cutoff, cutoff2)
    else
      cutoff <- cutoff2
  }
  
  
  if (length(cutoff) != length(prevalence)) {
    if (length(cutoff) > length(prevalence))
      bad <- length(cutoff) %% length(prevalence) != 0
    else
      bad <- length(prevalence) %% length(cutoff) != 0
    ##
    if (bad) {
      if (cutoff.given)
        stop("Lengths of arguments 'cutoff' and 'prevalence' do not match.",
             call. = FALSE)
      else
        stop("Number of cutoffs and prevalences (argument 'prevalence') ",
             "do not match.")
    }
  }
  
  cutoff.tr <- transf(cutoff, direction, x$log.cutoff, min.cutoff, max.cutoff)
  
  Sens <- 1 - pdiag(beta1 * cutoff.tr + alpha1, distr)
  Spec <- pdiag(beta0 * cutoff.tr + alpha0, distr)
  ##
  seSens <- sqrt(regr$var.alpha1 +
                 cutoff.tr^2 * regr$var.beta1 +
                 2 * cutoff.tr * regr$cov.alpha1.beta1 +
                 x$var.diseased)
  ##
  seSpec <- sqrt(regr$var.alpha0 +
                 cutoff.tr^2 * regr$var.beta0 +
                 2 * cutoff.tr * regr$cov.alpha0.beta0 +
                 x$var.nondiseased)
  ##
  ci1 <- ci(beta1 * cutoff.tr + alpha1, seSens, level = level)
  ci0 <- ci(beta0 * cutoff.tr + alpha0, seSpec, level = level)
  ##
  lower.Sens <- 1 - pdiag(ci1$upper, distr)
  upper.Sens <- 1 - pdiag(ci1$lower, distr)
  ##
  lower.Spec <- pdiag(ci0$lower, distr)
  upper.Spec <- pdiag(ci0$upper, distr)
  
  
  ## PPV, NPV and probability of disease (depending on cutoff, because
  ## Sens and Spec depend on cutoff)
  ##
  dens.nondiseased <- beta0 * ddiag(beta0 * cutoff.tr + alpha0, distr)
  dens.diseased    <- beta1 * ddiag(beta1 * cutoff.tr + alpha1, distr)
  ##
  PPV <- Sens * prevalence /
    (prevalence * Sens + (1 - prevalence) * (1 - Spec))
  varPPV <-
    Spec^2 * (var.alpha0 + cutoff.tr^2 * var.beta0 +
                2 * cutoff.tr * cov.alpha0.beta0 + x$var.nondiseased) +
    (1 - Sens)^2 * (var.alpha1 + cutoff.tr^2 * var.beta1 +
                      2 * cutoff.tr * cov.alpha1.beta1 + x$var.diseased) -
    2 * Spec * (1 - Sens) *
    (cov.alpha0.alpha1 + cutoff.tr * cov.alpha0.beta1 +
       cutoff.tr * cov.alpha1.beta0 + cutoff.tr^2 * cov.beta0.beta1)
  ciPPV <- ci(qdiag(PPV, distr), sqrt(varPPV), level = level)
  lower.PPV <- pdiag(ciPPV$lower, distr)
  upper.PPV <- pdiag(ciPPV$upper, distr)  
  #
  NPV <- Spec * (1 - prevalence) /
    (prevalence * (1 - Sens) + (1 - prevalence) * Spec)
  varNPV <-
    (1 - Spec)^2 * (var.alpha0 + cutoff.tr^2 * var.beta0 +
                      2 * cutoff.tr * cov.alpha0.beta0 + x$var.nondiseased) +
    Sens^2 * (var.alpha1 + cutoff.tr^2 * var.beta1 +
                2 * cutoff.tr * cov.alpha1.beta1 + x$var.diseased) -
    2 * Sens * (1 - Spec) *
    (cov.alpha0.alpha1 + cutoff.tr * cov.alpha0.beta1 +
       cutoff.tr * cov.alpha1.beta0 + cutoff.tr^2 * cov.beta0.beta1)
  ciNPV <- ci(qdiag(NPV, distr), sqrt(varNPV), level = level)
  lower.NPV <- pdiag(ciNPV$lower, distr)
  upper.NPV <- pdiag(ciNPV$upper, distr)  
  #
  PD <- prevalence * dens.diseased /
    (prevalence * dens.diseased + (1 - prevalence) * dens.nondiseased)
  varPD <- 
    (1 - 2 * Spec)^2 *
    (var.alpha0 + cutoff.tr^2 * var.beta0 +
       2 * cutoff.tr * cov.alpha0.beta0 + x$var.nondiseased) +
    (2 * Sens - 1)^2 *
    (var.alpha1 + cutoff.tr^2 * var.beta1 +
       2 * cutoff.tr * cov.alpha1.beta1 + x$var.diseased) -
    2 * (2 * Sens - 1) * (1 - 2 * Spec) *
    (cov.alpha0.alpha1 + cutoff.tr * cov.alpha0.beta1 +
       cutoff.tr * cov.alpha1.beta0 + cutoff.tr^2 * cov.beta0.beta1)
  ciPD <- ci(qdiag(PD, distr), sqrt(varPD), level = level)
  lower.PD <- pdiag(ciPD$lower, distr)
  upper.PD <- pdiag(ciPD$upper, distr)  
  #
  DOR <- Sens * Spec / (1 - Sens) / (1 - Spec)
  varDOR <- 
    var.alpha0 + cutoff.tr^2 * var.beta0 +
    2 * cutoff.tr * cov.alpha0.beta0 + x$var.nondiseased +
    var.alpha1 + cutoff.tr^2 * var.beta1 +
    2 * cutoff.tr * cov.alpha1.beta1 + x$var.diseased - 
    2 * cov.alpha0.alpha1 -
    2 * cutoff.tr * cov.alpha0.beta1 - 
    2 * cutoff.tr * cov.alpha1.beta0 -
    2 * cutoff.tr^2 * cov.beta0.beta1
  ciDOR <- ci(log(DOR), sqrt(varDOR), level = level)
  lower.DOR <- exp(ciDOR$lower)
  upper.DOR <- exp(ciDOR$upper)
  #
  LRpos <- Sens / (1 - Spec)
  ci.LRpos <- ci(log(LRpos), sqrt(varPPV), level = level)
  lower.LRpos <- exp(ci.LRpos$lower)
  upper.LRpos <- exp(ci.LRpos$upper)
  #
  LRneg <- (1 - Sens) / Spec
  ci.LRneg <- ci(log(LRneg), sqrt(varNPV), level = level)
  lower.LRneg <- exp(ci.LRneg$lower)
  upper.LRneg <- exp(ci.LRneg$upper)
  
  
  res <- data.frame(cutoff,
                    Sens, lower.Sens, upper.Sens,
                    Spec, lower.Spec, upper.Spec,
                    LRpos, lower.LRpos, upper.LRpos,
                    LRneg, lower.LRneg, upper.LRneg,
                    prevalence,
                    PPV,lower.PPV, upper.PPV,
                    NPV, lower.NPV, upper.NPV,
                    PD, lower.PD, upper.PD,
                    dens.nondiseased, dens.diseased)
  #
  class(res) <- c("diagstats", "data.frame")
  #
  res
}
