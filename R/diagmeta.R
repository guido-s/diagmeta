#' Meta-analysis of diagnostic test accuracy studies with the multiple
#' cutoffs model
#' 
#' Diagnostic tests may be based on an ordinal or continuous biomarker
#' or an ordinal score together with a cutoff. The decision whether
#' the target condition is present or not (positive or negative test
#' result) depends on whether the observed value is above or below the
#' cutoff. Here we assume that higher values of the biomarker indicate
#' a greater probability for the target condition (e.g., a disease);
#' otherwise, the poling must be changed by multiplying all values of
#' the biomarker by -1. Sensitivity and specificity of the test depend
#' on the chosen cutoff and vary with the cutoff. In meta-analysis of
#' diagnostic accuracy studies, results are often reported for
#' multiple cutoffs within a study, and the cutoffs may differ between
#' studies.  The multiple cutoffs model creates a link between the
#' range of cutoffs and the respective pairs of sensitivity and
#' specificity and thus allows identifying cutoffs at which the test
#' is likely to perform best (Steinhauser et al., 2016).
#' 
#' Each row of the data set provides at least a study label, a cutoff
#' and the numbers of true positives, false positives, true negatives
#' and false negatives. Different studies may contribute a varying
#' number of cutoffs, as well as different sets of cutoffs.
#' 
#' The multiple cutoffs model is a multi-level random effects
#' model. At the study level, for the group of patients without the
#' target condition (in short disease-free), the specificities at all
#' available cutoffs together provide an estimate of the cumulative
#' distribution function (cdf) of the test results within the
#' disease-free individuals. Likewise, for patients with the target
#' condition (in short diseased), via the observed sensitivities at
#' all observed cutoffs we obtain an estimate of the cdf of the test
#' results within the diseased patients. At the meta-analytic level,
#' the model fits the data for both groups and all available cutoffs
#' over all studies. Based on a parametric model, it provides
#' estimates of the two cdfs for the two groups across all studies,
#' accounting for the between-study heterogeneity and correlation
#' between groups.
#' 
#' Users have the choice between the normal (argument
#' \code{distr="normal"}) and the logistic distribution (argument
#' \code{distr="logistic"} which is the default). In addition, it is
#' possible to log-transform the cutoffs (argument \code{log.cutoff},
#' default is \code{FALSE}).
#' 
#' The cdf, transformed using the quantile function of the chosen
#' distribution, is modelled by one of eight mixed linear models
#' ("DIDS", "CIDS", "DICS", "CICS", "DS", "CS", "DI", "CI") as
#' described in Steinhauser et al. (2016).  The argument
#' \code{equalvar} indicates if the variances of the biomarker in both
#' groups are assumed to be equal (equalvar = TRUE) or unequal
#' (equalvar = FALSE).
#' 
#' The pooled sensitivity and specificity values can be obtained at
#' every cutoff; a multiple cutoffs summary ROC (sROC) naturally
#' follows while preserving cutoff information. The optimal cutoff is
#' defined as the cutoff where the maximum of a weighted sum of
#' sensitivity and specificity is obtained: lambda * sensitivity + (1
#' - lambda) * specificity. The 95\% confidence intervals of
#' sensitivities, specificities and the optimal cutoff are estimated
#' using the delta method (Steinhauser et al., 2016).
#' 
#' @param TP,FP,TN,FN Numeric vectors giving the number of true
#'   positives, false positives, true negatives and false negatives
#' @param cutoff A number vector indicating the cutoff values
#' @param studlab A numeric or a character vector with study labels
#' @param data An optional data frame containing the study information
#' @param distr A character indicating the distribution (see Details)
#' @param model A character indicating the model (see Details)
#' @param equalvar A logical indicating whether the variances of the
#'   biomarker in both groups are thought equal (see Details)
#' @param lambda A numeric between 0 and 1 indicating the weight of
#'   the sensitivity (such that specificity receives weight 1 -
#'   lambda)
#' @param log.cutoff A logical indicating whether the cutoffs should
#'   be log-transformed
#' @param method.weights A character indicating the method for
#'   weighting the studies: \code{invvar} (default) means inverse
#'   variance weighting, \code{size} means weighting by group sample
#'   size, \code{equal} means that all studies are equally weighted
#' @param incr A numeric between 0 and 1 that is added as a continuity
#'   correction
#' @param level A numeric indicating the significance level (1 -
#'   alpha) for tests (default is 0.95)
#' @param n.iter.max A numeric indicating the maximal number of fixed
#'   point iterations for finding the optimal cutoff
#' @param tol A numeric indicating the tolerance for convergence of
#'   the fixed point iteration
#' @param silent A logical indicating whether iterations should be
#'   suppressed
#' @param \dots additional arguments
#' 
#' @return An object of class "diagmeta" with corresponding print,
#'   summary, and plot function. The object is a list containing the
#'   following components
#'
#' \item{TP, FP, TN, FN}{As defined above.}
#' \item{cutoff, studlab}{As defined above.}
#' \item{Sens}{Sensitivity (original data).}
#' \item{Spec}{Specificity (original data).}
#' \item{distr, model, equalvar, lambda}{As defined above.}
#' \item{log.cutoff, method.weights}{As defined above.}
#' \item{level, incr}{As defined above.}
#' \item{k}{The number of studies in the meta-analysis.}
#' \item{optcut}{The optimal cutoff.}
#' \item{Sens.optcut}{The sensitivity at the optimal cutoff.}
#' \item{lower.Sens.optcut, upper.Sens.optcut}{Corresponding lower and
#'   upper confidence limits.}
#' \item{Spec.optcut}{The specificity at the optimal cutoff.}
#' \item{lower.Spec.optcut, upper.Spec.optcut}{Corresponding lower and
#'   upper confidence limits.}
#' \item{AUC}{Area under the curve (AUC)}
#' \item{AUCSens.lower, AUCSens.upper}{Corresponding lower and upper
#'   confidence limits (based on the confidence region for the
#'   sensitivity, given the specificity)}
#' \item{AUCSpec.lower, AUCSpec.upper}{Corresponding lower and upper
#'   confidence limits (based on the confidence region for the
#'   specificity, given the sensitivity)}
#' \item{var.diseased, var.nondiseased}{The within-study variance for
#'   the diseased and non-diseased group, respectively.}
#' \item{AIC}{The value of the Akaike information criterion of the
#'   lmer object.}
#' \item{BIC}{The value of the Bayesian information criterion of the
#'   lmer object.}
#' \item{data.lmer}{A list with elements Study (study labels), Group
#'   (group labels (0 or 1)), Cutoff, N (group sizes), Negative
#'   (number of negative test results), NN (frequencies of negative
#'   test results).}
#' \item{result.lmer}{An object of class \code{\link[lme4]{lmer}}.}
#'
#' \item{weights}{??????.}
#' 
#' \item{regr}{A list containing ??????.}
#' 
#' \item{dist}{A list containing ??????.}
#' 
#' \item{n.iter.max}{As defined above.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{diagmeta} used to create
#'   object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Susanne
#'   Steinhauser \email{susanne.steinhauser@@uni-koeln.de}, Srinath
#'   Kolampally \email{kolampal@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{plot.diagmeta}, \link{summary.diagmeta}}
#' 
#' @references Steinhauser S, Schumacher M, Rücker G (2016), Modelling
#'   multiple thresholds in meta-analysis of diagnostic test accuracy
#'   studies. \emph{BMC Med Res Methodol}, \bold{16}, 97.
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
#' summary(diag1)
#' plot(diag1)
#' 
#' @export
#'
#' @importFrom lme4 lmer
#' @importFrom stats AIC BIC coef cov logLik vcov
#' @importFrom utils packageDescription


diagmeta <- function(TP, FP, TN, FN, cutoff, studlab, data = NULL,
                     ##
                     distr = "logistic", model = "DICS", equalvar = FALSE,
                     lambda = 0.5, log.cutoff = FALSE,
                     method.weights = "invvar",
                     ##
                     level = 0.95, incr = 0.5,
                     ##
                     n.iter.max = 1000, tol = 1e-08, silent = TRUE,
                     ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chklength <- meta:::chklength
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknull <- meta:::chknull
  chknumeric <- meta:::chknumeric
  setchar <- meta:::setchar
  ##
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  is.logistic <- distr == "logistic"
  is.normal <- distr == "normal"
  ##
  model <- setchar(model,
                   c("CI", "DI", "CS", "DS", "CICS", "DICS", "CIDS", "DIDS"))
  ##
  chklogical(equalvar)
  chklevel(level)
  chklogical(log.cutoff)
  ##
  chknumeric(incr, min = 0, single = TRUE)
  ##
  method.weights <- setchar(method.weights,
                            c("equal", "size", "invvar"))
  ##
  chknumeric(lambda, min = 0, single = TRUE)
  chknumeric(n.iter.max, min = 0, single = TRUE)
  chknumeric(tol, min = 0, single = TRUE)
  ##
  chklogical(silent)
  ##
  ## Additional arguments / checks
  ##
  fun <- "diagmeta"
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'TP', 'FP', 'TN', 'FN', 'cutoff', and 'studlab'
  ##
  TP <- eval(mf[[match("TP", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(TP)
  k.All <- length(TP)
  ##
  FP <- eval(mf[[match("FP", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(FP)
  ##
  TN <- eval(mf[[match("TN", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(TN)
  ##
  FN <- eval(mf[[match("FN", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(FN)
  ##
  cutoff <- eval(mf[[match("cutoff", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(cutoff)
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(studlab)
  
  
  ##
  ##
  ## (3) Check of essential variables
  ##
  ##
  chklength(FP, k.All, fun, name = "TP")
  chklength(TN, k.All, fun, name = "TP")
  chklength(FN, k.All, fun, name = "TP")
  chklength(cutoff, k.All, fun, name = "TP")
  chklength(studlab, k.All, fun, name = "TP")
  ##
  if (length(unique(cutoff)) == 1)
    stop("Model cannot be use with a single cutoff. Consider using, e.g., madad() from R package mada.")
  
  
  ##
  ##
  ## (4) Auxiliary function
  ##     (to calculate weighted cut-off point of two logistic
  ##      distributions by an iterative fixpoint procedure)
  ##
  ##
  g <- function(x)
    mean1 - sd1 * acosh(lambda / (1 - lambda) * sd0 / sd1 *
                        (1 + cosh((x - mean0) / sd0)) - 1)
  ##
  ## Inverse function
  ##
  f <- function(x)
    mean0 + sd0 * acosh((1 - lambda) / lambda * sd1 / sd0 *
                        (1 + cosh((x - mean1) / sd1)) - 1)
  ##
  ## Error handling for iterate() - chooses f or g for iterations
  ##
  saveIterate <- function(x0, n.iter.max, tol, silent) {
    
    tryCatch({ # try iterating with function f
      x <- iterate(f, x0, n.iter.max, tol, !silent)
      if (!silent)
        cat("* Optimal cut-off iteration with f *\n")
      return(list(x = x, iter = "f"))
    },
    warning = function(w) {
      suppressWarnings(x <- iterate(f, x0, n.iter.max, tol, !silent))
      warning(w$message)
      return(list(x = x, iter = "f"))
    },
    ## if error occurs, iterate with function g
    error = function(e) {
      tryCatch({ # try iteration with function g
        x <- iterate(g, x0, n.iter.max, tol, !silent)
        if (!silent)
          cat("* Optimal cut-off iteration with g *\n")
        return(list(x = x, iter = "g"))
      },
      warning = function(wa) {
        suppressWarnings(x <- iterate(g, x0, n.iter.max, tol, !silent))
        warning(w$message)
        return(list(x = x, iter = "g"))
      },
      error = function(er) {
        stop("Optimal cutoff iteration didn't converge. Use argument distribution = \"normal\".")
        return(list(x = NULL, iter = NULL))
      })
    }
    )
  }
  
  
  ##
  ##
  ## (5) Assignments
  ##
  ##
  k <- length(unique(studlab))
  ##
  N0 <- FP + TN   # number of non-diseased patients
  N1 <- TP + FN   # number of diseased patients
  N  <- c(N0, N1) # number of all patients in one study IN ONE GROUP!!!(0 or 1)
  ##
  NN <- (c(TN, FN) + incr) / (N + 2 * incr)
  ##
  ## Inverse variance weights for variance within studies
  ##
  if (is.logistic) {
    w0.iv <- (TN + incr) * (FP + incr) / (N0 + 2 * incr)
    w1.iv <- (TP + incr) * (FN + incr) / (N1 + 2 * incr)
  }
  ##
  else if (is.normal) {
    w0.iv <- (N0 + 2 * incr)^3 * dnorm(qnorm((TN + incr) /
                                             (N0 + 2 * incr)))^2 / ((TN + incr) *
                                                                    (FP + incr))
    w1.iv <- (N1 + 2 * incr)^3 * dnorm(qnorm((TP + incr) /
                                             (N1 + 2 * incr)))^2 / ((TP + incr) *
                                                                    (FN + incr))
  }
  ##
  ## Variance component within studies for both models
  ##
  var.nondiseased <- k / sum(w0.iv)
  var.diseased <- k / sum(w1.iv)
  ##
  ## Weights
  ##
  if (method.weights == "equal")
    w <- rep(1, length(c(N0, N1)))
  else if (method.weights == "size")
    w <- (length(N) * N) / sum(N)
  else if (method.weights == "invvar") {
    w <- c(w0.iv, w1.iv)
    ## scaling
    w <-  length(w) * w / sum(w)
  }
  
  
  iter <- ""
  
  
  ##
  ##
  ## (6) Model fitting
  ##
  ##
  ## Data frame consisting of rows with data for each cutoff of each
  ## study, first for all non-diseased individuals and then everything
  ## again for the diseased individuals (each cutoff of each study is
  ## named twice)
  ##
  Group <- c(rep(0, length(studlab)), rep(1, length(studlab)))
  ##
  if (log.cutoff)
    Cutoff <- log(c(cutoff, cutoff))
  else
    Cutoff <- c(cutoff, cutoff)
  ##
  Study <- c(studlab, studlab)
  ##
  if (equalvar) {
    if (model == "CI")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (1 | Study),
                   weights = w, ...)
    else if (model == "DI")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (1 + Group | Study),
                   weights = w, ...)
    else if (model == "CS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (0 + Cutoff | Study),
                   weights = w, ...)
    else if (model == "DS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (0 + Cutoff + Group:Cutoff | Study),
                   weights = w, ...)
    else if (model == "CICS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (Cutoff | Study),
                   weights = w, ...)
    else if (model == "DICS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (Cutoff + Group | Study),
                   weights = w, ...)
    else if (model == "CIDS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (Cutoff +  Group:Cutoff | Study),
                   weights = w, ...)
    else if (model == "DIDS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group + Cutoff + (Group * Cutoff | Study),
                   weights = w, ...)
  }
  else {
    if (model == "CI")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (1 | Study),
                   weights = w, ...)
    else if (model == "DI")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (1 + Group | Study),
                   weights = w, ...)
    else if (model == "CS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (0 + Cutoff | Study),
                   weights = w, ...)
    else if (model == "DS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (0 + Cutoff + Group:Cutoff | Study),
                   weights = w, ...)
    else if (model == "CICS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (Cutoff | Study),
                   weights = w, ...)
    else if (model == "DICS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (Cutoff + Group | Study),
                   weights = w, ...)
    else if (model == "CIDS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (Cutoff +  Group:Cutoff | Study),
                   weights = w, ...)
    else if (model == "DIDS")
      lme1 <- lmer(qdiag(NN, distr) ~ Group * Cutoff + (Group * Cutoff | Study),
                   weights = w, ...)
  }
  ##
  slme1 <- summary(lme1)
  ##
  ## Fixed effects
  ##
  cf <- coef(slme1)
  vc <- vcov(slme1)
  ##
  ## Random effects: variance-covariance matrix
  ##
  V <- cov(slme1$varcor$Study)
  ##
  ## Extract regression coefficients:
  ## - alpha0, beta0 (non-diseased)
  ## - alpha1 and beta1 (diseased)
  ##
  alpha0 <- cf[1, 1]
  alpha1 <- alpha0 + cf[2, 1]
  beta0 <- cf[3, 1]
  ##
  var.alpha0 <- vc[1, 1]
  var.alpha1 <- var.alpha0 + vc[2, 2] + 2 * vc[1, 2]
  var.beta0 <- vc[3, 3]
  ##
  cov.alpha0.beta0 <- vc[1, 3]
  cov.alpha0.alpha1 <- var.alpha0 + vc[1, 2]
  ##
  if (equalvar) {
    beta1 <- beta0
    var.beta1 <- var.beta0
    ##
    cov.alpha0.beta1 <- cov.alpha0.beta0
    cov.alpha1.beta0 <- cov.alpha0.beta0 + vc[2, 3]
    cov.alpha1.beta1 <- cov.alpha1.beta0
    ##
    cov.beta0.beta1 <- var.beta0
  }
  else {
    beta1 <- beta0 + cf[4, 1]
    var.beta1 <- var.beta0 + vc[4, 4] + 2 * vc[3, 4]
    ##
    cov.alpha0.beta1 <- cov.alpha0.beta0 + vc[1, 4]
    cov.alpha1.beta0 <- cov.alpha0.beta0 + vc[2, 3]
    cov.alpha1.beta1 <- cov.alpha0.beta1 + vc[2, 3] + vc[2, 4]
    ##
    cov.beta0.beta1 <- var.beta0 + vc[3, 4]
  }
  ##
  ## Correlation of data must be positive
  ##
  if (beta0 <= 0 | beta1 <= 0)
    stop("Regression yields negative correlation. Try another model or get better data. :)")
  
  
  ##
  ##
  ## (7) Compute parameters of the biomarker distributions and their
  ##     variances
  ##
  ##
  mean0 <- - alpha0 / beta0  # Mean disease-free
  mean1 <- - alpha1 / beta1  # Mean diseased
  ##
  sd0 <- 1 / beta0           # Standard deviation disease-free
  sd1 <- 1 / beta1           # Standard deviation diseased
  ##
  var.mean0 <- (alpha0^2) / (beta0^4) * var.beta0 + var.alpha0 /
    (beta0^2) - 2 * alpha0 / (beta0^3) * cov.alpha0.beta0
  var.mean1 <- (alpha1^2) / (beta1^4) * var.beta1 + var.alpha1 /
    (beta1^2) - 2 * alpha1 / (beta1^3) * cov.alpha1.beta1
  ##
  var.sd0 <- var.beta0 / (beta0^4)
  var.sd1 <- var.beta1 / (beta1^4)
  ##
  if (mean1 < mean0)
    stop("Estimated distribution of diseased patients is left of ",
         "non-diseased ones. Check if higher values for biomarker ",
         "really indicate illness.")
  
  
  ##
  ##
  ## (8) Distributions
  ##
  ##
  if (is.logistic) {
    ## Cutoffs of two logistics, weighted with lambda and 1 - lambda
    wmean <- (1 - lambda) * mean0 + lambda * mean1
    ##
    if ((1 - lambda) * sd1 != lambda * sd0) {
      x0 <- wmean
      iterateResult <- saveIterate(x0, n.iter.max, tol, silent)
      iter <- iterateResult$iter
      optcut <- iterateResult$x
    }
    else
      optcut <- wmean
  }
  else if (is.normal) {
    ## Cutoffs of two normals, weighted with lambda and 1 - lambda
    turn <- (mean0 * sd1^2 - mean1 * sd0^2) / (sd1^2 - sd0^2)
    rad <- sqrt(sd0^2 * sd1^2 * (2 * (sd1^2 - sd0^2) *
                                 (log(sd1) - log(sd0) - qlogis(lambda)) +
                                 (mean1 - mean0)^2) / (sd1^2 - sd0^2)^2)
    x0 <- turn - rad
    x1 <- turn + rad
    ##
    if (sd0 < sd1)
      optcut <- x1
    else if (sd0 > sd1)
      optcut <- x0
    else
      optcut <- (-qlogis(lambda) * sd0^2 - 0.5 * (mean0^2 - mean1^2)) /
        (mean1 - mean0)
    ##
    if (sd1 != sd0) {
      ## Derivations of optimal cutoff function
      S <- sqrt(2 * (beta0^2 - beta1^2) * (log(beta0 / beta1) - qlogis(lambda)) +
                (alpha0 * beta1 - alpha1 * beta0)^2)
      ##
      dalpha0 <- (-beta0 + beta1 / S * (alpha0 * beta1 - alpha1 * beta0)) /
        (beta0^2 - beta1^2)
      ##
      dalpha1 <- (beta1 - beta0 / S * (alpha0 * beta1 - alpha1 * beta0)) /
        (beta0^2 - beta1^2)
      ##
      dbeta0 <- (- alpha0 + 1 / (beta0 * S) * (beta0^2 - beta1^2)) /
        (beta0^2 - beta1^2)+
        (4 * beta0 * (log(beta0 / beta1) - qlogis(lambda)) -
         2 * alpha1 * (alpha0 * beta1 - alpha1 * beta0)) /
        (2 * S * (beta0^2 - beta1^2)) -
        (2 * beta0 * (alpha1 * beta1 - alpha0 * beta0 + S)) /
        ((beta0^2 - beta1^2)^2)
      ##
      dbeta1 <- (alpha1 - 1 / (beta1 * S) * (beta0^2 - beta1^2)) /
        (beta0^2 - beta1^2) +
        (-4 * beta1 * (log(beta0 / beta1) - qlogis(lambda)) + 2 * alpha0 *
         (alpha0 * beta1 - alpha1 * beta0)) / (2 * S * (beta0^2 - beta1^2)) +
        (2 * beta1 * (alpha1 * beta1 - alpha0 * beta0 + S)) /
        ((beta0^2 - beta1^2)^2)
      ##
      ## Variance estimate of optimal cutoff
      ##
      var.optcut <- dalpha0^2 * var.alpha0 + dalpha1^2 * var.alpha1 +
                                                       dbeta0^2 * var.beta0 + dbeta1^2 * var.beta1 +
                                                                                       2 * dalpha0 * dalpha1 * cov.alpha0.alpha1 + 2 * dalpha0 * dbeta0 * cov.alpha0.beta0 +
                                                                                       2 * dalpha0 * dbeta1 * cov.alpha0.beta1 + 2 * dalpha1 * dbeta0 * cov.alpha1.beta0 +
                                                                                       2 * dalpha1 * dbeta1 * cov.alpha1.beta1 + 2 * dbeta0 * dbeta1 * cov.beta0.beta1
    }
    else {
      ##
      ## Derivations of optimal cutoff function
      ##
      S <- qlogis(lambda) + 0.5 * (alpha1^2 - alpha0^2)
      ##
      dalpha0 <- (alpha1 * (alpha1 - alpha0) - S / beta0) /
        ((alpha1 - alpha0)^2)
      ##
      dalpha1 <- (-alpha0 * (alpha1 - alpha0) + S / beta0) /
        ((alpha1 - alpha0)^2)
      ##
      dbeta0 <- (-S / beta0^2) / (alpha1 - alpha0)
      ##
      ## Variance estimate of optimal cutoff
      ##
      var.optcut <- dalpha0^2 * var.alpha0 + dalpha1^2 * var.alpha1 + dbeta0^2 * var.beta0 +
                                                                               2 * dalpha0 * dalpha1 * cov.alpha0.alpha1 + 2 * dalpha0 * dbeta0 * cov.alpha0.beta0 +
                                                                               2 * dalpha1 * dbeta0 * cov.alpha1.beta0
    }
  }
  
  
  ##
  ##
  ## (9) Calculate sensitivity and specificity at optimal cutpoint
  ##
  ##
  ci.regr1 <- ciRegr(optcut,
                     alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                     var.diseased,
                     level)
  ##
  ci.regr0 <- ciRegr(optcut,
                     alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                     var.nondiseased,
                     level)
  ##
  ## Calculate sensitivity and specificity at optimal cutpoint
  ##
  Se <- calcSens(ci.regr1$TE, distr)
  lower.Se <- calcSens(ci.regr1$upper, distr)  
  upper.Se <- calcSens(ci.regr1$lower, distr)  
  ##
  Sp <- calcSpec(ci.regr0$TE, distr)
  lower.Sp <- calcSpec(ci.regr0$lower, distr)  
  upper.Sp <- calcSpec(ci.regr0$upper, distr) 
  
  
  ##
  ## AUC
  ##
  index <- function(a, b, n) a + (b - a) * (1:(n - 1))/n
  t <- index((qdiag(0.000001, distr) - alpha1)/beta1,
  (qdiag(0.999999, distr) - alpha0)/beta0, 1000)
  
  t0 <- ciRegr(t, alpha0, var.alpha0, beta0, var.beta0,
               cov.alpha0.beta0, var.nondiseased, level)
  t1 <- ciRegr(t, alpha1, var.alpha1, beta1, var.beta1,
               cov.alpha1.beta1,var.diseased, level)
  Sens <- calcSens(t1$TE, distr)
  ## Note: ciRegr refers to qdiag(1 - Sens), therefore CI limits must
  ##       be exchanged to achieve upperSens > Sens > lowerSens
  upperSens <- calcSens(t1$lower, distr)
  lowerSens <- calcSens(t1$upper, distr)
  ## calcSpec calculates spec, given a cutoff (lowerSpec < Spec <
  ## upperSpec)
  Spec <- calcSpec(t0$TE, distr)
  lowerSpec <- calcSpec(t0$lower, distr)
  upperSpec <- calcSpec(t0$upper, distr)
  AUC <- trapz(Spec, Sens)
  AUCSens.lower <- trapz(Spec, lowerSens)
  AUCSens.upper <- trapz(Spec, upperSens)
  AUCSpec.lower <- trapz(lowerSpec, Sens)
  AUCSpec.upper <- trapz(upperSpec, Sens)
  
  
  ## ciRegr calculates confidence intervals for qdiag(spec) and qdiag(1 - sens),
  
  ##
  ##
  ## (10) List with results
  ##
  ##
  res <- list(TP = TP, FP = FP, TN = TN, FN = FN,
              cutoff = cutoff, studlab = studlab,
              ##
              Sens = 1 - (FN + incr) / (N1 + 2 * incr),
              Spec = (TN + incr) / (N0 + 2 * incr),
              ##
              distr = distr, model = model, equalvar = equalvar,
              lambda = lambda,
              ##
              level = level, log.cutoff = log.cutoff,
              incr = incr, method.weights = method.weights,
              ##
              k = k,
              ##
              optcut = if (log.cutoff) exp(optcut) else optcut,
              ##
              Sens.optcut = Se,
              lower.Sens.optcut = lower.Se,
              upper.Sens.optcut = upper.Se,
              Spec.optcut = Sp,
              lower.Spec.optcut = lower.Sp,
              upper.Spec.optcut = upper.Sp,
              AUC =  AUC,
              AUCSens.lower = AUCSens.lower,
              AUCSens.upper = AUCSens.upper,
              AUCSpec.lower = AUCSpec.lower,
              AUCSpec.upper = AUCSpec.upper,
              ##
              var.diseased = var.diseased,
              var.nondiseased = var.nondiseased,
              ##
              AIC = AIC(logLik(lme1)),
              BIC = BIC(lme1),
              ##
              data.lmer = list(Study = Study, Group = Group, Cutoff = Cutoff,
                               N = N, Negative = c(TN, FN), NN = NN),
              result.lmer = lme1,
              weights = w,
              ##
              regr = list(alpha0 = alpha0, var.alpha0 = var.alpha0,
                          beta0 = beta0, var.beta0 = var.beta0,
                          cov.alpha0.beta0 = cov.alpha0.beta0,
                          alpha1 = alpha1, var.alpha1 = var.alpha1,
                          beta1 = beta1, var.beta1 = var.beta1,
                          cov.alpha1.beta1 = cov.alpha1.beta1,
                          cov.alpha0.alpha1 = cov.alpha0.alpha1,
                          cov.alpha0.beta1 = cov.alpha0.beta1,
                          cov.alpha1.beta0 = cov.alpha1.beta0,
                          cov.beta0.beta1 = cov.beta0.beta1),
              ##
              dist = list(mean0 = mean0, var.mean0 = var.mean0,
                          sd0 = sd0, var.sd0 = var.sd0,
                          mean1 = mean1, var.mean1 = var.mean1,
                          sd1 = sd1, var.sd1 = var.sd1),
              ##
              n.iter.max = n.iter.max, tol = tol, iter = iter,
              call = match.call(),
              version = packageDescription("diagmeta")$Version
              )
  
  
  class(res) <- "diagmeta"
  
  res
}
