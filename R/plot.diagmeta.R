#' Plot for meta-analysis of diagnostic test accuracy studies with the
#' multiple cutoffs model
#' 
#' Provides several plots for meta-analysis of diagnostic test
#' accuracy studies with the multiple cutoffs model
#' 
#' The first argument of the plot function is an object of class
#' "diagmeta".
#' 
#' The second argument \code{which} indicates which sort of plot(s)
#' should be shown. For \code{which="reg"}, a scatter plot of the
#' quantil-transformed proportions of negative test results with two
#' regression lines is shown.  Points belonging to the same study are
#' marked with the same colour. For \code{which="cdf"}, the two
#' cumulative distribution functions are shown, corresponding to the
#' proportions of negative test results. For \code{which="survival"},
#' the survival functions are shown, corresponding to the proportions
#' of positive test results. For \code{which="Youden"}, the
#' (potentially weighted) sum of sensitivity and specificity minus 1
#' is shown; in case of \code{lambda=0.5} (the default) this is the
#' Youden index. For \code{which="ROC"}, study-specific ROC curves are
#' shown. For \code{which="SROC"}, the model-based summary ROC curve
#' is shown. For \code{which="density"}, the model-based densities of
#' both groups are shown. For \code{which="sensspec"}, the sensitivity
#' (decreasing with increasing cutoff) and the specificity (increasing
#' with increasing cutoff) are shown. Instead of character strings, a
#' numeric value or vector can be used to specify plots with numbers
#' corresponding to the following order of plots: "regression", "cdf",
#' "survival", "youden", "roc", "sroc", "density", and "sensspec".
#' 
#' Other arguments refer to further plot parameters, such as
#' \code{lines} (whether points belonging to the same study should be
#' joined), \code{rlines} (whether regression curves should be drawn),
#' \code{ci} / \code{ciSens} / \code{ciSpec} / \code{ellipse} (whether
#' confidence regions should be shown), \code{line.optcut} /
#' \code{mark.optcut} (whether the optimal cutoff should be
#' indicated), and additional plot parameters (see Arguments).
#' 
#' If no further arguments are provided, four standard plots
#' ("survival", "Youden", "ROC", and "SROC") are produced in a 2 x 2
#' format.
#' 
#' @param x An object of class \code{diagmeta}
#' @param which A character vector indicating the type of plot, either
#'   \code{"reg"} or\code{"cdf"} or \code{"survival"} or
#'   \code{"Youden"} or \code{"ROC"} or \code{"SROC"} or
#'   \code{"density"}
#' @param xlab An x axis label
#' @param main A logical indicating title to the plot
#' @param ci A logical indicating whether confidence intervals should
#'   be plotted for \code{"reg"}, \code{"cdf"}, \code{"survival"},
#'   \code{"Youden"}, and \code{"sensspec"}
#' @param ciSens A logical indicating whether confidence intervals
#'   should be plotted for sensitivity, given the specificity in
#'   \code{"SROC"} plot
#' @param ciSpec A logical indicating whether confidence intervals
#'   should be plotted for specificity, given the sensitivity in
#'   \code{"SROC"} plot
#' @param mark.optcut A logical indicating whether the optimal cutoff
#'   should be marked on \code{"SROC"} plot
#' @param mark.cutpoints A logical indicating whether the given
#'   cutoffs should be marked on \code{"SROC"} plot
#' @param points A logical indicating whether points should be plotted
#'   on \code{"ROC"} plot
#' @param lines A logical indicating whether polygonal lines
#'   connecting points belonging to the same study should be printed
#'   in plots \code{"reg"}, \code{"cdf"}, \code{"survival"},
#'   \code{"Youden"}, and \code{"sensspec"}
#' @param rlines A logical indicating whether regression lines or
#'   curves should be plotted for plots \code{"reg"}, \code{"cdf"},
#'   \code{"survival"}, \code{"Youden"}, and \code{"sensspec"}
#' @param line.optcut A logical indicating whether a vertical line
#'   should be plotted at the optimal cutoff line for plots
#'   \code{"cdf"}, \code{"survival"}, \code{"Youden"}, and
#'   \code{"density"}
#' @param col.points A character vector indicating color, either
#'   \code{"rainbow"} or \code{"topo"} or \code{"heat"} or
#'   \code{"terrain"} or \code{"cm"} or \code{"gray"} or
#'   \code{"black"}
#' @param cex A numeric indicating magnification to be used for
#'   plotting text and symbols.
#' @param pch.points A numeric indicating plot symbol(s) for points
#' @param cex.marks A numeric indicating magnification(s) to be used
#'   for marking cutoffs
#' @param lwd A numeric indicating line width
#' @param lwd.optcut A numeric indicating line width of optimal cutoff
#' @param shading A character indicating shading and hatching
#'   confidence region in \code{"SROC"} plot, either \code{"none"} or
#'   \code{"shade"} or \code{"hatch"}
#' @param ellipse A logical indicating whether a confidence ellipse
#'   should be drawn around the optimal cutoff
#' @param xlim A character or numerical vector indicating the minimum
#'   and maximum value for the horizontal axes.
#' @param \dots Additional graphical arguments
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Susanne
#'   Steinhauser \email{susanne.steinhauser@@uni-koeln.de}, Srinath
#'   Kolampally \email{kolampal@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}}
#' 
#' @references
#'
#' Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017), A
#' novel statistical model for analyzing data of a systematic review
#' generates optimal cutoff values for fractional exhaled nitric oxide
#' for asthma diagnosis \emph{Journal of Clinical Epidemiology},
#' \bold{92}, 69--78, doi: 10.1016/j.jclinepi.2017.09.001 .
#' 
#' Steinhauser S, Schumacher M, Rücker G (2016), Modelling multiple
#' thresholds in meta-analysis of diagnostic test accuracy
#' studies. \emph{BMC Med Res Methodol}, \bold{16}, 97.
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
#' 
#' # Regression plot with confidence intervals
#' #
#' plot(diag1, which = "reg", lines = FALSE, ci = TRUE)
#' 
#' # Cumulative distribution plot with optimal cutoff line and
#' # confidence intervals
#' #
#' plot(diag1, which = "cdf", line.optcut = TRUE, ci = TRUE)
#' 
#' # Survival plot with optimal cutoff line and confidence intervals
#' #
#' plot(diag1, which = "survival", line.optcut = TRUE, ci = TRUE)
#' 
#' # Youden plot of optimal cutoff line and confidence intervals
#' #
#' plot(diag1, which = "youden",
#'      lines = TRUE, line.optcut = TRUE, ci = TRUE)
#' 
#' # ROC plot of lines connecting points belonging to the same study
#' #
#' plot(diag1, which = "ROC", lines = TRUE)
#' 
#' # SROC plot of confidence regions for sensitivity and specificity
#' # with optimal cutoff mark
#' #
#' plot(diag1, which = "SROC",
#'      ciSens = TRUE, ciSpec = TRUE, mark.optcut = TRUE,
#'      shading = "hatch")
#' 
#' # Density plot of densities for both groups with optimal cutoff
#' # line
#' #
#' plot(diag1, which = "density", line.optcut = TRUE)
#' 
#' @export plot.diagmeta
#' @export
#'
#' @importFrom grDevices cm.colors heat.colors rainbow rgb
#'   terrain.colors topo.colors
#' @importFrom graphics abline curve par plot polygon text
#' @importFrom stats qchisq


plot.diagmeta <- function(x,
                          which = c("survival", "youden", "roc", "sroc"),
                          xlab = "threshold",
                          main,
                          ci = FALSE, ciSens = FALSE, ciSpec = FALSE,
                          mark.optcut = FALSE, mark.cutpoints = FALSE,
                          points = TRUE, lines = FALSE,
                          rlines = TRUE, line.optcut = TRUE,
                          col.points = "rainbow",
                          cex = 1, pch.points = 16,
                          cex.marks = 0.7 * cex,
                          lwd = 1, lwd.optcut = 2 * lwd,
                          shading = "none",
                          ellipse = FALSE,
                          xlim = NULL,
                          ...) {
  
  
  meta:::chkclass(x, "diagmeta")
  ##
  setchar <- meta:::setchar
  chklength <- meta:::chklength
  chklogical <- meta:::chklogical
  ##
  plot.types <- c("regression", "cdf", "survival", "youden",
                  "roc", "sroc", "density", "sensspec")
  
  
  if (is.character(which))
    which <- setchar(which, plot.types)
  else if (is.numeric(which)) {
    if (any(which < 1) | any(which > 8))
      stop("Numeric values of argument 'which' must be in 1:8")
    which <- plot.types[which]
  }
  else
    stop("Argument 'which' must be a character vector or ",
         "a numeric vector with values between 1 and 8.")
  ##
  which <- unique(which)
  ##
  mains <- c("Regression lines", "Cumulative density functions",
             "Survival curves", "Youden index", "ROC curves",
             "SROC curve", "Density functions",
             "Sensitivity and Specificity")[match(which, plot.types)]
  ##
  if (!missing(main))
    if (is.logical(main)) {
      if (!main)
        mains <- rep_len("", length(mains))
    }
    else
      mains <- main
  ##
  chklength(mains, length(which), "diagmeta",
            text = paste("Number of plots (argument 'which') and number of",
                         "plot titles (argument 'main') differ."))
  ##
  chklogical(ci)
  chklogical(ciSens)
  chklogical(ciSpec)
  chklogical(mark.optcut)
  chklogical(mark.cutpoints)
  chklogical(points)
  chklogical(lines)
  chklogical(rlines)
  chklogical(line.optcut)
  ##
  col.points <- setchar(col.points,
                        c("rainbow", "topo", "heat", "terrain",
                          "cm", "gray", "black"))
  ##
  shading <- setchar(shading,
                     c("none", "hatch", "shade"))
  ##
  is.logistic <- x$distr == "logistic"
  is.normal <- x$distr == "normal"
  
  
  ##
  ## Define plot layout (depending on number of plots)
  ##
  n.plots <- length(which)
  ##
  if (n.plots == 1)
    oldpar <- par(pty = "s")
  else if (n.plots == 2)
    oldpar <- par(mfrow = c(1, 2), pty = "s")
  else if (n.plots <  5)
    oldpar <- par(mfrow = c(2, 2), pty = "s")
  else if (n.plots <  7)
    oldpar <- par(mfrow = c(2, 3), pty = "s")
  else
    oldpar <- par(mfrow = c(3, 3), pty = "s")
  ##
  on.exit(par(oldpar))
  

  k <- x$k
  ##
  o <- order(x$studlab, x$cutoff)
  ##
  cutoff <- x$cutoff[o]
  studlab <- x$studlab[o]
  Spec <- x$Spec[o]
  Sens <- x$Sens[o]
  
  
  study.no <- as.numeric(as.factor(studlab))
  ##  
  if (col.points == "rainbow")
    cols <- rainbow(k)[study.no]
  else if (col.points == "topo")
    cols <- topo.colors(k)[study.no]
  else if (col.points == "heat")
    cols <- heat.colors(k)[study.no]
  else if (col.points == "terrain")
    cols <- terrain.colors(k)[study.no]
  else if (col.points == "cm")
    cols <- cm.colors(k)[study.no]
  else if (col.points == "gray")
    cols <- gray(1:k / (k + 1))[study.no]
  else if (col.points == "black")
    cols <- rep(1, k)[study.no]
  ##
  col.points <- cols
  ##
  gray <- gray(0.75)
  
  
  ## Plot y axis label
  ##
  if (is.logistic)
    ylab <- expression(atop("", "Logit(P(negative test result))"))
  else if (is.normal)
    ylab <- expression(atop("", Phi^{-1}~"(P(negative test result))"))
  else
    ylab <- ""
  
  
  distr <- x$distr
  level <- x$level
  lambda <- x$lambda
  optcut <- x$optcut
  ##
  alpha0 <- x$regr$alpha0
  var.alpha0 <- x$regr$var.alpha0
  beta0 <- x$regr$beta0
  var.beta0 <- x$regr$var.beta0
  cov.alpha0.beta0 <- x$regr$cov.alpha0.beta0
  alpha1 <- x$regr$alpha1
  var.alpha1 <- x$regr$var.alpha1
  beta1 <- x$regr$beta1
  var.beta1 <- x$regr$var.beta1
  cov.alpha1.beta1 <- x$regr$cov.alpha1.beta1
  cov.alpha0.alpha1 <- x$regr$cov.alpha0.alpha1
  cov.alpha0.beta1 <- x$regr$cov.alpha0.beta1
  cov.alpha1.beta0 <- x$regr$cov.alpha1.beta0
  cov.beta0.beta1 <- x$regr$cov.beta0.beta1
  ##
  var.nondiseased <- x$var.nondiseased
  var.diseased <- x$var.diseased
  ##
  NN <- x$data.lmer$NN
  Cutoff <- x$data.lmer$Cutoff
  ##
  log.cutoff <- x$log.cutoff
  ##
  log.axis <- if (log.cutoff) "x" else ""
  ##
  youden <- calcYouden.SeSp(Sens, Spec, lambda)
  
  
  ##
  ##
  ## Linear regression lines in logit / probit space
  ##
  ##
  if ("regression" %in% which) {
    ##
    plot(c(cutoff, cutoff), qdiag(c(Spec, 1 - Sens), distr),
         type = "n", las = 1, log = log.axis,
         ylab = ylab, xlab = xlab,
         main = mains[match("regression", which)],
         xlim = xlim, ...)
    ##
    ## Add data
    ##
    if (lines)
      for (s in studlab) {
        lines(cutoff[studlab == s],
              qdiag(Spec[studlab == s], distr),
              col = col.points[studlab == s], lwd = lwd, lty = 2)
        ##
        lines(cutoff[studlab == s],
              qdiag(1 - Sens[studlab == s], distr),
              col = col.points[studlab == s], lwd = lwd, lty = 1)
      }
    ##
    points(cutoff, qdiag(Spec, distr),
           pch = 1, cex = cex, col = col.points)
    ##
    points(cutoff, qdiag(1 - Sens, distr),
           pch = pch.points, cex = cex, col = col.points)      
    ##
    ## Add linear regression lines
    ##
    if (log.cutoff) {
      if (rlines) {
        curve(alpha0 + beta0 * log(x),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        curve(alpha1 + beta1 * log(x),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        ##
        curve(ciRegr(log(x),
                     alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                     var.nondiseased,
                     level)$lower,
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(log(x),
                     alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                     var.nondiseased,
                     level)$upper,
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(log(x),
                     alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                     var.diseased,
                     level)$lower,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(log(x),
                     alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                     var.diseased,
                     level)$upper,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      if (rlines) {
        curve(alpha0 + beta0 * x,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        curve(alpha1 + beta1 * x,
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        ##
        curve(ciRegr(x,
                     alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                     var.nondiseased,
                     level)$lower,
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(x,
                     alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                     var.nondiseased,
                     level)$upper,
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(x,
                     alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                     var.diseased,
                     level)$lower,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciRegr(x,
                     alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                     var.diseased,
                     level)$upper,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
  }
  
  
  ##
  ##
  ## Data and biomarker distributions functions (empirical
  ## distribution function)
  ##
  ##
  if ("cdf" %in% which) {
    ##
    plot(c(cutoff, cutoff), c(Spec, 1 - Sens),
         type = "n", las = 1, log = log.axis,
         xlab = xlab, ylab = "Prob(negative test)",
         main = mains[match("cdf", which)],
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##
    if (lines)
      for (s in studlab) {
        lines(cutoff[studlab == s], 
              Spec[studlab == s],
              col = col.points[studlab == s], lwd = lwd, lty = 2)
        ##
        lines(cutoff[studlab == s], 
              1 - Sens[studlab == s],
              col = col.points[studlab == s], lwd = lwd, lty = 1)
      }
    ##
    ## Add data
    ##
    points(cutoff, 1 - Sens,
           pch = pch.points, cex = cex, col = col.points)
    ##
    points(cutoff, Spec,
           pch = 1, cex = cex, col = col.points)
    ##
    ## Add regression curves
    ##
    if (log.cutoff) {
      ##
      if (rlines) {
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      ##
      if (rlines) {
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
  }
  
  
  ##
  ##
  ## Data and biomarker distributions functions (survival function)
  ##
  ##
  if ("survival" %in% which) {
    ##
    plot(c(cutoff, cutoff), 1 - c(Spec, 1 - Sens),
         type = "n", las = 1, log = log.axis,
         xlab = xlab, ylab = "Prob(positive test)",
         main = mains[match("survival", which)],
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##
    if (lines) {
      for (s in studlab)
        lines(cutoff[studlab == s], 
              Sens[studlab == s],
              col = col.points[studlab == s], lwd = lwd, lty = 2)
      ##
      for (s in studlab)
        lines(cutoff[studlab == s], 
              1 - Spec[studlab == s],
              col = col.points[studlab == s], lwd = lwd, lty = 1)
    }
    ##
    ## Add data
    ##
    points(cutoff, Sens,
           pch = pch.points, cex = cex, col = col.points)
    ##
    points(cutoff, 1 - Spec, pch = 1, cex = cex, col = col.points)
    ##
    ## Add regression curves
    ##
    if (log.cutoff) {
      ##
      if (rlines) {
        curve(calcSens(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSens(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      ##
      if (rlines) {
        curve(calcSens(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSens(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
  }
  
  
  ##
  ##
  ## Plot of Youden index = Sens + Spec - 1 ~ TNR - FNR
  ## not for weighted Youden index!
  ## (-> need to weight the data, too)
  ##
  ##
  if ("youden" %in% which) {
    ##
    plot(cutoff, youden,
         type = "n",
         las = 1, log = log.axis,
         ylab = "(Weighted) Youden index", xlab = xlab,
         main = mains[match("youden", which)],
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##
    if (lines) {
      for (s in studlab)
        lines(cutoff[studlab == s],
              youden[studlab == s],
              col = col.points[studlab == s], lwd = lwd)
    }
    ##
    ## Add data
    ##
    points(cutoff, youden,
           pch = pch.points, col = col.points, cex = cex)
    ##
    ## Plot Youden index
    ##
    if (log.cutoff) {
      if (ci) {
        curve(ciYouden(log(x), distr, lambda,
                       alpha0, var.alpha0, beta0, var.beta0,
                       cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                       alpha1, var.alpha1, beta1, var.beta1,
                       cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                       var.nondiseased, var.diseased,
                       level)$lower,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciYouden(log(x), distr, lambda,
                       alpha0, var.alpha0, beta0, var.beta0,
                       cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                       alpha1, var.alpha1, beta1, var.beta1,
                       cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                       var.nondiseased, var.diseased,
                       level)$upper,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
      ##
      if (rlines)
        curve(calcYouden(log(x), distr, lambda,
                         alpha0, beta0, alpha1, beta1),
              col = 1, add = TRUE)     
    }
    else {
      if (ci) {
        curve(ciYouden(x, distr, lambda,
                       alpha0, var.alpha0, beta0, var.beta0,
                       cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                       alpha1, var.alpha1, beta1, var.beta1,
                       cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                       var.nondiseased, var.diseased,
                       level)$lower,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(ciYouden(x, distr, lambda,
                       alpha0, var.alpha0, beta0, var.beta0,
                       cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                       alpha1, var.alpha1, beta1, var.beta1,
                       cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                       var.nondiseased, var.diseased,
                       level)$upper,
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
      ##
      if (rlines)
        curve(calcYouden(x, distr, lambda,
                         alpha0, beta0, alpha1, beta1),
              col = 1, add = TRUE)     
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
  }
  
  
  ##
  ##
  ## Study-specific ROC curves
  ##
  ##
  if ("roc" %in% which) {
    ##
    plot(1 - Spec, Sens,
         type = "n", las = 1,
         xlab = "1 - Specificity", ylab = "Sensitivity",
         main = mains[match("roc", which)],
         xlim = c(0, 1), ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##
    for (s in studlab)
      lines(c(1, 1 - Spec[studlab == s], 0),
            c(1, Sens[studlab == s], 0),
            col = col.points[studlab == s], lwd = lwd)
    ##
    ## Add data
    ##
    if (points)
      points(1 - Spec, Sens,
             pch = pch.points, col = col.points, cex = cex)
  }
  
  
  ##
  ##
  ## SROC curve
  ##
  ##
  if ("sroc" %in% which) {
    ##
    plot(1 - Spec, Sens,
         type = "n", las = 1,
         xlab = "1 - Specificity", ylab = "Sensitivity",
         main = mains[match("sroc", which)],
         xlim = c(0, 1), ylim = c(0, 1), ...)
    ##
    if (log.cutoff)
      cuts <- log(unique(cutoff))
    else
      cuts <- unique(cutoff)
    ##
    cuts <- cuts[order(cuts)]
    ##
    tcuts0 <- ciRegr(cuts,
                     alpha0, var.alpha0, beta0, var.beta0,
                     cov.alpha0.beta0, var.nondiseased,
                     level)
    ##
    tcuts1 <- ciRegr(cuts,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha1.beta1, var.diseased,
                     level)
    ##
    Sens.cuts <- calcSens(tcuts1$TE, distr)
    lowerSens.cuts <- calcSens(tcuts1$lower, distr)
    upperSens.cuts <- calcSens(tcuts1$upper, distr)
    ##
    Spec.cuts <- calcSpec(tcuts0$TE, distr)
    lowerSpec.cuts <- calcSpec(tcuts0$lower, distr)
    upperSpec.cuts <- calcSpec(tcuts0$upper, distr)
    ##
    y.lower.Se <- lowerSens.cuts
    y.lower.Sp <- Sens.cuts
    ##
    y.upper.Se <- upperSens.cuts
    y.upper.Sp <- Sens.cuts
    ##
    x.lower.Se <- 1 - Spec.cuts
    x.lower.Sp <- 1 - lowerSpec.cuts
    ##
    x.upper.Se <- 1 - Spec.cuts
    x.upper.Sp <- 1 - upperSpec.cuts
    ##
    if (log.cutoff)
      ocut <- log(optcut)
    else
      ocut <- optcut
    ##
    ce <- ciEllipse(ocut, 
                    alpha0, var.alpha0, beta0, var.beta0,
                    alpha1, var.alpha1, beta1, var.beta1,
                    cov.alpha0.beta0, cov.alpha1.beta1,
                    cov.alpha0.alpha1, cov.alpha0.beta1,
                    cov.alpha1.beta0, cov.beta0.beta1, 
                    var.nondiseased, var.diseased, level)
    
    ##
    if (ciSens) {
      if(shading == "shade")
        polygon(c(x.upper.Se,x.lower.Se[order(x.lower.Se)]),
                c(y.upper.Se,y.lower.Se[order(x.lower.Se)]),
                col = rgb(0.5, 0.5, 0.5, alpha = 0.2), border = NA)
      else if (shading == "hatch")
        polygon(c(x.upper.Se, x.lower.Se[order(x.lower.Se)]),
                c(y.upper.Se, y.lower.Se[order(x.lower.Se)]),
                density = 20, angle = 90, col = "gray", border = "gray")
      ##
      lines(x.upper.Se, y.upper.Se, col = "darkgray")
      lines(x.lower.Se, y.lower.Se, col = "darkgray")
    }
    ##
    if (ciSpec) {
      if (shading == "shade")
        polygon(c(x.upper.Sp, x.lower.Sp[order(x.lower.Sp)]),
                c(y.upper.Sp, y.lower.Sp[order(x.lower.Sp)]),
                col = rgb(0.2, 0.2, 0.2, alpha = 0.2), border = NA)
      else if (shading == "hatch")
        polygon(c(x.upper.Sp, x.lower.Sp[order(x.lower.Sp)]),
                c(y.upper.Sp, y.lower.Sp[order(x.lower.Sp)]),
                density = 20, angle = 0, col = "gray", border = "gray")
      ##
      lines(x.upper.Sp, y.upper.Sp, col = "darkgray", lty = 2)
      lines(x.lower.Sp, y.lower.Sp, col = "darkgray", lty = 2)
    }
    ##
    ## Add ROC curve
    ##
    curve(pdiag(alpha1 + beta1 *
                (qdiag(1 - x, distr) - alpha0) / beta0,
                distr, FALSE),
          lwd = lwd, col = 1, add = TRUE)
    ##
    if (mark.optcut) {
      if (log.cutoff)
        ocut <- log(optcut)
      else
        ocut <- optcut
      ##
      points(pdiag(alpha0 + beta0 * ocut, distr, FALSE),
             pdiag(alpha1 + beta1 * ocut, distr, FALSE),
             lwd = lwd.optcut, cex = 2, pch = 3, col = 1)
    }
    ##
    ## Add data
    ##
    points(1 - Spec, Sens,
           pch = pch.points, col = col.points, cex = cex)
    ##
    ## Add text
    ##
    if (mark.cutpoints) {
      if (log.cutoff)
        cuts <- log(unique(cutoff))
      else
        cuts <- unique(cutoff)
      ##
      text.cuts <- as.character(round(unique(cutoff), 2))
      ##
      points(pdiag(alpha0 + beta0 * cuts, distr, FALSE),
             pdiag(alpha1 + beta1 * cuts, distr, FALSE),
             pch = 3, cex = cex)
      ##
      text(1.02 - pdiag(alpha0 + beta0 * cuts, distr),
           0.98 - pdiag(alpha1 + beta1 * cuts, distr),
           text.cuts, cex = cex.marks)
    }
    ##
    ## Add ellipse
    ##
    if (ellipse) {
      n <- 2 * pi * 100 
      xx <- yy <- rep(0, n) 
      q <- sqrt(qchisq(0.05, 2, lower.tail = FALSE))
      ##
      for (t in 1:n) {
        xx[t] <- pdiag(-ce$logit.spec +
                       q * ce$se.y0 * cos(t / 100 + acos(ce$r)), distr)
        yy[t] <- pdiag(ce$logit.sens +
                       q * ce$se.y1 * cos(t / 100), distr)
      }
      lines(xx, yy)
    }
  }
  
  
  ##
  ##
  ## Biomarker distributions (densities)
  ##
  ##
  if ("density" %in% which) {
    ##
    if (log.cutoff) {
      curve(beta0 * ddiag(beta0 * log(x) + alpha0, distr),
            lty = 2, log = "x", las = 1,
            from = min(exp(Cutoff)), to = max(exp(Cutoff)),
            xlab = xlab, ylab = "",
            main = mains[match("density", which)])
      ##
      curve(beta1 * ddiag(beta1 * log(x) + alpha1, distr),
            lty = 1, add = TRUE)
    }
    else {
      curve(beta0 * ddiag(beta0 * x + alpha0, distr),
            lty = 2, las = 1,
            from = min(Cutoff), to = max(Cutoff),
            xlab = xlab, ylab = "",
            main = mains[match("density", which)])
      ##
      curve(beta1 * ddiag(beta1 * x + alpha1, distr),
            lty = 1, add = TRUE)
    }
    ##
    ## draw optimal cut-off
    ##
    if (line.optcut)
      abline(v = optcut)
  }
  
  
  ##
  ##
  ## Sensitivity and specificity
  ##
  ##
  if ("sensspec" %in% which) {
    ##
    plot(c(cutoff, cutoff), c(Spec, Sens),
         type = "n",
         las = 1, log = log.axis,
         ylab = "Sensitivity / Specificity", xlab = xlab,
         main = mains[match("sensspec", which)],
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##    
    if (lines) {
      for (s in studlab) {
        lines(cutoff[studlab == s], Spec[studlab == s],
              col = col.points[studlab == s], lwd = lwd,
              lty = 2)
        lines(cutoff[studlab == s], Sens[studlab == s],
              col = col.points[studlab == s], lwd = lwd,
              lty = 1)
      }
    }
    ##
    ## Add data
    ##
    points(cutoff, Sens, pch = pch.points, cex = cex,
           col = col.points)
    points(cutoff, Spec, pch = 1, cex = cex, col = col.points)
    ##
    ## Add curves
    ##
    if (log.cutoff) {
      ##
      if (rlines) {
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(log(x),
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(log(x),
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      ##
      if (rlines) {
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$TE,
                       distr),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$TE,
                       distr),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$lower,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ciRegr(x,
                              alpha0, var.alpha0,
                              beta0, var.beta0,
                              cov.alpha0.beta0,
                              var.nondiseased,
                              level)$upper,
                       distr),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$lower,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ciRegr(x,
                              alpha1, var.alpha1,
                              beta1, var.beta1,
                              cov.alpha1.beta1,
                              var.diseased,
                              level)$upper,
                       distr),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut) 
      abline(v = optcut)
  }  
  
  invisible(NULL)
}
