#' Plot for meta-analysis of diagnostic test accuracy studies with the
#' multiple cutoffs model
#' 
#' @description
#' Provides several plots for meta-analysis of diagnostic test
#' accuracy studies with the multiple cutoffs model
#' 
#' @param x An object of class \code{diagmeta}
#' @param which A character vector indicating the type of plot, either
#'   \code{"regression"} or \code{"cdf"} or \code{"survival"} or
#'   \code{"Youden"} or \code{"ROC"} or \code{"SROC"} or
#'   \code{"density"} or \code{"sensspec"}, can be abbreviated
#' @param xlab An x axis label
#' @param main A logical indicating title to the plot
#' @param ci A logical indicating whether confidence intervals should
#'   be plotted for \code{"regression"}, \code{"cdf"},
#'   \code{"survival"}, \code{"Youden"}, and \code{"sensspec"}
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
#'   in plots \code{"regression"}, \code{"cdf"}, \code{"survival"},
#'   \code{"Youden"}, \code{"ROC"}, and \code{"sensspec"}
#' @param lines A logical indicating whether polygonal lines
#'   connecting points belonging to the same study should be printed
#'   in plots \code{"regression"}, \code{"cdf"}, \code{"survival"},
#'   \code{"Youden"}, and \code{"sensspec"}
#' @param rlines A logical indicating whether regression lines or
#'   curves should be plotted for plots \code{"regression"},
#'   \code{"cdf"}, \code{"survival"}, \code{"Youden"}, and
#'   \code{"sensspec"}
#' @param line.optcut A logical indicating whether a vertical line
#'   should be plotted at the optimal cutoff line for plots
#'   \code{"cdf"}, \code{"survival"}, \code{"Youden"}, and
#'   \code{"density"}
#' @param col.points A character string indicating color of points,
#'   either \code{"rainbow"}, \code{"topo"}, \code{"heat"},
#'   \code{"terrain"}, \code{"cm"}, \code{"grayscale"}, or any color
#'   defined in \code{\link[grDevices]{colours}}
#' @param cex A numeric indicating magnification to be used for
#'   plotting text and symbols
#' @param pch.points A numeric indicating plot symbol(s) for points
#' @param col A character string indicating color of lines
#' @param col.ci A character string indicating color of confidence
#'   lines
#' @param col.optcut A character string indicating color of optimal
#'   cutoff line
#' @param cex.marks A numeric indicating magnification(s) to be used
#'   for marking cutoffs
#' @param lwd A numeric indicating line width
#' @param lwd.ci A numeric indicating line width of confidence lines
#' @param lwd.optcut A numeric indicating line width of optimal cutoff
#' @param lwd.study A numeric indicating line width of individual
#'   studies
#' @param shading A character indicating shading and hatching
#'   confidence region in \code{"SROC"} plot, either \code{"none"} or
#'   \code{"shade"} or \code{"hatch"}
#' @param col.hatching A character string indicating color used in
#'   hatching of the confidence region
#' @param lwd.hatching A numeric indicating line width used in hatching
#'   of the confidence region
#' @param ellipse A logical indicating whether a confidence ellipse
#'   should be drawn around the optimal cutoff
#' @param xlim A character or numerical vector indicating the minimum
#'   and maximum value for the horizontal axes
#' @param ylim A numerical vector indicating the minimum and maximum value for
#'   the vertical axes
#' @param \dots Additional graphical arguments
#' 
#' @details
#' The first argument of the plot function is an object of class
#' "diagmeta".
#' 
#' The second argument \code{which} indicates which sort of plot(s)
#' should be shown. For \code{which="regression"}, a scatter plot of
#' the quantile-transformed proportions of negative test results with
#' two regression lines is shown. Points belonging to the same study
#' are marked with the same colour. For \code{which="cdf"}, the two
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
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de},
#' Susanne Steinhauser \email{susanne.steinhauser@@uni-koeln.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de},
#' Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}}
#' 
#' @keywords hplot
#' 
#' @references
#' Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017):
#' A novel statistical model for analyzing data of a systematic review
#' generates optimal cutoff values for fractional exhaled nitric oxide
#' for asthma diagnosis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{92}, 69--78
#' 
#' Steinhauser S, Schumacher M, Rücker G (2016):
#' Modelling multiple thresholds in meta-analysis of diagnostic test
#' accuracy studies.
#' \emph{BMC Medical Research Methodology},
#' \bold{16}, 97
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
#' 
#' # Regression plot with confidence intervals
#' #
#' plot(diag1, which = "regr", lines = FALSE, ci = TRUE)
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
#' @method plot diagmeta
#' @export
#' @export plot.diagmeta
#'
#' @importFrom grDevices colours cm.colors heat.colors rainbow rgb
#'   terrain.colors topo.colors
#' @importFrom graphics abline curve par plot polygon text
#' @importFrom stats qchisq


plot.diagmeta <- function(x,
                          which = c("regression", "cdf", "sensspec",
                                    "youden", "roc", "sroc"),
                          xlab = "Threshold",
                          main,
                          ci = FALSE, ciSens = FALSE, ciSpec = FALSE,
                          mark.optcut = FALSE, mark.cutpoints = FALSE,
                          points = TRUE, lines = FALSE,
                          rlines = TRUE, line.optcut = TRUE,
                          col.points = "rainbow",
                          cex = 1, pch.points = 16,
                          col = "black",
                          col.ci = "gray",
                          col.optcut = "black",
                          cex.marks = 0.7 * cex,
                          lwd = 1, lwd.ci = lwd, lwd.optcut = 2 * lwd,
                          lwd.study = lwd,
                          shading = "none",
                          col.hatching = col.ci, lwd.hatching = lwd.ci,
                          ellipse = FALSE,
                          xlim = NULL, ylim = NULL,
                          ...) {
  
  
  chkclass(x, "diagmeta")
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
  mains <- c("Regression lines", "Cumulative distribution functions",
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
  chkchar(col.points)
  col.points <- setchar(col.points,
                        c("rainbow", "topo", "heat", "terrain",
                          "cm", "grayscale", colours()),
                        text = paste0("should be 'rainbow', 'topo', 'heat', ",
                                      "'terrain', 'cm', 'grayscale', ",
                                      "or any color defined in colours()"))
  ##
  col <- setchar(col, c("transparent", colours()),
                 text = paste0("should be 'transparent' or ",
                               "any color defined in colours()"))
  ##
  col.ci <- setchar(col.ci, c("transparent", colours()),
                    text = paste0("should be 'transparent' or ",
                                  "any color defined in colours()"))
  ##
  col.optcut <- setchar(col.optcut, c("transparent", colours()),
                        text = paste0("should be 'transparent' or ",
                                      "any color defined in colours()"))
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
  else if (col.points == "grayscale")
    cols <- gray(1:k / (k + 1))[study.no]
  else
    cols <- rep(col.points, k)[study.no]
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
  
  
  direction <- replaceNULL(x$direction, "increasing")
  #
  min.cutoff <- replaceNULL(x$min.cutoff, min(x$data.lmer$Cutoff, na.rm = TRUE))
  max.cutoff <- replaceNULL(x$max.cutoff, max(x$data.lmer$Cutoff, na.rm = TRUE))
  #
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
  ##
  log.cutoff <- x$log.cutoff
  ##
  log.axis <- if (log.cutoff) "x" else ""
  ##
  youden <- calcYouden.SeSp(Sens, Spec, lambda)
    
  
  for (i in which)
    do.call(i,
            list(cutoff = cutoff,
                 Sens = Sens, Spec = Spec,
                 studlab = studlab,
                 direction = direction, distr = distr, log.cutoff = log.cutoff,
                 xlab = xlab, ylab = ylab, xlim = xlim, log.axis = log.axis,
                 ylim = ylim,
                 mains = mains, which = which,
                 lines = lines, lwd.study = lwd.study, rlines = rlines,
                 alpha0 = alpha0, var.alpha0 = var.alpha0,
                 beta0 = beta0, var.beta0 = var.beta0,
                 alpha1 = alpha1, var.alpha1 = var.alpha1,
                 beta1 = beta1, var.beta1 = var.beta1,
                 cov.alpha0.alpha1 = cov.alpha0.alpha1,
                 cov.alpha0.beta0 = cov.alpha0.beta0,
                 cov.alpha0.beta1 = cov.alpha0.beta1,
                 cov.alpha1.beta0 = cov.alpha1.beta0,
                 cov.alpha1.beta1 = cov.alpha1.beta1,
                 cov.beta0.beta1 = cov.beta0.beta1,
                 var.nondiseased = var.nondiseased, var.diseased = var.diseased,
                 lambda = lambda,
                 cex = cex, col = col, lwd = lwd,
                 points = points, col.points = col.points,
                 pch.points = pch.points,
                 ci = ci, ciSens = ciSens, ciSpec = ciSpec,
                 level = level, col.ci = col.ci, lwd.ci = lwd.ci,
                 min.cutoff = min.cutoff, max.cutoff = max.cutoff,
                 mark.cutpoints = mark.cutpoints,
                 optcut = optcut,
                 mark.optcut = mark.optcut,
                 line.optcut = line.optcut, col.optcut = col.optcut,
                 lwd.optcut = lwd.optcut,
                 cex.marks = cex.marks,
                 ellipse = ellipse, shading = shading,
                 col.hatching = col.hatching, lwd.hatching = lwd.hatching,
                 youden = youden,
                 ...))
  
  
  invisible(NULL)
}
