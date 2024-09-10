#
# Linear regression lines in logit / probit space
#
regression <- function(cutoff, Sens, Spec, studlab,
                       direction, distr, log.cutoff,
                       xlab, ylab, xlim, log.axis, ylim,
                       mains, which,
                       lines, lwd.study, rlines,
                       alpha0, var.alpha0, beta0, var.beta0,
                       alpha1, var.alpha1, beta1, var.beta1,
                       cov.alpha0.alpha1, cov.alpha0.beta0,
                       cov.alpha0.beta1, cov.alpha1.beta0,
                       cov.alpha1.beta1, cov.beta0.beta1,
                       var.nondiseased, var.diseased,
                       lambda,
                       col, lwd,
                       points, cex, col.points, pch.points,
                       ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                       min.cutoff, max.cutoff, mark.cutpoints,
                       optcut, mark.optcut,
                       line.optcut, col.optcut,
                       lwd.optcut,
                       cex.marks,
                       ellipse, shading,
                       col.hatching, lwd.hatching,
                       youden,
                       x,
                       ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  plot(c(cutoff, cutoff), qdiag(c(Spec, 1 - Sens), distr),
       type = "n", las = 1, log = log.axis,
       ylab = ylab, xlab = xlab,
       main = mains[match("regression", which)],
       xlim = xlim, ...)
  #
  # Add data
  #
  if (lines)
    for (s in studlab) {
      lines(cutoff[studlab == s], qdiag(Spec[studlab == s], distr),
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
      #
      lines(cutoff[studlab == s], qdiag(1 - Sens[studlab == s], distr),
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    }
  #
  if (points) {
    points(cutoff, qdiag(Spec, distr), pch = 1, cex = cex, col = col.points)
    ##
    points(cutoff, qdiag(1 - Sens, distr), pch = pch.points, cex = cex,
           col = col.points)
  }
  #
  # Add linear regression lines
  #
  if (rlines) {
    lines(xvals, alpha0 + beta0 * xvals.tr, lty = 2, col = col, lwd = lwd)
    lines(xvals, alpha1 + beta1 * xvals.tr, lty = 1, col = col, lwd = lwd)
  }
  #
  if (ci) {
    #
    lines(xvals,
          ciRegr(xvals.tr,
                 alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                 var.nondiseased,
                 level)$lower,
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          ciRegr(xvals.tr,
                 alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                 var.nondiseased,
                 level)$upper,
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          ciRegr(xvals.tr,
                 alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                 var.diseased,
                 level)$lower,
          lty = 1, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          ciRegr(xvals.tr,
                 alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                 var.diseased,
                 level)$upper,
          lty = 1, col = col.ci, lwd = lwd.ci)
  }
  #
  invisible(NULL)
}


##
## Data and biomarker distributions functions (empirical
## distribution function)
##
cdf <- function(cutoff, Sens, Spec, studlab,
                direction, distr, log.cutoff,
                xlab, ylab, xlim, log.axis, ylim,
                mains, which,
                lines, lwd.study, rlines,
                alpha0, var.alpha0, beta0, var.beta0,
                alpha1, var.alpha1, beta1, var.beta1,
                cov.alpha0.alpha1, cov.alpha0.beta0,
                cov.alpha0.beta1, cov.alpha1.beta0,
                cov.alpha1.beta1, cov.beta0.beta1,
                var.nondiseased, var.diseased,
                lambda,
                col, lwd,
                points, cex, col.points, pch.points,
                ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                min.cutoff, max.cutoff, mark.cutpoints,
                optcut, mark.optcut,
                line.optcut, col.optcut,
                lwd.optcut,
                cex.marks,
                ellipse, shading,
                col.hatching, lwd.hatching,
                youden,
                x,
                ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  if (direction == "decreasing")
    xlim <- c(xlim[2], xlim[1])
  #
  plot(c(cutoff, cutoff), c(Spec, 1 - Sens),
       type = "n", las = 1, log = log.axis,
       xlab = xlab, ylab = "Prob(negative test)",
       main = mains[match("cdf", which)],
       xlim = xlim, ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  #
  # Add lines
  #
  if (lines)
    for (s in studlab) {
      lines(cutoff[studlab == s], Spec[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
      #
      lines(cutoff[studlab == s], 1 - Sens[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    }
  #
  # Add data
  #
  if (points) {
    points(cutoff, 1 - Sens, pch = pch.points, cex = cex, col = col.points)
    points(cutoff, Spec, pch = 1, cex = cex, col = col.points)
  }
  #
  # Add regression curves
  #
  if (rlines) {
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$TE,
                   distr),
          lty = 2, col = col, lwd = lwd)
    ##
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$TE,
                   distr),
          lty = 1, col = col, lwd = lwd)
  }
  #
  if (ci) {
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$lower,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$upper,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$lower,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$upper,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
  }
  ##
  ## Draw line for optimal cutoff
  ##
  if (line.optcut)
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}


##
## Data and biomarker distributions functions (survival function)
##
survival <- function(cutoff, Sens, Spec, studlab,
                     direction, distr, log.cutoff,
                     xlab, ylab, xlim, log.axis, ylim,
                     mains, which,
                     lines, lwd.study, rlines,
                     alpha0, var.alpha0, beta0, var.beta0,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha0.alpha1, cov.alpha0.beta0,
                     cov.alpha0.beta1, cov.alpha1.beta0,
                     cov.alpha1.beta1, cov.beta0.beta1,
                     var.nondiseased, var.diseased,
                     lambda,
                     col, lwd,
                     points, cex, col.points, pch.points,
                     ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                     min.cutoff, max.cutoff, mark.cutpoints,
                     optcut, mark.optcut,
                     line.optcut, col.optcut,
                     lwd.optcut,
                     cex.marks,
                     ellipse, shading,
                     col.hatching, lwd.hatching,
                     youden,
                     x,
                     ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  if (direction == "decreasing")
    xlim <- c(xlim[2], xlim[1])
  #
  plot(c(cutoff, cutoff), 1 - c(Spec, 1 - Sens),
       type = "n", las = 1, log = log.axis,
       xlab = xlab, ylab = "Prob(positive test)",
       main = mains[match("survival", which)],
       xlim = xlim, ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  #
  # Add lines
  #
  if (lines) {
    for (s in studlab)
      lines(cutoff[studlab == s], Sens[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    #
    for (s in studlab)
      lines(cutoff[studlab == s], 1 - Spec[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
  }
  #
  # Add data
  #
  if (points) {
    points(cutoff, Sens, pch = pch.points, cex = cex, col = col.points)
    ##
    if (points)
      points(cutoff, 1 - Spec, pch = 1, cex = cex, col = col.points)
  }
  #
  # Add regression curves
  #
  if (rlines) {
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$TE,
                   distr),
          lty = 2, col = col, lwd = lwd)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$TE,
                   distr),
          lty = 1, col = col, lwd = lwd)
  }
  #
  if (ci) {
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$lower,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$upper,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$lower,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$upper,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
  }
  ##
  ## Draw line for optimal cutoff
  ##
  if (line.optcut)
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}


##
## Plot of Youden index = Sens + Spec - 1 ~ TNR - FNR
## not for weighted Youden index!
## (-> need to weight the data, too)
##
youden <- function(cutoff, Sens, Spec, studlab,
                   direction, distr, log.cutoff,
                   xlab, ylab, xlim, log.axis, ylim,
                   mains, which,
                   lines, lwd.study, rlines,
                   alpha0, var.alpha0, beta0, var.beta0,
                   alpha1, var.alpha1, beta1, var.beta1,
                   cov.alpha0.alpha1, cov.alpha0.beta0,
                   cov.alpha0.beta1, cov.alpha1.beta0,
                   cov.alpha1.beta1, cov.beta0.beta1,
                   var.nondiseased, var.diseased,
                   lambda,
                   col, lwd,
                   points, cex, col.points, pch.points,
                   ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                   min.cutoff, max.cutoff, mark.cutpoints,
                   optcut, mark.optcut,
                   line.optcut, col.optcut,
                   lwd.optcut,
                   cex.marks,
                   ellipse, shading,
                   col.hatching, lwd.hatching,
                   youden,
                   x,
                   ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  plot(cutoff, youden,
       type = "n",
       las = 1, log = log.axis,
       ylab = "(Weighted) Youden index", xlab = xlab,
       main = mains[match("youden", which)],
       xlim = xlim, ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  #
  # Add lines
  #
  if (lines) {
    for (s in studlab)
      lines(cutoff[studlab == s], youden[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study)
  }
  #
  # Add data
  #
  if (points)
    points(cutoff, youden, pch = pch.points, col = col.points, cex = cex)
  #
  # Plot Youden index
  #
  if (ci) {
    lines(xvals,
          ciYouden(xvals.tr,
                   distr, lambda,
                   alpha0, var.alpha0, beta0, var.beta0,
                   cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                   alpha1, var.alpha1, beta1, var.beta1,
                   cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                   var.nondiseased, var.diseased,
                   level)$lower,
          lty = 1, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          ciYouden(xvals.tr,
                   distr, lambda,
                   alpha0, var.alpha0, beta0, var.beta0,
                   cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                   alpha1, var.alpha1, beta1, var.beta1,
                   cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                   var.nondiseased, var.diseased,
                   level)$upper,
          lty = 1, col = col.ci, lwd = lwd.ci)
  }
  #
  if (rlines)
    lines(xvals,
          calcYouden(xvals.tr, distr, lambda, alpha0, beta0, alpha1, beta1),
          col = col, lwd = lwd)
  #
  # Draw line for optimal cutoff
  #
  if (line.optcut)
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}


##
## Study-specific ROC curves
##
roc <- function(cutoff, Sens, Spec, studlab,
                direction, distr, log.cutoff,
                xlab, ylab, xlim, log.axis, ylim,
                mains, which,
                lines, lwd.study, rlines,
                alpha0, var.alpha0, beta0, var.beta0,
                alpha1, var.alpha1, beta1, var.beta1,
                cov.alpha0.alpha1, cov.alpha0.beta0,
                cov.alpha0.beta1, cov.alpha1.beta0,
                cov.alpha1.beta1, cov.beta0.beta1,
                var.nondiseased, var.diseased,
                lambda,
                col, lwd,
                points, cex, col.points, pch.points,
                ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                min.cutoff, max.cutoff, mark.cutpoints,
                optcut, mark.optcut,
                line.optcut, col.optcut,
                lwd.optcut,
                cex.marks,
                ellipse, shading,
                col.hatching, lwd.hatching,
                youden,
                x,
                ...) {
  #
  if (direction == "increasing") {
    xy1 <- 1
    xy2 <- 0
  }
  else {
    xy1 <- 0
    xy2 <- 1
  }
  #
  plot(1 - Spec, Sens,
       type = "n", las = 1,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = mains[match("roc", which)],
       xlim = c(0, 1), ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  ##
  ## Add lines
  ##
  for (s in studlab)
    lines(c(xy1, 1 - Spec[studlab == s], xy2),
          c(xy1, Sens[studlab == s], xy2),
          col = col.points[studlab == s], lwd = lwd.study)
  ##
  ## Add data
  ##
  if (points)
    points(1 - Spec, Sens, pch = pch.points, col = col.points, cex = cex)
  ##
  invisible(NULL)
}


##
## SROC curve
##
sroc <- function(cutoff, Sens, Spec, studlab,
                 direction, distr, log.cutoff,
                 xlab, ylab, xlim, log.axis, ylim,
                 mains, which,
                 lines, lwd.study, rlines,
                 alpha0, var.alpha0, beta0, var.beta0,
                 alpha1, var.alpha1, beta1, var.beta1,
                 cov.alpha0.alpha1, cov.alpha0.beta0,
                 cov.alpha0.beta1, cov.alpha1.beta0,
                 cov.alpha1.beta1, cov.beta0.beta1,
                 var.nondiseased, var.diseased,
                 lambda,
                 col, lwd,
                 points, cex, col.points, pch.points,
                 ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                 min.cutoff, max.cutoff, mark.cutpoints,
                 optcut, mark.optcut,
                 line.optcut, col.optcut,
                 lwd.optcut,
                 cex.marks,
                 ellipse, shading,
                 col.hatching, lwd.hatching,
                 youden,
                 x,
                 ...) {
  #
  plot(1 - Spec, Sens,
       type = "n", las = 1,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = mains[match("sroc", which)],
       xlim = c(0, 1), ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  #
  cuts <- unique(cutoff)
  cuts <- cuts[order(cuts)]
  #
  cuts.tr <- transf(cuts, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  tcuts0 <- ciRegr(cuts.tr,
                   alpha0, var.alpha0, beta0, var.beta0,
                   cov.alpha0.beta0, var.nondiseased,
                   level)
  ##
  tcuts1 <- ciRegr(cuts.tr,
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
  #
  ocut <- transf(optcut, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  ce <- ciEllipse(ocut, 
                  alpha0, var.alpha0, beta0, var.beta0,
                  alpha1, var.alpha1, beta1, var.beta1,
                  cov.alpha0.beta0, cov.alpha1.beta1,
                  cov.alpha0.alpha1, cov.alpha0.beta1,
                  cov.alpha1.beta0, cov.beta0.beta1, 
                  var.nondiseased, var.diseased, level)
  
  ##
  if (ciSens) {
    xvals <- c(x.upper.Se, x.lower.Se[order(x.lower.Se)])
    yvals <- c(y.upper.Se, y.lower.Se[order(y.lower.Se)])
    #
    hull_indices <- chull(cbind(xvals, yvals))
    hull_indices <- c(hull_indices, hull_indices[1])
    #
    if (shading == "shade")
      polygon(xvals[hull_indices], yvals[hull_indices],
              col = rgb(0.5, 0.5, 0.5, alpha = 0.2), border = NA)
    else if (shading == "hatch")
      polygon(xvals[hull_indices], yvals[hull_indices],
              density = 20, angle = 90,
              col = col.hatching, border = col.hatching, lwd = lwd.hatching)
    #
    lines(x.upper.Se, y.upper.Se, col = col.ci, lwd = lwd.ci)
    lines(x.lower.Se, y.lower.Se, col = col.ci, lwd = lwd.ci)
  }
  ##
  if (ciSpec) {
    xvals <- c(x.upper.Sp, x.lower.Sp[order(y.lower.Sp)])
    yvals <- c(y.upper.Sp, y.lower.Sp[order(y.lower.Sp)])
    #
    hull_indices <- chull(cbind(xvals, yvals))
    hull_indices <- c(hull_indices, hull_indices[1])
    #
    if (shading == "shade")
      polygon(xvals[hull_indices], yvals[hull_indices],
              col = rgb(0.2, 0.2, 0.2, alpha = 0.2), border = NA)
    else if (shading == "hatch")
      polygon(xvals[hull_indices], yvals[hull_indices],
              density = 20, angle = 0,
              col = col.hatching, border = col.hatching, lwd = lwd.hatching)
    ##
    lines(x.upper.Sp, y.upper.Sp, col = col.ci, lwd = lwd.ci, lty = 2)
    lines(x.lower.Sp, y.lower.Sp, col = col.ci, lwd = lwd.ci, lty = 2)
  }
  ##
  ## Add ROC curve
  ##
  curve(pdiag(alpha1 + beta1 *
              (qdiag(1 - x, distr) - alpha0) / beta0,
              distr, FALSE),
        lwd = lwd, col = col, add = TRUE)
  ##
  if (mark.optcut)
    points(pdiag(alpha0 + beta0 * ocut, distr, FALSE),
           pdiag(alpha1 + beta1 * ocut, distr, FALSE),
           lwd = lwd.optcut, cex = 2, pch = 3, col = col.optcut)
  ##
  ## Add data
  ##
  points(1 - Spec, Sens,
         pch = pch.points, col = col.points, cex = cex)
  ##
  ## Add text
  ##
  if (mark.cutpoints) {
    text.cuts <- as.character(round(cuts, 2))
    ##
    points(pdiag(alpha0 + beta0 * cuts.tr, distr, FALSE),
           pdiag(alpha1 + beta1 * cuts.tr, distr, FALSE),
           pch = 3, cex = cex)
    ##
    text(1.02 - pdiag(alpha0 + beta0 * cuts.tr, distr),
         0.98 - pdiag(alpha1 + beta1 * cuts.tr, distr),
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
    for (t in seq_len(n)) {
      xx[t] <- pdiag(-ce$logit.spec +
                     q * ce$se.y0 * cos(t / 100 + acos(ce$r)), distr)
      yy[t] <- pdiag(ce$logit.sens +
                     q * ce$se.y1 * cos(t / 100), distr)
    }
    lines(xx, yy)
  }
  ##
  invisible(NULL)
}


##
## Biomarker distributions (densities)
##
density <- function(cutoff, Sens, Spec, studlab,
                    direction, distr, log.cutoff,
                    xlab, ylab, xlim, log.axis, ylim,
                    mains, which,
                    lines, lwd.study, rlines,
                    alpha0, var.alpha0, beta0, var.beta0,
                    alpha1, var.alpha1, beta1, var.beta1,
                    cov.alpha0.alpha1, cov.alpha0.beta0,
                    cov.alpha0.beta1, cov.alpha1.beta0,
                    cov.alpha1.beta1, cov.beta0.beta1,
                    var.nondiseased, var.diseased,
                    lambda,
                    col, lwd,
                    points, cex, col.points, pch.points,
                    ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                    min.cutoff, max.cutoff, mark.cutpoints,
                    optcut, mark.optcut,
                    line.optcut, col.optcut,
                    lwd.optcut,
                    cex.marks,
                    ellipse, shading,
                    col.hatching, lwd.hatching,
                    youden,
                    x,
                    ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  ymax <- max(c(beta0 * ddiag(alpha0 + beta0 * xvals.tr, distr),
                beta1 * ddiag(alpha1 + beta1 * xvals.tr, distr)))
  #
  plot(xvals, beta0 * ddiag(alpha0 + beta0 * xvals.tr, distr),
       xlim = xlim, ylim = if (is.null(ylim)) c(0, ymax) else ylim,
       las = 1, type = "l",
       xlab = xlab, ylab = "",
       main = mains[match("density", which)],
       lty = 2, col = col, lwd = lwd)
  #
  lines(xvals, beta1 * ddiag(alpha1 + beta1 * xvals.tr, distr),
        lty = 1, col = col, lwd = lwd)
  #
  # draw optimal cut-off
  #
  if (line.optcut)
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}


##
## Sensitivity and specificity
##
sensspec <- function(cutoff, Sens, Spec, studlab,
                     direction, distr, log.cutoff,
                     xlab, ylab, xlim, log.axis, ylim,
                     mains, which,
                     lines, lwd.study, rlines,
                     alpha0, var.alpha0, beta0, var.beta0,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha0.alpha1, cov.alpha0.beta0,
                     cov.alpha0.beta1, cov.alpha1.beta0,
                     cov.alpha1.beta1, cov.beta0.beta1,
                     var.nondiseased, var.diseased,
                     lambda,
                     col, lwd,
                     points, cex, col.points, pch.points,
                     ci, ciSens, ciSpec, level, col.ci, lwd.ci,
                     min.cutoff, max.cutoff, mark.cutpoints,
                     optcut, mark.optcut,
                     line.optcut, col.optcut,
                     lwd.optcut,
                     cex.marks,
                     ellipse, shading,
                     col.hatching, lwd.hatching,
                     youden,
                     x,
                     ...) {
  #
  if (is.null(xlim))
    xlim <- c(min.cutoff, max.cutoff)
  #
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  xvals.tr <- transf(xvals, direction, log.cutoff, min.cutoff, max.cutoff)
  #
  plot(c(cutoff, cutoff), c(Spec, Sens),
       type = "n",
       las = 1, log = log.axis,
       ylab = "Sensitivity / Specificity", xlab = xlab,
       main = mains[match("sensspec", which)],
       xlim = xlim, ylim = if (is.null(ylim)) c(0, 1) else ylim, ...)
  #
  # Add lines
  #   
  if (lines) {
    for (s in studlab) {
      lines(cutoff[studlab == s], Spec[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
      lines(cutoff[studlab == s], Sens[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    }
  }
  #
  # Add data
  #
  if (points) {
    points(cutoff, Sens, pch = pch.points, cex = cex, col = col.points)
    points(cutoff, Spec, pch = 1, cex = cex, col = col.points)
  }
  #
  # Add curves
  #
  if (rlines) {
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$TE,
                   distr),
          lty = 2, col = col, lwd = lwd)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$TE,
                   distr),
          lty = 1, col = col, lwd = lwd)
  }
  #
  if (ci) {
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$lower,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSpec(ciRegr(xvals.tr,
                          alpha0, var.alpha0,
                          beta0, var.beta0,
                          cov.alpha0.beta0,
                          var.nondiseased,
                          level)$upper,
                   distr),
          lty = 2, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$lower,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
    #
    lines(xvals,
          calcSens(ciRegr(xvals.tr,
                          alpha1, var.alpha1,
                          beta1, var.beta1,
                          cov.alpha1.beta1,
                          var.diseased,
                          level)$upper,
                   distr),
          lty = 1, col = col.ci, lwd = lwd.ci)
  }
  #
  # Draw line for optimal cutoff
  #
  if (line.optcut) 
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}
