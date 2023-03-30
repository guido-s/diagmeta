##
## Linear regression lines in logit / probit space
##
regression <- function(cutoff, Sens, Spec, studlab,
                       distr, log.cutoff,
                       xlab, ylab, xlim, log.axis,
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
                       Cutoff, mark.cutpoints,
                       optcut, mark.optcut,
                       line.optcut, col.optcut,
                       lwd.optcut,
                       cex.marks,
                       ellipse, shading,
                       col.hatching, lwd.hatching,
                       youden,
                       x,
                       ...) {
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
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
      ##
      lines(cutoff[studlab == s],
            qdiag(1 - Sens[studlab == s], distr),
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    }
  ##
  if (points) {
    points(cutoff, qdiag(Spec, distr),
           pch = 1, cex = cex, col = col.points)
    ##
    points(cutoff, qdiag(1 - Sens, distr),
           pch = pch.points, cex = cex, col = col.points)
  }
  ##
  ## Add linear regression lines
  ##
  if (log.cutoff) {
    if (rlines) {
      curve(alpha0 + beta0 * log(x),
            lty = 2, col = col, lwd = lwd, add = TRUE)
      curve(alpha1 + beta1 * log(x),
            lty = 1, col = col, lwd = lwd, add = TRUE)
    }
    ##
    if (ci) {
      ##
      curve(ciRegr(log(x),
                   alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                   var.nondiseased,
                   level)$lower,
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(log(x),
                   alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                   var.nondiseased,
                   level)$upper,
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(log(x),
                   alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                   var.diseased,
                   level)$lower,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(log(x),
                   alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                   var.diseased,
                   level)$upper,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
  }
  else {
    if (rlines) {
      curve(alpha0 + beta0 * x,
            lty = 2, col = col, lwd = lwd, add = TRUE)
      curve(alpha1 + beta1 * x,
            lty = 1, col = col, lwd = lwd, add = TRUE)
    }
    ##
    if (ci) {
      ##
      curve(ciRegr(x,
                   alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                   var.nondiseased,
                   level)$lower,
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(x,
                   alpha0, var.alpha0, beta0, var.beta0, cov.alpha0.beta0,
                   var.nondiseased,
                   level)$upper,
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(x,
                   alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                   var.diseased,
                   level)$lower,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciRegr(x,
                   alpha1, var.alpha1, beta1, var.beta1, cov.alpha1.beta1,
                   var.diseased,
                   level)$upper,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
  }
  ##
  invisible(NULL)
}


##
## Data and biomarker distributions functions (empirical
## distribution function)
##
cdf <- function(cutoff, Sens, Spec, studlab,
                distr, log.cutoff,
                xlab, ylab, xlim, log.axis,
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
                Cutoff, mark.cutpoints,
                optcut, mark.optcut,
                line.optcut, col.optcut,
                lwd.optcut,
                cex.marks,
                ellipse, shading,
                col.hatching, lwd.hatching,
                youden,
                x,
                ...) {
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
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
      ##
      lines(cutoff[studlab == s], 
            1 - Sens[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    }
  ##
  ## Add data
  ##
  if (points) {
    points(cutoff, 1 - Sens,
           pch = pch.points, cex = cex, col = col.points)
    ##
    points(cutoff, Spec,
           pch = 1, cex = cex, col = col.points)
  }
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSpec(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(log(x),
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSpec(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(x,
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
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
                     distr, log.cutoff,
                     xlab, ylab, xlim, log.axis,
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
                     Cutoff, mark.cutpoints,
                     optcut, mark.optcut,
                     line.optcut, col.optcut,
                     lwd.optcut,
                     cex.marks,
                     ellipse, shading,
                     col.hatching, lwd.hatching,
                     youden,
                     x,
                     ...) {
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
            col = col.points[studlab == s], lwd = lwd.study, lty = 1)
    ##
    for (s in studlab)
      lines(cutoff[studlab == s], 
            1 - Spec[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study, lty = 2)
  }
  ##
  ## Add data
  ##
  if (points) {
    points(cutoff, Sens,
           pch = pch.points, cex = cex, col = col.points)
    ##
    if (points)
      points(cutoff, 1 - Spec, pch = 1, cex = cex, col = col.points)
  }
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
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
                   distr, log.cutoff,
                   xlab, ylab, xlim, log.axis,
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
                   Cutoff, mark.cutpoints,
                   optcut, mark.optcut,
                   line.optcut, col.optcut,
                   lwd.optcut,
                   cex.marks,
                   ellipse, shading,
                   col.hatching, lwd.hatching,
                   youden,
                   x,
                   ...) {
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
            col = col.points[studlab == s], lwd = lwd.study)
  }
  ##
  ## Add data
  ##
  if (points)
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
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciYouden(log(x), distr, lambda,
                     alpha0, var.alpha0, beta0, var.beta0,
                     cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                     var.nondiseased, var.diseased,
                     level)$upper,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
    ##
    if (rlines)
      curve(calcYouden(log(x), distr, lambda,
                       alpha0, beta0, alpha1, beta1),
            col = col, lwd = lwd, add = TRUE)     
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
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(ciYouden(x, distr, lambda,
                     alpha0, var.alpha0, beta0, var.beta0,
                     cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                     var.nondiseased, var.diseased,
                     level)$upper,
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
    ##
    if (rlines)
      curve(calcYouden(x, distr, lambda,
                       alpha0, beta0, alpha1, beta1),
            col = col, add = TRUE)     
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
## Study-specific ROC curves
##
roc <- function(cutoff, Sens, Spec, studlab,
                distr, log.cutoff,
                xlab, ylab, xlim, log.axis,
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
                Cutoff, mark.cutpoints,
                optcut, mark.optcut,
                line.optcut, col.optcut,
                lwd.optcut,
                cex.marks,
                ellipse, shading,
                col.hatching, lwd.hatching,
                youden,
                x,
                ...) {
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
          col = col.points[studlab == s], lwd = lwd.study)
  ##
  ## Add data
  ##
  if (points)
    points(1 - Spec, Sens,
           pch = pch.points, col = col.points, cex = cex)
  ##
  invisible(NULL)
}


##
## SROC curve
##
sroc <- function(cutoff, Sens, Spec, studlab,
                 distr, log.cutoff,
                 xlab, ylab, xlim, log.axis,
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
                 Cutoff, mark.cutpoints,
                 optcut, mark.optcut,
                 line.optcut, col.optcut,
                 lwd.optcut,
                 cex.marks,
                 ellipse, shading,
                 col.hatching, lwd.hatching,
                 youden,
                 x,
                 ...) {
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
              density = 20, angle = 90,
              col = col.hatching, border = col.hatching, lwd = lwd.hatching)
    ##
    lines(x.upper.Se, y.upper.Se, col = col.ci, lwd = lwd.ci)
    lines(x.lower.Se, y.lower.Se, col = col.ci, lwd = lwd.ci)
  }
  ##
  if (ciSpec) {
    if (shading == "shade")
      polygon(c(x.upper.Sp, x.lower.Sp[order(y.lower.Sp)]),
              c(y.upper.Sp, y.lower.Sp[order(y.lower.Sp)]),
              col = rgb(0.2, 0.2, 0.2, alpha = 0.2), border = NA)
    else if (shading == "hatch")
      polygon(c(x.upper.Sp, x.lower.Sp[order(y.lower.Sp)]),
              c(y.upper.Sp, y.lower.Sp[order(y.lower.Sp)]),
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
  if (mark.optcut) {
    if (log.cutoff)
      ocut <- log(optcut)
    else
      ocut <- optcut
    ##
    points(pdiag(alpha0 + beta0 * ocut, distr, FALSE),
           pdiag(alpha1 + beta1 * ocut, distr, FALSE),
           lwd = lwd.optcut, cex = 2, pch = 3, col = col.optcut)
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
  ##
  invisible(NULL)
}


##
## Biomarker distributions (densities)
##
density <- function(cutoff, Sens, Spec, studlab,
                    distr, log.cutoff,
                    xlab, ylab, xlim, log.axis,
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
                    Cutoff, mark.cutpoints,
                    optcut, mark.optcut,
                    line.optcut, col.optcut,
                    lwd.optcut,
                    cex.marks,
                    ellipse, shading,
                    col.hatching, lwd.hatching,
                    youden,
                    x,
                    ...) {
  ##
  if (log.cutoff) {
    curve(beta0 * ddiag(beta0 * log(x) + alpha0, distr),
          log = "x", las = 1,
          from = min(exp(Cutoff)), to = max(exp(Cutoff)),
          xlab = xlab, ylab = "",
          main = mains[match("density", which)],
          lty = 2, col = col, lwd = lwd)
    ##
    curve(beta1 * ddiag(beta1 * log(x) + alpha1, distr),
          lty = 1, col = col, lwd = lwd, add = TRUE)
  }
  else {
    curve(beta0 * ddiag(beta0 * x + alpha0, distr),
          las = 1,
          from = min(Cutoff), to = max(Cutoff),
          xlab = xlab, ylab = "",
          main = mains[match("density", which)],
          lty = 2, col = col, lwd = lwd)
    ##
    curve(beta1 * ddiag(beta1 * x + alpha1, distr),
          lty = 1, col = col, lwd = lwd, add = TRUE)
  }
  ##
  ## draw optimal cut-off
  ##
  if (line.optcut)
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}


##
## Sensitivity and specificity
##
sensspec <- function(cutoff, Sens, Spec, studlab,
                     distr, log.cutoff,
                     xlab, ylab, xlim, log.axis,
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
                     Cutoff, mark.cutpoints,
                     optcut, mark.optcut,
                     line.optcut, col.optcut,
                     lwd.optcut,
                     cex.marks,
                     ellipse, shading,
                     col.hatching, lwd.hatching,
                     youden,
                     x,
                     ...) {
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
            col = col.points[studlab == s], lwd = lwd.study,
            lty = 2)
      lines(cutoff[studlab == s], Sens[studlab == s],
            col = col.points[studlab == s], lwd = lwd.study,
            lty = 1)
    }
  }
  ##
  ## Add data
  ##
  if (points) {
    points(cutoff, Sens, pch = pch.points, cex = cex,
           col = col.points)
    points(cutoff, Spec, pch = 1, cex = cex, col = col.points)
  }
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(log(x),
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(log(x),
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
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
            lty = 2, col = col, lwd = lwd, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$TE,
                     distr),
            lty = 1, col = col, lwd = lwd, add = TRUE)
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
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSpec(ciRegr(x,
                            alpha0, var.alpha0,
                            beta0, var.beta0,
                            cov.alpha0.beta0,
                            var.nondiseased,
                            level)$upper,
                     distr),
            lty = 2, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$lower,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
      ##
      curve(calcSens(ciRegr(x,
                            alpha1, var.alpha1,
                            beta1, var.beta1,
                            cov.alpha1.beta1,
                            var.diseased,
                            level)$upper,
                     distr),
            lty = 1, col = col.ci, lwd = lwd.ci, add = TRUE)
    }
  }
  ##
  ## Draw line for optimal cutoff
  ##
  if (line.optcut) 
    abline(v = optcut, col = col.optcut, lwd = lwd.optcut)
  ##
  invisible(NULL)
}  
