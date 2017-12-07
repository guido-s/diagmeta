plot.diagmeta <- function(x,
                          which = c("survival", "youden", "roc", "sroc"),
                          xlab = "threshold",
                          ci = FALSE, ciSens = FALSE, ciSpec = FALSE,
                          mark.optcut = FALSE, mark.cutpoints = FALSE,
                          points = TRUE, lines = FALSE,
                          rlines = TRUE, line.optcut = TRUE,
                          col.points = "rainbow",
                          cex = 1, pch.points = 16,
                          cex.marks = 0.7 * cex,
                          lwd = 1, lwd.optcut = 2 * lwd,
                          shading = "none",
                          xlim = NULL,
                          ...) {
  
  
  meta:::chkclass(x, "diagmeta")
  ##
  setchar <- meta:::setchar
  chklogical <- meta:::chklogical
  ##
  plot.types <- c("regression", "cdf", "survival", "youden",
                  "roc", "sroc", "density")
  
  
  if (is.character(which))
    which <- setchar(which, plot.types)
  if (is.numeric(which)) {
    if (any(which < 1) | any(which > 7))
      stop("Numeric values of argument 'which' must be in 1:7")
    which <- plot.types[which]
  }
  which <- unique(which)
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
  
  
  study.no <- as.numeric(as.factor(x$studlab))
  ##  
  if (col.points == "rainbow")
    cols <- rainbow(x$k)[study.no]
  else if (col.points == "topo")
    cols <- topo.colors(x$k)[study.no]
  else if (col.points == "heat")
    cols <- heat.colors(x$k)[study.no]
  else if (col.points == "terrain")
    cols <- terrain.colors(x$k)[study.no]
  else if (col.points == "cm")
    cols <- cm.colors(x$k)[study.no]
  else if (col.points == "gray")
    cols <- gray(1:x$k / (x$k + 1))[study.no]
  else if (col.points == "black")
    cols <- rep(1, x$k)
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
  
  
  cutoff <- x$cutoff
  studlab <- x$studlab
  ##
  distr <- x$distr
  level <- x$level
  lambda <- x$lambda
  optcut <- x$optcut
  ##
  alpha0 <- x$regr$alpha0
  varalpha0 <- x$regr$varalpha0
  beta0 <- x$regr$beta0
  varbeta0 <- x$regr$varbeta0
  covalpha0beta0 <- x$regr$covalpha0beta0
  alpha1 <- x$regr$alpha1
  varalpha1 <- x$regr$varalpha1
  beta1 <- x$regr$beta1
  varbeta1 <- x$regr$varbeta1
  covalpha1beta1 <- x$regr$covalpha1beta1
  covalpha0alpha1 <- x$regr$covalpha0alpha1
  covalpha0beta1 <- x$regr$covalpha0beta1
  covalpha1beta0 <- x$regr$covalpha1beta0
  covbeta0beta1 <- x$regr$covbeta0beta1
  ##
  var.nondiseased <- x$var.nondiseased
  var.diseased <- x$var.diseased
  ##
  mean0 <- x$dist$mean0
  var.mean0 <- x$dist$var.mean0
  sd0 <- x$dist$sd0
  var.sd0 <- x$dist$var.sd0
  mean1 <- x$dist$mean1
  var.mean1 <- x$dist$var.mean1
  sd1 <- x$dist$sd1
  var.sd1 <- x$dist$var.sd1
  ##
  NN0 <- x$NN0
  NN1 <- x$NN1
  ##
  NN <- x$data.lmer$NN
  Cutoff <- x$data.lmer$Cutoff
  ##
  Se.optcut <- x$Se.optcut
  Sp.optcut <- x$Sp.optcut
  ##
  log.cutoff <- x$log.cutoff
  log.axis <- if (log.cutoff) "x" else ""
  ##
  youden <- calcYouden(1 - NN1, NN0, lambda)
  
  
  ##
  ##
  ## Linear regression lines in logit / probit space
  ##
  ##
  if ("regression" %in% which) {
    ##
    plot(c(cutoff, cutoff), rescale(c(NN0, NN1), distr),
         type = "n", las = 1, log = log.axis,
         ylab = ylab, xlab = xlab,
         xlim = xlim, ...)
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
        curve(ci.y(log(x),
                   alpha0, varalpha0, beta0, varbeta0, covalpha0beta0,
                   var.nondiseased,
                   level)$lower,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(log(x),
                   alpha0, varalpha0, beta0, varbeta0, covalpha0beta0,
                   var.nondiseased,
                   level)$upper,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(log(x),
                   alpha1, varalpha1, beta1, varbeta1, covalpha1beta1,
                   var.diseased,
                   level)$lower,
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(log(x),
                   alpha1, varalpha1, beta1, varbeta1, covalpha1beta1,
                   var.diseased,
                   level)$upper,
              lty = 1, col = 1, lwd = lwd, add = TRUE)
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
        curve(ci.y(x,
                   alpha0, varalpha0, beta0, varbeta0, covalpha0beta0,
                   var.nondiseased,
                   level)$lower,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(x,
                   alpha0, varalpha0, beta0, varbeta0, covalpha0beta0,
                   var.nondiseased,
                   level)$upper,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(x,
                   alpha1, varalpha1, beta1, varbeta1, covalpha1beta1,
                   var.diseased,
                   level)$lower,
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ci.y(x,
                   alpha1, varalpha1, beta1, varbeta1, covalpha1beta1,
                   var.diseased,
                   level)$upper,
              lty = 1, col = 1, lwd = lwd, add = TRUE)
      }
    }
    ##
    if (lines)
      for (s in studlab)
        lines(cutoff[which(studlab == s)],
              rescale(NN0[which(studlab == s)], distr),
              col = col.points[studlab == s], lwd = lwd, lty = 2)
    ##
    ## Add data
    ##
    points(cutoff, rescale(NN1, distr),
           pch = pch.points, cex = cex, col = col.points)      
    ##
    if (lines)
      for (s in studlab)
        lines(cutoff[which(studlab == s)],
              rescale(NN1[which(studlab == s)], distr),
              col = col.points[studlab == s], lwd = lwd)
    ##
    points(cutoff, rescale(NN0, distr),
           pch = 1, cex = cex, col = col.points)
  }
  
  
  ##
  ##
  ## Data and biomarker distributions functions (empirical
  ## distribution function)
  ##
  ##
  if ("cdf" %in% which) {
    ##
    plot(c(cutoff, cutoff), c(NN0, NN1),
         type = "n", las = 1, log = log.axis,
         xlab = xlab, ylab = "Prob(negative test)",
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add regression curves
    ##
    if (log.cutoff) {
      ##
      if (rlines) {
        curve(calcSpec(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$TE),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$TE),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$lower),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$upper),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$lower),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$upper),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      ##
      if (rlines) {
        curve(calcSpec(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$TE),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$TE),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSpec(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$lower),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$upper),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$lower),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSpec(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$upper),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
    ##
    ## Add lines
    ##
    if (lines)
      for (s in studlab)
        lines(cutoff[which(studlab == s)], 
              NN0[which(studlab == s)],
              col = col.points[studlab == s], lwd = lwd, lty = 2)
    ##
    ## Add data
    ##
    points(cutoff, NN1,
           pch = pch.points, cex = cex, col = col.points)
    ##
    if (lines)
      for (s in studlab)
        lines(cutoff[which(studlab == s)], 
              NN1[which(studlab == s)],
              col = col.points[studlab == s], lwd = lwd)
    ##
    points(cutoff, NN0,
           pch = 1, cex = cex, col = col.points)
  }
  
  
  ##
  ##
  ## Data and biomarker distributions functions (survival function)
  ##
  ##
  if ("survival" %in% which) {
    ##
    plot(c(cutoff, cutoff), 1 - c(NN0, NN1),
         type = "n", las = 1, log = log.axis,
         xlab = xlab, ylab = "Prob(positive test)",
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Add regression curves
    ##
    if (log.cutoff) {
      ##
      if (rlines) {
        curve(calcSens(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$TE),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$TE),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSens(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$lower),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(log(x),
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$upper),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$lower),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(log(x),
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$upper),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    else {
      ##
      if (rlines) {
        curve(calcSens(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$TE),
              lty = 1, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$TE),
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (ci) {
        curve(calcSens(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$lower),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(x,
                            alpha0, varalpha0,
                            beta0, varbeta0,
                            covalpha0beta0,
                            var.nondiseased,
                            level)$upper),
              lty = 1, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$lower),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
        ##
        curve(calcSens(ci.y(x,
                            alpha1, varalpha1,
                            beta1, varbeta1,
                            covalpha1beta1,
                            var.diseased,
                            level)$upper),
              lty = 2, col = "gray", lwd = lwd, add = TRUE)
      }
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
    ##
    ## Add lines
    ##
    if (lines) {
      for (s in studlab)
        lines(cutoff[which(studlab == s)], 
              1 - NN1[which(studlab == s)],
              col = col.points[studlab == s], lwd = lwd)
      ##
      for (s in studlab)
        lines(cutoff[which(studlab == s)], 
              1 - NN0[which(studlab == s)],
              col = col.points[studlab == s], lwd = lwd, lty = 2)
    }
    ##
    ## Add data
    ##
    points(cutoff, 1 - NN1,
           pch = pch.points, cex = cex, col = col.points)
    ##
    points(cutoff, 1 - NN0, pch = 1, cex = cex, col = col.points)
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
         xlim = xlim, ylim = c(0, 1), ...)
    ##
    ## Plot Youden index
    ##
    if (log.cutoff) {
      if (ci) {
        curve(ciYouden(log(x), distr, lambda,
                       alpha0, varalpha0, beta0, varbeta0,
                       covalpha0alpha1, covalpha0beta0, covalpha0beta1,
                       alpha1, varalpha1, beta1, varbeta1,
                       covalpha1beta0, covalpha1beta1, covbeta0beta1,
                       var.nondiseased, var.diseased,
                       level)$lower,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ciYouden(log(x), distr, lambda,
                       alpha0, varalpha0, beta0, varbeta0,
                       covalpha0alpha1, covalpha0beta0, covalpha0beta1,
                       alpha1, varalpha1, beta1, varbeta1,
                       covalpha1beta0, covalpha1beta1, covbeta0beta1,
                       var.nondiseased, var.diseased,
                       level)$upper,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (rlines)
        curve(calcYouden2(log(x), distr, lambda,
                          alpha0, beta0, alpha1, beta1,
                          mean0, sd0, mean1, sd1),
              col = 1, add = TRUE)     
    }
    else {
      if (ci) {
        curve(ciYouden(x, distr, lambda,
                       alpha0, varalpha0, beta0, varbeta0,
                       covalpha0alpha1, covalpha0beta0, covalpha0beta1,
                       alpha1, varalpha1, beta1, varbeta1,
                       covalpha1beta0, covalpha1beta1, covbeta0beta1,
                       var.nondiseased, var.diseased,
                       level)$lower,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
        ##
        curve(ciYouden(x, distr, lambda,
                       alpha0, varalpha0, beta0, varbeta0,
                       covalpha0alpha1, covalpha0beta0, covalpha0beta1,
                       alpha1, varalpha1, beta1, varbeta1,
                       covalpha1beta0, covalpha1beta1, covbeta0beta1,
                       var.nondiseased, var.diseased,
                       level)$upper,
              lty = 2, col = 1, lwd = lwd, add = TRUE)
      }
      ##
      if (rlines)
        curve(calcYouden2(x, distr, lambda,
                          alpha0, beta0, alpha1, beta1,
                          mean0, sd0, mean1, sd1),
              col = 1, add = TRUE)     
    }
    ##
    ## Draw line for optimal cutoff
    ##
    if (line.optcut)
      abline(v = optcut)
    ##
    ## Add lines
    ##
    if (lines) {
      for (s in studlab)
        lines(cutoff[which(studlab == s)],
              youden[which(studlab == s)],
              col = col.points[studlab == s], lwd = lwd)
    }
    ##
    ## Add data
    ##
    points(cutoff, youden,
           pch = pch.points, col = col.points, cex = cex)
  }
  
  
  ##
  ##
  ## Study-specific ROC curves
  ##
  ##
  if ("roc" %in% which) {
    ##
    plot(1 - NN0, 1 - NN1,
         type = "n", las = 1,
         xlab = "1 - Specificity", ylab = "Sensitivity",
         xlim = c(0, 1), ylim = c(0, 1), ...)
    ##
    ## Add lines
    ##
    for (s in studlab)
      lines(c(1, 1 - NN0[studlab == s], 0),
            c(1, 1 - NN1[studlab == s], 0),
            col = col.points[studlab == s], lwd = lwd)
    ##
    ## Add data
    ##
    if (points)
      points(1 - NN0, 1 - NN1,
             pch = pch.points, col = col.points, cex = cex)
  }
  
  
  ##
  ##
  ## SROC curve
  ##
  ##
  if ("sroc" %in% which) {
    ##
    plot(1 - NN0, 1 - NN1,
         type = "n", las = 1,
         xlab = "1 - Specificity", ylab = "Sensitivity",
         xlim = c(0, 1), ylim = c(0, 1), ...)
    ##
    if (log.cutoff)
      cuts <- log(unique(cutoff))
    else
      cuts <- unique(cutoff)
    ##
    cuts <- cuts[order(cuts)]
    ##
    tcuts0 <- ci.y(cuts,
                   alpha0, varalpha0, beta0, varbeta0,
                   covalpha0beta0, var.nondiseased,
                   level)
    ##
    tcuts1 <- ci.y(cuts,
                   alpha1, varalpha1, beta1, varbeta1,
                   covalpha1beta1, var.diseased,
                   level)
    ##
    Sens <- calcSens(tcuts1$TE, distr)
    lowerSens <- calcSens(tcuts1$lower, distr)
    upperSens <- calcSens(tcuts1$upper, distr)
    ##
    Spec <- calcSpec(tcuts0$TE, distr)
    lowerSpec <- calcSpec(tcuts0$lower, distr)
    upperSpec <- calcSpec(tcuts0$upper, distr)
    ##
    y.lower.Se <- lowerSens
    y.lower.Sp <- Sens
    ##
    y.upper.Se <- upperSens
    y.upper.Sp <- Sens
    ##
    x.lower.Se <- 1 - Spec
    x.lower.Sp <- 1 - lowerSpec
    ##
    x.upper.Se <- 1 - Spec
    x.upper.Sp <- 1 - upperSpec
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
    if (is.logistic) {
      curve(1 - expit(alpha1 + beta1 * (logit(1 - x) - alpha0) / beta0),
            lwd = lwd, col = 1, add = TRUE)
      ##
      if (mark.optcut)
        points(1 - expit(alpha0 + beta0 * optcut),
               1 - expit(alpha1 + beta1 * optcut),
               lwd = lwd.optcut, cex = 2, pch = 3, col = 1)
    }
    else {
      curve(1 - pnorm(qnorm(1 - x, mean0, sd0), mean1, sd1),
            lwd = lwd, col = 1, add = TRUE)
      ##
      if (mark.optcut)
        points(1 - pnorm(optcut, mean0, sd0),
               1 - pnorm(optcut, mean1, sd1),
               lwd = lwd.optcut, cex = 2, pch = 3, col = 1)
    }
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
      if (is.logistic) { 
        points(1 - expit(alpha0 + beta0 * cuts),
               1 - expit(alpha1 + beta1 * cuts),
               pch = 3, cex = cex)
        ##
        text(1.02 - expit(alpha0 + beta0 * cuts),
             0.98 - expit(alpha1 + beta1 * cuts),
             text.cuts, cex = cex.marks)
      }
      else {
        points(1 - pnorm(cuts, mean0, sd0),
               1 - pnorm(cuts, mean1, sd1),
               pch = 3, cex = cex)
        ##
        text(1.02 - pnorm(cuts, mean0, sd0),
             0.98 - pnorm(cuts, mean1, sd1),
             text.cuts, cex = cex.marks)
      }
    }
    ##
    ## Add data
    ##
    points(1 - NN0, 1 - NN1,
           pch = pch.points, col = col.points, cex = cex)

  }
  
  
  ##
  ##
  ## Biomarker distributions (densities)
  ##
  ##
  if ("density" %in% which) {
    ##
    if (is.logistic) {
      if (log.cutoff) {
        curve(beta0 * expit(beta0 * log(x) + alpha0) /
              (1 + exp(beta0 * log(x) + alpha0)),
              lty = 2, log = "x", las = 1,
              from = min(Cutoff), to = max(Cutoff),
              xlab = xlab, ylab = "")
        ##
        curve(beta1 * expit(beta1 * log(x) + alpha1) /
              (1 + exp(beta1 * log(x) + alpha1)),
              lty = 1, add = TRUE)
      }
      else {
        curve(beta0 * expit(beta0 * x + alpha0) /
              (1 + exp(beta0 * x + alpha0)),
              lty = 2, las = 1,
              from = min(Cutoff), to = max(Cutoff),
              xlab = xlab, ylab = "")
        ##
        curve(beta1 * expit(beta1 * x + alpha1) /
              (1 + exp(beta1 * x + alpha1)),
              lty = 1, add = TRUE)
      }
    }
    else if (is.normal) {
      if (log.cutoff) {
        curve(dnorm(log(x), mean0, sd0),
              lty = 2, log = "x",
              from = min(Cutoff), to = max(Cutoff),
              xlab = xlab)
        ##
        curve(dnorm(log(x), mean1, sd1),
              lty = 1, add = TRUE)
      }
      else {
        curve(dnorm(x, mean0, sd0),
              lty = 2,
              from = min(Cutoff), to = max(Cutoff),
              xlab = xlab, ...) 
        ##
        curve(dnorm(x, mean1, sd1),
              lty = 1, add = TRUE, ...)
      }
    }
    ##
    ## draw optimal cut-off
    ##
    if (line.optcut)
      abline(v = optcut)
  }
  
  
  invisible(NULL)
}
