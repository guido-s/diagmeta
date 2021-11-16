#' @importFrom meta ci gs
#' @importFrom stats dlogis dnorm plogis pnorm qlogis qnorm
#' @importFrom utils packageDescription


.onAttach <- function (libname, pkgname) {
  msg <- paste("Loading 'diagmeta' package (version ",
               utils::packageDescription("diagmeta")$Version,
               ").",
               "\nType 'help(\"diagmeta-package\")' for a brief overview.",
               sep = "")
  packageStartupMessage(msg)
}


## Quantil function
##
qdiag <- function (x, distr, lower.tail = TRUE) {
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic") 
    res <- qlogis(x, lower.tail = lower.tail)
  else if (distr == "normal") 
    res <- qnorm(x, lower.tail = lower.tail)
  ##
  res
}


## Cumulative distribution function
##
pdiag <- function (x, distr, lower.tail = TRUE) {
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic") 
    res <- plogis(x, lower.tail = lower.tail)
  else if (distr == "normal") 
    res <- pnorm(x, lower.tail = lower.tail)
  ##
  res
}


## Density function
##
ddiag <- function (x, distr) {
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic") 
    res <- dlogis(x)
  else if (distr == "normal") 
    res <- dnorm(x)
  ##
  res
}


ciRegr <- function(x,
                   alpha, var.alpha, beta, var.beta,
                   cov.alpha.beta, var.group,
                   level = 0.95) {
  ##
  y <- alpha + beta * x
  se.y <- sqrt(var.alpha + x^2 * var.beta + 2 * x * cov.alpha.beta + var.group)
  res <- ci(y, se.y, level)
  ##
  res
}


ciEllipse <- function(x,
                      alpha0, var.alpha0, beta0, var.beta0,
                      alpha1, var.alpha1, beta1, var.beta1,
                      cov.alpha0.beta0, cov.alpha1.beta1,
                      cov.alpha0.alpha1, cov.alpha0.beta1, cov.alpha1.beta0, cov.beta0.beta1, 
                      var.group0, var.group1,
                      level = 0.95) {
  ##
  y0 <- alpha0 + beta0 * x
  y1 <- alpha1 + beta1 * x
  se.y0 <- sqrt(var.alpha0 + x^2 * var.beta0 + 2 * x * cov.alpha0.beta0 + var.group0)
  se.y1 <- sqrt(var.alpha1 + x^2 * var.beta1 + 2 * x * cov.alpha1.beta1 + var.group1)
  covar <- cov.alpha0.alpha1 + x * cov.alpha0.beta1 + x *cov.alpha1.beta0 + x^2 *cov.beta0.beta1
  r <- covar/se.y0/se.y1
  res <- list(logit.spec = y0, logit.sens = -y1, se.y0 = se.y0, se.y1 = se.y1, covar = covar, r = r)
  ##
  res
}


calcSens <- function(x, distr = "logistic") {
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  ##
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  res <- pdiag(x, distr, FALSE)
  ##
  res
}


calcSpec <- function(x, distr = "logistic") {
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  ##
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  res <- pdiag(x, distr)
  ##
  res
}


## Calculate Youden index (from parameters)
##
calcYouden <- function(x, distr, lambda,
                       alpha0, beta0, alpha1, beta1) {
  ##
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  res <- 2 * (1 - lambda) * pdiag(alpha0 + beta0 * x, distr) +
    2 * lambda * pdiag(alpha1 + beta1 * x, distr, FALSE) - 1
  ##
  res
}


## Function to compute the variance and confidence intervals of the
## Youden index (normal distribtion)
##
ciYouden <- function(x, distr, lambda,
                     alpha0, var.alpha0, beta0, var.beta0,
                     cov.alpha0.alpha1, cov.alpha0.beta0, cov.alpha0.beta1,
                     alpha1, var.alpha1, beta1, var.beta1,
                     cov.alpha1.beta0, cov.alpha1.beta1, cov.beta0.beta1,
                     var.nondiseased, var.diseased,
                     level = 0.95) {
  
  distr <- setchar(distr, c("logistic", "normal"))
  ##
  youden <- calcYouden.SeSp(calcSens(alpha1 + beta1 * x, distr),
                            calcSpec(alpha0 + beta0 * x, distr),
                            lambda)
  ##
  p <- 2 * lambda       * ddiag(alpha1 + beta1 * x, distr)
  q <- 2 * (1 - lambda) * ddiag(alpha0 + beta0 * x, distr)
  ##
  var.youden <- q^2 *
    (var.alpha0 + x^2 * var.beta0 + 2 * x * cov.alpha0.beta0 + var.nondiseased) +
    p^2 * (var.alpha1 + x^2 * var.beta1 + 2 * x * cov.alpha1.beta1 + var.diseased) -
    2 * p * q * (cov.alpha0.alpha1 + x * cov.alpha0.beta1 +
                   x * cov.alpha1.beta0 + x^2 * cov.beta0.beta1)
  ##
  res <- ci(youden, sqrt(var.youden), level)
  ##
  res
}


## Calculate Youden index (from sensitivity and specificity)
##
calcYouden.SeSp <- function(sens, spec, lambda)
  2 * lambda * sens + 2 * (1 - lambda) * spec - 1


## Function for fixed point iteration
##
iterate <- function(f, x0, nmax, tol, print = FALSE) {
  x <- x0
  n <- 0
  while (abs(f(x) - x) > tol & n < nmax) {
    n <- n + 1
    x <- f(x)
    if (print) print(x)
  }
  if(n < nmax)
    return(x)
  else 
    stop("Iteration of maximal Youden index reached number of maximal iterations without converging.")
}


trapz <- function (x, y) {
  ##
  ## R function taken R package caTools, version 1.17.1
  ##
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1])) / 2)
}
