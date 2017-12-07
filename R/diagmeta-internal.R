.onAttach <- function (libname, pkgname) {
  msg <- paste("Loading 'diagmeta' package (version ",
               utils::packageDescription("diagmeta")$Version,
               ").",
               "\nType 'help(diagmeta-package)' for a brief overview.",
               sep = "")
  packageStartupMessage(msg)
}


## Define logit function
##
logit <- function(x)
  log(x) - log(1 - x)


## Define expit function
##
expit <- function(x)
  (1 + exp(-x))^(-1)


ci.y <- function(x,
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


calcSens <- function(x, distr = "logistic") {
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  ##
  distr <- meta:::setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic")
    x <- 1 - expit(x)
  else
    x <- 1 - pnorm(x)
  ##
  x
}


calcSpec <- function(x, distr = "logistic") {
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  ##
  distr <- meta:::setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic")
    x <- expit(x)
  else
    x <- pnorm(x)
  ##
  x
}


## Calculate Youden index
##
calcYouden <- function(sens, spec, lambda)
  2 * lambda * sens + 2 * (1 - lambda) * spec - 1


calcYouden2 <- function(x, distr, lambda,
                        alpha0, beta0, alpha1, beta1,
                        mean0, sd0, mean1, sd1) {
  ##
  distr <- meta:::setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic")
    res <- 2 * (1 - lambda) * expit(alpha0 + beta0 * x) +
      2 * lambda * (1 - expit(alpha1 + beta1 * x)) - 1
  else
    res <- 2 * (1 - lambda) * pnorm(log(x), mean0, sd0) +
      2 * lambda * (1 - pnorm(log(x), mean1, sd1)) - 1
  ##
  res
}


## Function to transform / rescale the x-values, either with logit
## oder probit function.
##
rescale <- function(x, distr) {
  if (distr == "logistic")
    res <- logit(x)
  else if (distr == "normal")
    res <- qnorm(x)
  ##
  res
}


## Function for fixed point iteration
##
iterate <- function(f, x0, nmax, tol, print = FALSE) {
  x <- x0
  n <- 0
  while (abs(f(x) - x) > tol & n < nmax) {
    n <- n+1
    x <- f(x)
    if (print) print(x)
  }
  if(n < nmax)
    return(x)
  else 
    stop("Iteration of maximal Youden index reached number of maximal iterations without converging.")
}


## Function to compute the variance and confidence intervals of the
## Youden index (normal distribtion)
##
ciYouden <- function(x, distr, lambda,
                     alpha0, varalpha0, beta0, varbeta0,
                     covalpha0alpha1, covalpha0beta0, covalpha0beta1,
                     alpha1, varalpha1, beta1, varbeta1,
                     covalpha1beta0, covalpha1beta1, covbeta0beta1,
                     var.nondiseased, var.diseased,
                     level = 0.95) {
  
  sp <- calcSpec(alpha0 + beta0 * x, distr)
  se <- calcSens(alpha1 + beta1 * x, distr)
  youden <- calcYouden(se, sp, lambda)
  ##
  q <- 2 * (1 - lambda) * sp * (1 - sp)
  ##
  distr <- meta:::setchar(distr, c("logistic", "normal"))
  ##
  if (distr == "logistic")
    p <- 2 * lambda * se * (1 - se)
  else
    p <- 2 * lambda * dnorm(beta1 * x + alpha1)
  ##
  var.youden <- q^2 *
                    (varalpha0 + x^2 * varbeta0 + 2 * x * covalpha0beta0 + var.nondiseased) +
                    p^2 * (varalpha1 + x^2 * varbeta1 + 2 * x * covalpha1beta1 + var.diseased) -
                        2 * p * q * (covalpha0alpha1 + x * covalpha0beta1 +
                                     x * covalpha1beta0 + x^2 * covbeta0beta1)
  ##
  se.youden <- sqrt(var.youden)
  ##
  res <- ci(youden, se.youden, level)
  ##
  res
}
