## This function is to extract values of sensitivity, specificity
## (with confidence intervals), PPV, NPV, PD (probability of disease) 
## and AUC, as well as the probability functions (densities) 
## for a given model (r), cutoff (x) and prevalence (pr). 
## r is obligatory, x and pr can be NA.

## Arguments: r = a diagmeta object, pr = prevalence, x = cutoff

values <- function(r = NULL, x = NA, pr = NA) {
  if (!inherits(r, "diagmeta")) 
    stop("Argument 'r' must be an object of class \"diagmeta\"")
  log <- r@log
  cf <- r@regressionCoefficients
  v1 <- r@v.within1
  v0 <- r@v.within0
  if (!r@logistic) {
    if (!log) {
      se <- 1 - pnorm(cf$beta1 * x + cf$alpha1)
      sp <- pnorm(cf$beta0 * x + cf$alpha0)
      SE.se <- sqrt(cf$varalpha1 + x^2 * cf$varbeta1 + 2 * x * cf$covalpha1beta1 + v1)
      CI.se <- c(se, 1-pnorm(cf$beta1 * x + cf$alpha1 + 1.96 * SE.se),
                 1-pnorm(cf$beta1 * x + cf$alpha1 - 1.96 * SE.se))
      SE.sp <- sqrt(cf$varalpha0 + x^2 * cf$varbeta0 + 2 * x * cf$covalpha0beta0 + v0)
      CI.sp <- c(sp, pnorm(cf$beta0 * x + cf$alpha0 - 1.96 * SE.sp),
                 pnorm(cf$beta0 * x + cf$alpha0 + 1.96 * SE.sp))
      dens.g <- cf$beta0  *  dnorm(cf$beta0 * x + cf$alpha0)
      dens.f <- cf$beta1  *  dnorm(cf$beta1 * x + cf$alpha1)
    }
    if (log) {
      se <- 1 - pnorm(cf$beta1 * log(x) + cf$alpha1)
      sp <- pnorm(cf$beta0 * log(x) + cf$alpha0)
      SE.se <- sqrt(cf$varalpha1 + log(x)^2 * cf$varbeta1 + 2 * log(x) * cf$covalpha1beta1 + v1)
      CI.se <- c(se, 1-pnorm(cf$beta1 * log(x) + cf$alpha1 + 1.96 * SE.se),
                 1-pnorm(cf$beta1 * log(x) + cf$alpha1 - 1.96 * SE.se))
      SE.sp <- sqrt(cf$varalpha0 + log(x)^2 * cf$varbeta0 + 2 * log(x) * cf$covalpha0beta0 + v0)
      CI.sp <- c(sp, pnorm(cf$beta0 * log(x) + cf$alpha0 - 1.96 * SE.sp),
                 pnorm(cf$beta0 * log(x) + cf$alpha0 + 1.96 * SE.sp))
      dens.g <- cf$beta0  *  dnorm(cf$beta0 * log(x) + cf$alpha0)
      dens.f <- cf$beta1  *  dnorm(cf$beta1 * log(x) + cf$alpha1)
    }
  }
  if (r@logistic) {
    if (!log) {
      se <- 1 - expit(cf$beta1 * x + cf$alpha1)
      sp <- expit(cf$beta0 * x + cf$alpha0)
      SE.se <- sqrt(cf$varalpha1 + x^2 * cf$varbeta1 + 2 * x * cf$covalpha1beta1 + v1)
      CI.se <- c(se, 1-expit(cf$beta1 * x + cf$alpha1 + 1.96 * SE.se),
                 1-expit(cf$beta1 * x + cf$alpha1 - 1.96 * SE.se))
      SE.sp <- sqrt(cf$varalpha0 + x^2 * cf$varbeta0 + 2 * x * cf$covalpha0beta0 + v0)
      CI.sp <- c(sp, expit(cf$beta0 * x + cf$alpha0 - 1.96 * SE.sp),
                 expit(cf$beta0 * x + cf$alpha0 + 1.96 * SE.sp))
    }
    if (log) {
      se <- 1 - expit(cf$beta1 * log(x) + cf$alpha1)
      sp <- expit(cf$beta0 * log(x) + cf$alpha0)
      SE.se <- sqrt(cf$varalpha1 + log(x)^2 * cf$varbeta1 + 2 * log(x) * cf$covalpha1beta1 + v1)
      CI.se <- c(se, 1-expit(cf$beta1 * log(x) + cf$alpha1 + 1.96 * SE.se),
                 1-expit(cf$beta1 * log(x) + cf$alpha1 - 1.96 * SE.se))
      SE.sp <- sqrt(cf$varalpha0 + log(x)^2 * cf$varbeta0 + 2 * log(x) * cf$covalpha0beta0 + v0)
      CI.sp <- c(sp, expit(cf$beta0 * log(x) + cf$alpha0 - 1.96 * SE.sp),
                 expit(cf$beta0 * log(x) + cf$alpha0 + 1.96 * SE.sp))
    }
    dens.g <- cf$beta0  *  sp  *  (1 - sp)
    dens.f <- cf$beta1  *  se  *  (1 - se)
  }
  ## PPV, NPV and probability of disease (depending on x, because se and sp depend on x)
  dp <- se  *  pr + (1 - sp)  *  (1 - pr)
  dn <- (1 - pr)  *  sp  + pr  *  (1 - se)
  PPV <- se  *  pr / dp
  NPV <- sp  *  (1 - pr) / dn
  dd <- pr  *  dens.f + (1 - pr)  *  dens.g
  PD <- pr  *  dens.f / dd
  ## AUC with CI (not depending on whether x or log(x) is modelled)
  ## Grid for numerical integration
  t <- (1:9999)/10000
  a0 <- (cf$beta1 * cf$alpha0 - cf$beta0 * cf$alpha1)/cf$beta0 ## a0 intercept on the transformed scale
  b0 <- cf$beta1/cf$beta0                                  ## b0 slope on the transformed scale
  if (!r@logistic) {
    AUC <- trapz(t, pnorm(a0 + b0 * qnorm(t)))
  }
  if (r@logistic) {
    AUC <- trapz(t, expit(a0 + b0 * logit(t)))
  }
  ## Result
  res <- list(cut = x, pr = pr,
              se = se, SE.se = SE.se, CI.se = CI.se,
              sp = sp, SE.sp = SE.sp, CI.sp = CI.sp,
              dens.g = dens.g, dens.f = dens.f,
              PPV = PPV, NPV = NPV, PD = PD, AUC = AUC)


  res
}
