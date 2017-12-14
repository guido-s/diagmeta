\name{diagmeta}
\alias{diagmeta}

\title{
Meta-analysis of diagnostic accuracy studies with the multiple cutoffs model
}

\description{
Diagnostic accuracy tests may be based on an ordinal or continuous
biomarker or an ordinal score together with a cutoff. The decision
whether the target condition is present or not (positive or negative
test result) depends on whether the observed value is above or below
the cutoff. Sensitivity and specificity of the test depend on the
chosen cutoff and vary with the cutoff. In meta-analysis of diagnostic
accuracy studies, results are often reported for multiple cutoffs
within a study, and the cutoffs may differ between studies. The
multiple cutoffs model creates a link between the range of cutoffs and
the respective pairs of sensitivity and specificity and thus allows
identifying cutoffs at which the test is likely to perform best
[Rücker 2016].
}

\usage{
diagmeta(TP, FP, TN, FN, cutoff, studlab, data = NULL, distr = "logistic",
         model = "DICS", equalvar = FALSE, lambda = 0.5, log.cutoff = FALSE,
         method.weights = "invvar", level = 0.95, incr = 0.5, n.iter.max = 1000,
         tol = 1e-08, silent = TRUE, ...)
}

\arguments{
  \item{TP, FP, TN, FN}{Numeric vectors giving the number of true
    positives, false positives, true negatives and false negatives}
  \item{cutoff}{A number vector indicating the cutoff values}
  \item{studlab}{A numeric or a character vector with study labels}
  \item{data}{An optional data frame containing the study information}
  \item{distr}{A character indicating the distribution(see details)}
  \item{model}{A character indicating the model (see Details)}
  \item{equalvar}{A logical indicating whether the variances of the
    biomarker in both groups are thought equal (see Details)}
  \item{lambda}{A numeric between 0 and 1 indicating the weight of the
    sensitivity}
  \item{log.cutoff}{A logical indicating whether the cutoffs should be
    log-transformed}
  \item{method.weights}{A character indicating the weighting method
    (see Details)}
   \item{incr}{A numeric between 0 and 1 that is added as a continuity
     correction}
   \item{level}{A numeric indicating the significance levels for
     tests} 
  \item{n.iter.max}{A numeric indicating the maximal number of fixed
    point iterations for finding the optimal cutoff}
  \item{tol}{A numeric indicating the tolerance for convergence of the
    fixed point iteration}
  \item{silent}{A logical indicating whether iterations should be
    printed}
  \item{...}{additional arguments}
}

\details{
Each row of the data set provides at least a study label, a cutoff and
the numbers of true positives, false positives, true negatives and
false negatives. Different studies may contribute a varying number of
cutoffs, as well as different sets of cutoffs.

The multiple cutoffs model is a multilevel random effects model. At
the study level, for the group of patients without the target
condition (in short disease-free), the specificities at all available
cutoffs together provide an estimate of the cumulative distribution
function (cdf) of the test results within the disease-free
individuals. Likewise, for patients with the target condition (in
short diseased), via the observed sensitivities at all observed
cutoffs we obtain an estimate of the cdf of the test results within
the diseased patients. At the meta-analytic level, the model fits the
data for both groups and all available cutoffs over all studies. Based
on a parametric model, it provides estimates of the two cdfs for the
two groups across all studies, accounting for the between-study
heterogeneity and correlation between groups.

Users have the choice between the normal (argument distr = “normal”)
and the logistic distribution (argument distr = “logistic” which is
the default). In addition, the cutoffs may be log-transformed
(argument log, default is FALSE).

The ‘pooled’ sensitivity and specificity values can be obtained at
every cutoff; and a multiple cutoffs summary ROC (sROC) naturally
follows while preserving cutoff information. The optimal cutoff is
defined as the cutoff where the maximum of a weighted sum of
sensitivity and specificity is obtained. The 95\% confidence intervals
(95\%CI) of sensitivities, specificities and the optimal cutoff are
estimated using the delta method [16].
}

\value{
An object of class “netmeta” with corresponding print, summary, and
plot function. The object is a list containing the following
components
  \item{TP}{As defined above}
  \item{FP}{As defined above}
  \item{TN}{As defined above}
  \item{FN}{As defined above}
  \item{cutoff}{As defined above}
  \item{studlab}{As defined above}
  \item{data}{As defined above}
  \item{distr}{As defined above}
  \item{model}{As defined above}
  \item{equalvar}{As defined above}
  \item{lambda}{As defined above}
  \item{log.cutoff}{As defined above}
  \item{method.weights}{As defined above}
  \item{k}{The number of studies in the meta-analysis}
  \item{optcut}{The optimal cutoff, defined as the point where lambda
    * sensitivity + (1 – lambda) * specificity is maximized}
  \item{Se.optcut}{The sensitivity at the optimal cutoff}
  \item{lower.Se.optcut}{The lower limit of the confidence interval
    (CI) of the sensitivity at the optimal cutoff}
  \item{upper.Se.optcut}{The upper limit of the confidence interval
    (CI) of the sensitivity at the optimal cutoff}
  \item{Sp.optcut}{The specificity at the optimal cutoff}
  \item{lower.Sp.optcut}{The lower limit of the confidence interval
    (CI) of the specificity at the optimal cutoff}
  \item{upper.Sp.optcut}{The upper limit of the confidence interval
    (CI) of the specificity at the optimal cutoff}
  \item{var.diseased}{The within-study variance for the diseased
    group} 
  \item{var.nondiseased}{The within-study variance for the
    non-diseased group}
  \item{NN0}{The total number of individuals in the diseased group}
  \item{NN1}{The total number of individuals in the non-diseased
    group}  
  \item{AIC}{The value of the Akaike information criterion of the lmer
    object}
  \item{BIC}{The value of the Bayesian information criterion of the
    lmer object}
  \item{..$ Study}{A vector of study labels}
  \item{..$ Group}{A vector of group labels (0 or 1)}
  \item{..$ Cutoff}{A vector of cutoffs}
  \item{..$ N}{A vector of group sizes}
  \item{..$ Negative}{A vector of number of negative test results}
  \item{..$ NN}{A vector of frequencies of negative test results}
  \item{result.lmer}{An object of class \code{\link[lme4]{lmer}}
    with 13 slots}
}

\references{
Steinhauser S, Schumacher M and Rücker G (2016),
Modelling multiple thresholds in meta-analysis of diagnostic test accuracy
studies.
\emph{BMC Med Res Methodol}, \bold{16}, 97.
}

\author{
Susanne Steinhauser \url{susanne.steinhauser@uni-koeln.de}, Srinath Kolampally  \url{kolampal@imbi.uni-freiburg.de}, Gerta Rücker \url{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \url{sc@imbi.uni-freiburg.de}
}


\seealso{
\code{\link{plot.diagmeta}}}
%%, \code{\link{print.diagmeta}},\code{\link{summary.print.diagmeta}},\code{\link{lmer}}}

\examples{
# Fenomeda data
#
data(fenomeda)

str(fenomeda)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint, study,
                  data = fenomeda, log.cutoff = TRUE,
                  model = "DIDS", distr = "l")

plot(diag1)
}
