\name{diagmeta}
\alias{diagmeta}

\title{
Meta-analysis of diagnostic test accuracy studies with the multiple cutoffs model
}

\description{ Diagnostic tests may be based on an ordinal or
continuous biomarker or an ordinal score together with a cutoff. The
decision whether the target condition is present or not (positive or
negative test result) depends on whether the observed value is above
or below the cutoff. Sensitivity and specificity of the test depend on
the chosen cutoff and vary with the cutoff. In meta-analysis of
diagnostic accuracy studies, results are often reported for multiple
cutoffs within a study, and the cutoffs may differ between
studies. The multiple cutoffs model creates a link between the range
of cutoffs and the respective pairs of sensitivity and specificity and
thus allows identifying cutoffs at which the test is likely to perform
best (Steinhauser et al., 2016).
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
  \item{distr}{A character indicating the distribution(see Details)}
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

Users have the choice between the normal (argument
\code{distr="normal"}) and the logistic distribution (argument
\code{distr="logistic"} which is the default). In addition, the
cutoffs may be log-transformed (argument \code{log.cutoff}, default is
\code{FALSE}).

The cdf, transformed using the quantile function of the chosen
distribution, is modelled by one of eight mixed linear models ("DIDS",
"CIDS", "DICS", "CICS", "DS", "CS", "DI", "CI") as described in
Steinhauser et al. (2016). The argument \code{equalvar} indicates if
the variances of the biomarker in both groups are assumed to be equal
(equalvar = TRUE) or unequal (equalvar = FALSE).

The ‘pooled’ sensitivity and specificity values can be obtained at
every cutoff; a multiple cutoffs summary ROC (sROC) naturally follows
while preserving cutoff information. The optimal cutoff is defined as
the cutoff where the maximum of a weighted sum of sensitivity and
specificity is obtained. The 95\% confidence intervals of
sensitivities, specificities and the optimal cutoff are estimated
using the delta method (Steinhauser et al., 2016).
}

\value{
An object of class "diagmeta" with corresponding print, summary, and
plot function. The object is a list containing the following
components
  \item{TP}{As defined above.}
  \item{FP}{As defined above.}
  \item{TN}{As defined above.}
  \item{FN}{As defined above.}
  \item{cutoff}{As defined above.}
  \item{studlab}{As defined above.}
  \item{Sens}{Sensitivity (original data).}
  \item{Spec}{Specificity (original data).}  
  \item{data}{As defined above.}
  \item{distr}{As defined above.}
  \item{model}{As defined above.}
  \item{equalvar}{As defined above.}
  \item{lambda}{As defined above.}
  \item{log.cutoff}{As defined above.}
  \item{method.weights}{As defined above.}
  \item{k}{The number of studies in the meta-analysis.}
  \item{optcut}{The optimal cutoff, defined as the point where lambda
    * sensitivity + (1 – lambda) * specificity is maximized.}
  \item{Sens.optcut}{The sensitivity at the optimal cutoff.}
  \item{lower.Sens.optcut, upper.Sens.optcut}{The lower and upper
    limit of the confidence interval of the sensitivity at the optimal
    cutoff.}
  \item{Spec.optcut}{The specificity at the optimal cutoff.}
  \item{lower.Spec.optcut, upper.Spec.optcut}{The lower and upper
    limit of the confidence interval of the specificity at the optimal
    cutoff.}
  \item{var.diseased, var.nondiseased}{The within-study variance for
    the diseased and non-diseased group, respectively.}
  \item{AIC}{The value of the Akaike information criterion of the lmer
    object.}
  \item{BIC}{The value of the Bayesian information criterion of the
    lmer object.}
  \item{data.lmer}{A list with objects Study (study labels), Group
    (group labels (0 or 1)), Cutoff, N (group sizes), Negative (number
    of negative test results), NN (frequencies of negative test
    results).}
  \item{result.lmer}{An object of class \code{\link[lme4]{lmer}}.}
  \item{n.iter.max}{As defined above.}
  \item{call}{Function call.}
  \item{version}{Version of R package \bold{diagmeta} used to create object.}
}

\references{
Steinhauser S, Schumacher M, Rücker G (2016),
Modelling multiple thresholds in meta-analysis of diagnostic test accuracy
studies.
\emph{BMC Med Res Methodol}, \bold{16}, 97.
}

\author{
Susanne Steinhauser \email{susanne.steinhauser@uni-koeln.de}, Srinath Kolampally  \email{kolampal@imbi.uni-freiburg.de}, Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}
}


\seealso{
\code{\link{plot.diagmeta}}}
%%, \code{\link{print.diagmeta}},\code{\link{summary.print.diagmeta}},\code{\link{lmer}}}

\examples{
# Fenomeda data
#
data(fenomeda)

str(fenomeda)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
                  studlab = paste(author, year, group),
                  data = fenomeda, log.cutoff = TRUE,
                  model = "DIDS", distr = "l")

plot(diag1)
}
