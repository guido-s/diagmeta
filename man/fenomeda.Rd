\name{fenomeda}
\alias{fenomeda}

\title{Meta-analysis of diagnostic accuracy studies with the multiple cutoffs model in asthma}

\description{
Meta-analysis of diagnostic accuracy studies with the optimal cutoff
values for fractional exhaled nitric oxide in asthma.

These data are used as an example in Schneider et al. (2017) and have
been preprocessed for use in R package diagmeta.
}

\usage{
data(fenomeda)
}


\format{
  \tabular{rll}{
    \tab \bold{study}    \tab First author and year \cr
    \tab \bold{study_id} \tab Numeric study ID \cr
    \tab \bold{n}        \tab Number of participants \cr
    \tab \bold{prev}     \tab Prevalence \cr
    \tab \bold{cutpoint} \tab Cutpoint \cr
    \tab \bold{tpos}     \tab True positives\cr
    \tab \bold{fneg}     \tab False Negatives (degrees F)\cr
    \tab \bold{fpos}     \tab False positives (lang)\cr
    \tab \bold{tneg}     \tab True negatives (mph)\cr
  }
}

\references{
Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017),
A novel statistical model for analyzing data of a systematic review generates optimal cutoff values for fractional exhaled nitric oxide for asthma diagnosis
\emph{Journal of Clinical Epidemiology}, \bold{92}, 69--78,
doi: 10.1016/j.jclinepi.2017.09.001 .
}

\author{
Srinath Kolampally \url{kolampal@imbi.uni-freiburg.de}, Gerta Rücker \url{ruecker@imbi.uni-freiburg.de}
}

\seealso{\code{\link{diagmeta}}, \code{\link{plot.diagmeta}}}

\examples{
# Feno data
#
data(fenomeda)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint, study,
                  data = fenomeda, log.cutoff = TRUE,
                  model = "DIDS", distr = "l")

plot(diag1)
}

\keyword{datasets}
