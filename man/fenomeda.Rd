\name{fenomeda}
\alias{fenomeda}

\title{Meta-analysis of studies of the diagnostic test accuracy of FENO for diagnosis of asthma}

\description{
Meta-analysis of studies of the diagnostic test accuracy of fractional
exhaled nitric oxide (FENO) for diagnosis of asthma.

The data were collected for a systematic review by Karrasch et
al. (2017) and are published as a supplement (Appendix 1) to Schneider
et al. (2017). They have been preprocessed for use in R package
diagmeta.
}

\usage{
data(fenomeda)
}


\format{A data frame with the following columns:
  \tabular{ll}{
    \bold{study_id} \tab Numeric study ID \cr
    \bold{author} \tab First author \cr
    \bold{year} \tab Year of publication \cr
    \bold{group} \tab Information on subgroup \cr
    \bold{cutpoint} \tab Cutpoint (FeNO [ppb]) \cr
    \bold{tpos} \tab Number of true positives \cr
    \bold{fneg} \tab Number of false negatives \cr
    \bold{fpos} \tab Number of false positives \cr
    \bold{tneg} \tab Number of true negatives
 }
}

\source{
Karrasch S, Linde K, Rücker G, Sommer H, Karsch-Volk M, Kleijnen J, Jörres RA, Schneider A (2017),
Accuracy of FENO for diagnosing asthma: a systematic review.
\emph{Thorax}, \bold{72}, 109e16.

Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017),
A novel statistical model for analyzing data of a systematic review
generates optimal cutoff values for fractional exhaled nitric oxide
for asthma diagnosis
\emph{Journal of Clinical Epidemiology}, \bold{92}, 69--78.
doi: 10.1016/j.jclinepi.2017.09.001
}

\seealso{\code{\link{diagmeta}}, \code{\link{plot.diagmeta}}}

\examples{
# Feno data
#
data(fenomeda)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
                  studlab = paste(author, year, group),
                  data = fenomeda, log.cutoff = TRUE,
                  model = "DIDS", distr = "l")

plot(diag1)
}

\keyword{datasets}
