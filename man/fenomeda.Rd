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


\format{A data frame with the following columns:
  \tabular{ll}{
    \bold{study_id} \tab Numeric study ID \cr
    \bold{author} \tab First author \cr
    \bold{year} \tab Year of publication \cr
    \bold{group} \tab Information on subgroup \cr
    \bold{cutpoint} \tab Cutpoint \cr
    \bold{tpos} \tab Number of true positives \cr
    \bold{fneg} \tab Number of false negatives \cr
    \bold{fpos} \tab Number of false positives \cr
    \bold{tneg} \tab Number of true negatives
 }
}

\source{
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
