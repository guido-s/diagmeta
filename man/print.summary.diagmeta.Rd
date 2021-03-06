% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.summary.diagmeta.R
\name{print.summary.diagmeta}
\alias{print.summary.diagmeta}
\title{Print method for summary of diagmeta objects}
\usage{
\method{print}{summary.diagmeta}(x, digits = 3, ...)
}
\arguments{
\item{x}{An object of class \code{summary.diagmeta}.}

\item{digits}{Number of significant digits for printing.}

\item{\dots}{Additional arguments.}
}
\description{
Print method for objects of class \code{summary.diagmeta}.
}
\examples{
# FENO dataset
#
data(Schneider2017)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
                  studlab = paste(author, year, group),
                  data = Schneider2017,
                  log.cutoff = TRUE)

summary(diag1)
print(summary(diag1), digits = 2)

}
\seealso{
\code{\link{diagmeta}} \code{\link{summary.diagmeta}}
}
\author{
Gerta Rücker \email{ruecker@imbi.uni-freiburg.de},
Susanne Steinhauser \email{susanne.steinhauser@uni-koeln.de},
Srinath Kolampally \email{kolampal@imbi.uni-freiburg.de},
Guido Schwarzer \email{sc@imbi.uni-freiburg.de}
}
\keyword{print}
