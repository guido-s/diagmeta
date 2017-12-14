\name{plot.diagmeta}
\alias{plot.diagmeta}

\title{
Plot for meta-analysis of diagnostic accuracy studies with the
multiple cutoffs
}

\description{
Draws a plot for meta-analysis of diagnostic accuracy studies with the
multiple cutoffs
}

\usage{
\method{plot}{diagmeta}(x, which = c("survival", "youden", "roc", "sroc"), xlab = "threshold", main, ci = FALSE,
              ciSens = FALSE, ciSpec = FALSE, mark.optcut = FALSE, mark.cutpoints = FALSE,
              points = TRUE, lines = FALSE, rlines = TRUE, line.optcut = TRUE,
              col.points =  "rainbow",cex = 1, pch.points = 16, cex.marks = 0.7 * cex,
              lwd = 1, lwd.optcut = 2 * lwd, shading = "none", xlim = NULL, ...)
}


\arguments{
  \item{x}{An object of class \code{diagmeta}}
  \item{which}{A character vector indicating the type of plot, either
    \code{"reg"} or\code{"cdf"} or \code{"survival"} or
    \code{"Youden"} or \code{"ROC"} or \code{"SROC"} or
    \code{"density"}}
  \item{xlab}{An x axis label}
  \item{main}{A logical indicating title to the plot}
  \item{ci}{A logical indicating whether confidence intervals should
    be plotted for \code{"cdf"}, \code{"survival"} and \code{"Youden"}}
  \item{ciSens}{A logical indicating whether confidence intervals
    should be plotted for sensitivity, given the specificity in
    \code{"SROC"} plot}
  \item{ciSpec}{A logical indicating whether confidence intervals
    should be plotted for specificity, given the sensitivity in
    \code{"SROC"} plot}
  \item{mark.optcut}{A logical indicating whether the optimal cutoff
    should be marked on \code{"SROC"} plot}
  \item{mark.cutpoints}{A logical indicating whether the given cutoffs
    should be marked on \code{"SROC"} plot}
  \item{points}{A logical indicating whether points should be plotted
    on \code{"ROC"} plot}
  \item{lines}{A logical indicating polygonal lines connecting points
    belonging to the same study should be printed in plots
    \code{"reg"}, \code{"cdf"}, \code{"survival"} and \code{"Youden"}} 
  \item{rlines}{A logical indicating whether regression lines or
    curves should be plotted for plots \code{"reg"}, \code{"cdf"},
    \code{"survival"} and \code{"Youden"}}
  \item{line.optcut}{A logical indicating whether a vertical line
    should be plotted at optimal cutoff line}
  \item{col.points}{A character vector indicating color, either
    \code{"rainbow"} or \code{"topo"} or \code{"heat"} or
    \code{"terrain"} or \code{"cm"} or \code{"gray"} or
    \code{"black"}}
  \item{cex}{A numeric indicating magnification to be used for
    plotting text and symbols.}
  \item{pch.points}{A numeric inicating plot symbol(s) for points}
  \item{cex.marks}{A numeric indicating magnification(s) to be used
    for marking cutoffs}
  \item{lwd}{A numeric indicating line width}
  \item{lwd.optcut}{A numeric indicating line width of optimal cutoff}
  \item{shading}{A character indicating shading and hatching
    confidence region in \code{"SROC"} plot, either \code{"none"} or
    \code{"shade"} or \code{"hatch"}}
  \item{xlim}{A character or numerical vector indicating the minimum
    and maximum value for the horizontal axes.}
  \item{\dots}{Additional graphical arguments}
}

\details{
To be written
}

\references{
Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017),
A novel statistical model for analyzing data of a systematic review generates optimal cutoff values for fractional exhaled nitric oxide for asthma diagnosis
\emph{Journal of Clinical Epidemiology}, \bold{92}, 69--78,
doi: 10.1016/j.jclinepi.2017.09.001 .
}


\seealso{
\code{\link{diagmeta}}
}

\examples{
# Fenomeda data
#
data(fenomeda)

diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint, study,
                  data = fenomeda, log.cutoff = TRUE,
                  model = "DIDS", distr = "l")


# Regression plot with confidence intervals
plot(diag1, which="reg", lines = FALSE, ci = TRUE)

# Cumulative distribution plot with optimal cutoff line and confidence intervals
plot(diag1, which="cdf", line.optcut = TRUE, ci= TRUE)

# Survival plot with optimal cutoff line and confidence intervals
plot(diag1, which="survival", line.optcut = TRUE, ci = TRUE)

# Youden plot of optimal cutoff line and confidence intervals
plot(diag1, which="youden", lines = TRUE, line.optcut = TRUE, ci = TRUE)

# ROC plot of lines connecting points belonging to the same study
plot(diag1, which="ROC", lines = TRUE)

# SROC plot of confidence regions for sensitivity and specificity with optimal cutoff mark
plot(diag1, which="SROC", ciSens = TRUE, ciSpec = TRUE, mark.optcut = TRUE, shading = "hatch")
}
