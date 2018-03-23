#' Meta-analysis of studies of the diagnostic test accuracy of FENO
#' for diagnosis of asthma
#'
#' Meta-analysis of studies of the diagnostic test accuracy of fractional
#' exhaled nitric oxide (FENO) for diagnosis of asthma.
#' 
#' The data were collected for a systematic review by Karrasch et
#' al. (2017) and are published as a supplement (Appendix 1) to
#' Schneider et al. (2017). They have been preprocessed for use in R
#' package \bold{diagmeta}.
#'
#' @docType data
#' 
#' @usage Schneider2017
#'
#' @format A data frame with the following columns:
#' 
#' \itemize{
#'   \item study_id. Numeric study ID
#'   \item author. First author
#'   \item year. Year of publication
#'   \item group. Information on subgroup
#'   \item cutpoint. Cutpoint (FeNO [ppb])
#'   \item tpos. Number of true positives
#'   \item fneg. Number of false negatives
#'   \item fpos. Number of false positives
#'   \item tneg. Number of true negatives
#' }
#' 
#' @source
#'
#'   Karrasch S, Linde K, Rücker G, Sommer H, Karsch-Volk M,
#'   Kleijnen J, Jörres RA, Schneider A (2017), Accuracy of FENO for
#'   diagnosing asthma: a systematic review.  \emph{Thorax},
#'   \bold{72}, 109e16.
#' 
#'   Schneider A, Linde K, Reitsma JB, Steinhauser S, Rücker G (2017),
#'   A novel statistical model for analyzing data of a systematic
#'   review generates optimal cutoff values for fractional exhaled
#'   nitric oxide for asthma diagnosis \emph{Journal of Clinical
#'   Epidemiology}, \bold{92}, 69--78.  doi:
#'   10.1016/j.jclinepi.2017.09.001
#'
#' @examples
#'
#' # FENO dataset
#' #
#' data(Schneider2017)
#' 
#' diag1 <- diagmeta(tpos, fpos, tneg, fneg, cutpoint,
#'                   studlab = paste(author, year, group),
#'                   data = Schneider2017,
#'                   model = "DIDS", log.cutoff = TRUE)
#' 
#' plot(diag1)


"Schneider2017"