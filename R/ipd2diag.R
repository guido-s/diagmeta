#' Individual participant data to enter them into diagmeta
#' 
#' Function to transform individual patient data (IPD) to enter them
#' into \code{diagmeta}
#' 
#' @param studlab A vector with study labels
#' @param value A vector with individual patients' measurements of a
#'   discrete or continuous variable
#' @param status A vector with information of the individual's status
#'   (0 = non-diseased, 1 = diseased)
#' 
#' @return A data frame with values that can be entered into
#'   \code{\link{diagmeta}}.
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Srinath
#'   Kolampally \email{kolampal@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}, \link{plot.diagmeta},
#'   \link{print.diagmeta}, \link{summary.diagmeta}}
#' 
#' @examples
#' 
#' # Simulate IPD data for three studies, each with 30 patients based
#' # on normally distributed marker values
#' #
#' k <- 3
#' n <- 30
#' m <- c(20, 23, 26)
#' d <- 10
#' s <- 5
#' studlab <- c(rep(1, n), rep(2, n), rep(3, n))
#' status <- rep(c(rep(0, n / 2), rep(1, n / 2)), k)
#' measurement <- c(rnorm(n / 2, m[1], s), rnorm(n/2, m[1] + d, s),
#'                  rnorm(n / 2, m[2], s), rnorm(n/2, m[2] + d, s), 
#'                  rnorm(n / 2, m[3], s), rnorm(n/2, m[3] + d, s))
#' #
#' IPDdata <- data.frame(studlab, measurement, status)
#' str(IPDdata)
#' 
#' # Transform these data using ipd2diag()
#' #
#' diagdata <- ipd2diag(studlab, value = measurement, status = status)
#' str(diagdata)
#' 
#' # Run diagmeta()
#' #
#' diag1 <- diagmeta(TP, FP, TN, FN, cutoff, studlab,
#'                   data = diagdata, 
#'                   model = "DIDS", distr = "normal")
#' summary(diag1)
#' plot(diag1)
#' par(mfrow = c(1, 2))
#' plot(diag1, which = "ROC", lines = TRUE)
#' plot(diag1, which = "SROC", ciSens = TRUE,
#'      ciSpec = TRUE, lines = TRUE, shading = "hatch")
#'
#' @export


ipd2diag <- function(studlab, value, status) {
  
  cutoffs <- unique(value)
  cutoffs <- cutoffs[order(cutoffs)]
  ##
  diagdata <- data.frame(studlab = NA, cutoff = NA,
                         TP = NA, FP = NA, FN = NA, TN = NA)[-1,]
  ##
  for (s in unique(studlab)) {
    for (c in cutoffs) {
      TP <- sum(studlab == s & status == 1 & value >= c)
      FP <- sum(studlab == s & status == 0 & value >= c)
      FN <- sum(studlab == s & status == 1 & value < c)
      TN <- sum(studlab == s & status == 0 & value < c)
      diagdata <- rbind(diagdata,
                        data.frame(studlab = s, cutoff = c,
                                   TP = TP, FP = FP, FN = FN, TN = TN))
    }
  }
  ##
  diagdata
}
