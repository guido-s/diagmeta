#' Individual participant data to enter them into diagmeta
#' 
#' @description
#' Function to transform individual patient data (IPD) to enter them
#' into \code{diagmeta}
#' 
#' @param studlab A vector with study labels
#' @param value A vector with individual patients' measurements of a
#'   discrete or continuous variable
#' @param status A vector with information of the individual's status
#'   (0 = non-diseased, 1 = diseased)
#' @param data An optional data frame containing the study information
#' @param direction A character string specifying whether the probability of
#'   the target condition (e.g., a disease) is \code{"increasing"} or
#'   \code{"decreasing"} with higher values of the biomarker, can be
#'   abbreviated (see \code{\link{diagmeta}}).
#' 
#' @return
#' A data frame with additional class 'ipd2diag' containing the following
#' variables:
#' \item{studlab}{As defined above.}
#' \item{cutoff}{Cutoff values.}
#' \item{TP, FP, TN, FN}{Number of true positives, false positives,
#'   true negatives and false negatives.}
#' 
#' @author
#' Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de},
#' Srinath Kolampally \email{kolampal@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{diagmeta}, \link{plot.diagmeta},
#'   \link{print.diagmeta}, \link{summary.diagmeta}}
#' 
#' @examples
#' # Simulate IPD data for three studies, each with 30 patients based
#' # on normally distributed marker values
#' #
#' set.seed(20)
#' k <- 3
#' n <- 60
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
#'                   distr = "normal")
#' summary(diag1)
#' plot(diag1)
#' par(mfrow = c(1, 2))
#' plot(diag1, which = "ROC", lines = TRUE)
#' plot(diag1, which = "SROC", ciSens = TRUE,
#'      ciSpec = TRUE, lines = TRUE, shading = "hatch")
#' 
#' @export


ipd2diag <- function(studlab, value, status, data = NULL,
                     direction = "increasing") {
  
  #
  #
  # (1) Read data
  #
  #
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (nulldata)
    data <- sfsp
  #
  # Catch 'studlab', 'value', and 'status'
  #
  studlab <- catch("studlab", mc, data, sfsp)
  chknull(studlab)
  k.All <- length(studlab)
  #
  value <- catch("value", mc, data, sfsp)
  chknull(value)
  chknumeric(value)
  chklength(value, k.All, "value", name = "studlab")
  #
  status <- catch("status", mc, data, sfsp)
  chknull(status)
  chklength(status, k.All, "status", name = "studlab")
  #
  if (is.logical(status))
    status <- as.numeric(status)
  #
  direction <- setchar(direction, c("increasing", "decreasing"))
  
  
  #
  #
  # (2) Generate data set as input to diagmeta()
  #
  #
  
  if (direction == "decreasing")
    value <- -value
  #
  cutoffs <- unique(value)
  #
  if (direction == "decreasing")
    cutoffs <- cutoffs[rev(order(cutoffs))]
  else
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
      #
      diagdata <- rbind(diagdata,
                        data.frame(studlab = s, cutoff = c,
                                   TP = TP, FP = FP, FN = FN, TN = TN))
    }
  }
  #
  if (direction == "decreasing")
    diagdata$cutoff <- -diagdata$cutoff
  #
  attr(diagdata, "direction") <- direction
  class(diagdata) <- c("ipd2diag", class(diagdata))
  #
  diagdata
}
