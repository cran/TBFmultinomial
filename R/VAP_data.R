#' Data on VAP acquistion in one ICU
#'
#' It is a simulated dataset with information very to the inform OUTCOMEREA database whose only purpose will be to
#' test an illustrate the functions of this package.
#'
#' @format A data frame with 1640 rows and 7 variables on 90 distinct patients:
#'  \describe{
#'   \item{ID}{distinct ID for each patient}
#'   \item{day}{day of ventilation, day = 1 is the first day of ventilation}
#'   \item{type}{is it a medical or a surgical patient}
#'   \item{gender}{gender of the patient, 1 = male, 0 = female}
#'   \item{SAPSadmission}{the SAPS 2 score at admission to the ICU}
#'   \item{SOFA}{the daily SOFA score}
#'   \item{outcome}{final outcome after the first observation period}
#' }
#' @docType data
#' @name VAP_data
#' @usage data(VAP_data)
#' @format A data frame with 1640 rows and 7 variables
NULL
