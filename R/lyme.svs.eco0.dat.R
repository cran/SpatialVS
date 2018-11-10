#' The Lyme disease dataset with Eco id=0
#'
#' The Lyme disease dataset contains case data from 2006 to 2011 demographic
#' data and land cover data in Virginia. Lyme disease case data were collected by Virginia Department of Health.
#' Eco id = 0 represents northern/western subregion, which includes Northern Piedmont, Blue Ridge, Ridge and Valley, and Central Appalachian.
#'
#' @docType data
#'
#' @usage lyme.svs.eco0.dat
#'
#' @format
#' \describe{
#'   \item{\code{y}}{integer vector, output count, each element represents the disease count in one area.}
#'   \item{\code{X}}{Numeric matrix, matrix of covariates, includes percentage of developed land in each census tract, percentage of forest in each census tract, etc.}
#'   \item{\code{offset}}{Numeric vector, vector for offset values, each element represents the population in one area.}
#'   \item{\code{location}}{Numeric matrix, location for each census tract.}
#'   \item{\code{geoid}}{Numeric vecotr, geo id for each census tract.}
#' }
#'
#' @keywords datasets
#'
#'
#'
#'
#'
#' @examples
#' data("lyme.svs.eco0")
"lyme.svs.eco0.dat"
