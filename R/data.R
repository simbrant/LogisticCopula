#' Example data set
#'
#' This radar data was collected by a system in Goose Bay, Labrador. 
#' This system consists of a phased array of 16 high-frequency antennas with a 
#' total transmitted power on the order of 6.4 kilowatts. 
#' See Sigillito, V. G., Wing, S. P., Hutton, L. V., & Baker, K. B. (1989) 
#' for more details. The targets were free electrons in the ionosphere. "Good"
#' radar returns are those showing evidence of some type of structure in the
#' ionosphere.  "Bad" returns are those that do not; their signals pass through
#' the ionosphere.

#'
#' @name Ionosphere
#' @docType data
#' @usage data(Ionosphere)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{351 by 34 matrix of numeric values.}
#'   \item{Class}{Character vector of length 351 containing 126 entries labeled "bad" and 225 labeled "good".}
#' }
"Ionosphere"
