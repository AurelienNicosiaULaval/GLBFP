#' River Ashuapmushuan daily flow and level data
#'
#' Daily observations of river flow and level for the Ashuapmushuan river.
#'
#' @format A data frame with 4,389 rows and 3 variables:
#' \describe{
#'   \item{flow}{Flow rate in cubic meters per second.}
#'   \item{level}{Water level in meters.}
#'   \item{day}{Day code as integer in `YYYYDDD` format.}
#' }
#'
#' @details
#' Data cover 22 March 1992 to 30 September 2007 with a small number of missing
#' calendar days.
#'
#' @source Environment and Climate Change Canada, Historical Hydrometric Data.
#'   The exact extraction query still needs to be documented in `data-raw/`.
#'
#' @examples
#' data(ashua)
#' summary(ashua)
"ashua"
