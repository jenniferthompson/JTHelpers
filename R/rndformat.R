#' Round and format a number to the same number of digits.
#'
#' @param x Numeric; value to be rounded and formatted
#' @param digits Numeric; number of digits to round and format x to
#'
#' @return Character string of x, rounded and formatted to the same number of digits.
#'
#' @export
#'
#' @examples
#' rndformat(1.9727, digits = 3)
#' rndformat(1.2, digits = 3)
#'

## -- Round and format a number to the same number of digits (default = 2) -------------------------
rndformat <- function(x, digits = 2){
  format(round(x, digits = digits), nsmall = digits)
}

