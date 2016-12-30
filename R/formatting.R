#' Formatting functions: Round and format a number to the same number of digits; format a p-value.
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

#' @param p Numeric; p-value to be rounded and formatted
#'
#' @return Character string; <0.0001, <0.001, or rndformat(p), depending on the value of p.
#'
#' @export
#'
#' @examples
#' formatp(0.000214)
#' formatp(0.28723)
#'

## -- Format p-values to 3 decimal places, <0.001, or <0.0001 --------------------------------------
formatp <- function(p){
  ifelse(p < 0.0001, '<0.0001',
  ifelse(p < 0.001, '<0.001',
         rndformat(p, digits = 3)))
}
