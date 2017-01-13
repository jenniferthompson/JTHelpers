#' Format a p-value, somewhat per NEJM style.
#'
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
