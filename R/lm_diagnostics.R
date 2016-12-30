#' Check assumptions of heteroscedasticity and normally distributed errors for a linear regression
#' model fit.
#'
#' @param lmObj Model fit of class `lm`.
#' @param outcomeString String to identify outcome in X axis label. Defaults to "Outcome."
#' @param titleString String to add to plot titles.
#'
#' @export
#'
#' @examples
#'
#' ## Linear regression model
#' mymod <- lm(Sepal.Width ~ Sepal.Length, data = iris)
#'
#' lm_diagnostics(mymod, outcomeString = 'Sepal.Width', titleString = 'Sepal.Width vs Sepal.Length')
#'

lm_diagnostics <- function(lmObj,
                           outcomeString = 'Outcome',
                           titleString){
  if(!(inherits(lmObj, 'lm'))){
    stop('lmObj must be a model fit of type lm', call. = FALSE)
  }

  if(missing(titleString)){
    qq.title <- 'Q-Q of residuals'
    rp.title <- 'RP plot'
  } else{
    qq.title <- paste('Q-Q of residuals,', titleString)
    rp.title <- paste('RP plot,', titleString)
  }

  ## fit.mult.impute objects work differently than non-imputed objects; use residuals and fitted
  ##  values from model object directly rather than resid() and fitted()
  if(inherits(lmObj, 'fit.mult.impute')){
    plot(lmObj$residuals ~ lmObj$fitted.values,
         xlab = paste('Predicted', outcomeString),
         ylab = paste('Model residual'),
         main = rp.title,
         col = 'turquoise4')
    abline(h = 0)
    qqnorm(lmObj$residuals, datax = TRUE, main = qq.title)
  } else{
    plot(resid(lmObj) ~ fitted(lmObj),
         xlab = paste('Predicted', outcomeString),
         ylab = paste('Model residual'),
         main = rp.title,
         col = 'turquoise4')
    abline(h = 0)
    qqnorm(resid(lmObj), datax = TRUE, main = qq.title)
  }
}
