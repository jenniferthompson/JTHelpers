#' Likelihood ratio p-values for full model vs model with one term removed
#' Runs model without `removeTerm` and performs likelihood ratio test for full model vs nested model.
#' `function(variable)` - eg, `rcs(age, 3)` - is considered one term, as are all levels of a
#' categorical variable.
#'
#' @param orgModel Model object containing full model.
#' @param removeTerm Character string of term(s) to be removed.
#'
#' @importFrom lmtest lrtest
#'
#' @return Numeric value; p-value from `lrtest()`.
#'
#' @export
#'
#' @seealso \code{\link[lmtest]{lrtest}}.
#'
#' @examples
#'
#' ## Fit full model using iris data
#' fullModel <- lm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris)
#'
#' ## What about just the effect of Petal.Length?
#' model_lrpval(fullModel, 'Petal.Length')
#'
#' ## What about the effect of both Petal variables?
#' model_lrpval(fullModel, 'Petal.Length + Petal.Width')
#'

model_lrpval <- function(orgModel, removeTerm){
  modcall <- as.character(orgModel$call)

  newform <- sprintf('%s(%s, data = %s)',
                     modcall[1],
                     gsub(' \\+ $', '',
                          gsub('(?<=[_+]) +\\+', '',
                               gsub(removeTerm, '', modcall[2], fixed = TRUE),
                                    perl = TRUE),
                               perl = TRUE),
                          modcall[3])

  newModel <- eval(parse(text = newform))
  return(as.data.frame(lmtest::lrtest(orgModel, newModel))[2, 'Pr(>Chisq)'])
}
