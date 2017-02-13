#' Likelihood ratio p-values for full model vs model with one term (or consecutive terms) removed
#'
#' Runs model without `removeTerm` and performs likelihood ratio test for full model vs nested model.
#' `function(variable)` - eg, `rcs(age, 3)` - is considered one term, as are all levels of a
#' categorical variable.
#'
#' @param orgModel Model object containing full model.
#' @param removeTerm Character string of term(s) to be removed.
#' @param ... Additional arguments specific to method.
#'
#' @importFrom lmtest lrtest
#'
#' @return List with following values from lrtest() or pool.compare(): test degrees of freedom,
#' X^2^ value, and p-value.
#'
#' @export
#'
#' @seealso \code{lrtest_removeTerm.default, lrtest_removeTerm.mira}.

lrtest_removeTerm <- function(orgModel, removeTerm, ...) UseMethod("lrtest_removeTerm", orgModel)

#' @describeIn lrtest_removeTerm Default method.
#' @param orgModel Model object containing full model.
#' @param removeTerm Character string of term(s) to be removed.
#'
#' @importFrom lmtest lrtest
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
#' lrtest_removeTerm(fullModel, 'Petal.Length')
#'
#' ## What about the effect of both Petal variables?
#' lrtest_removeTerm(fullModel, 'Petal.Length + Petal.Width')
#'
lrtest_removeTerm.default <- function(orgModel, removeTerm){
  modcall <- as.character(orgModel$call)

  if(length(grep(removeTerm, modcall[2], fixed = TRUE)) == 0){
    stop('removeTerm must be on the righthand side of the model', call. = FALSE)
  }

  newform <- sprintf('%s(%s, data = %s)',
                     modcall[1],
                     gsub(' \\+ $', '',
                          gsub('(?<=[_+]) +\\+', '',
                               gsub(removeTerm, '', modcall[2], fixed = TRUE),
                                    perl = TRUE),
                               perl = TRUE),
                          modcall[3])

  newModel <- eval(parse(text = newform))

  lrtestObj <- lmtest::lrtest(orgModel, newModel)

  return(list("df" = lrtestObj[1, 1] - lrtestObj[1, 2],
              "chisq" = lrtestObj[2, 'Chisq'],
              "pvalue" = lrtestObj[2, 'Pr(>Chisq)']))
}

#' @describeIn lrtest_removeTerm Method for models fit with the mice package.
#'
#' @param orgModel Model object containing full model.
#' @param removeTerm Character string of term(s) to be removed.
#' @param miceObjName Character string; name of `mice` object.
#'
#' @import mice
#'
#' @export
#'
#' @seealso \code{\link[mice]{pool.compare}}.
#'
#' @examples
#'
#' ## Add missingness to iris data
#' iris_missing <- iris
#' iris_missing[sample(1:nrow(iris), size = round(nrow(iris_missing) / 10)), 'Sepal.Width'] <- NA
#'
#' ## Create mice object
#' iris_mice <- mice(iris_missing)
#'
#' #' ## Fit full model using iris data
#' fullModel <- with(iris_mice, lm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width))
#'
#' ## What about just the effect of Petal.Length?
#' lrtest_removeTerm(fullModel, 'Petal.Length', 'iris_mice')
#'
#' ## What about the effect of both Petal variables?
#' lrtest_removeTerm(fullModel, 'Petal.Length + Petal.Width', 'iris_mice')

lrtest_removeTerm.mira <- function(orgModel, removeTerm, miceObjName){
  ## Extract original call from orgModel
  modcall <- as.character(orgModel$analyses[[1]]$call)

  ## Test: is removeTerm included in the original model?
  if(length(grep(removeTerm, modcall[2], fixed = TRUE)) == 0){
    stop('removeTerm must be on the righthand side of the model', call. = FALSE)
  }

  ## Take removeTerm out of model RHS and recreate formula argument
  newmodform <- gsub(' \\+ $', '',
                     gsub('(?<=[_+]) +\\+', '',
                          gsub(removeTerm, '', modcall[2], fixed = TRUE),
                          perl = TRUE),
                     perl = TRUE)

  ## Recreate entire model call
  newmodcall <- sprintf('with(%s, %s(%s))', miceObjName, modcall[1], newmodform)

  ## pool.compare() requires term(s) that are different to be at the end of model formula;
  ## may need to refit original model.

  ## What is the original model fit with removeTerm at the end?
  orgmodform <- paste(newmodform, removeTerm, sep = ' + ')
  ## Is the original model's RHS already in this order? If not, refit.
  if(!(orgmodform == modcall[2])){
    orgmodcall <- sprintf('with(%s, %s(%s))', miceObjName, modcall[1], orgmodform)
    orgModel <- eval(parse(text = orgmodcall))
  }

  ## Fit model without removeTerm
  newModel <- eval(parse(text = newmodcall))

  lrtestObj <- pool.compare(orgModel, newModel, data = get(miceObjName), method = 'likelihood')

  return(list("df" = lrtestObj$df1,
              "chisq" = lrtestObj$Dm,
              "pvalue" = lrtestObj$pvalue))
}
