#' Function to get estimates, on original or ratio scale, for a main variable at a given level of
#' another variable using summary.rms(). Returns a tidy data frame with columns for main and
#' interacting variable names and levels and numeric columns for effects and CIs.
#'
#' This will often be used in conjunction with \code{lapply} or \code{purrr::pmap}. For an example,
#' see chunk \code{death_inhosp_plot} in \href{https://github.com/jenniferthompson/MotoricSubtypes/blob/master/brain_motoric_subtypes.Rmd}{this Rnw file}
#' (function here is called \code{estimate.with.interaction.rms}).
#'
#' @param rmsObjName Character string; name of model fit object of class rms
#' @param estVar Character string; name of variable for which we want estimate/CI
#' @param estVals Vector of length 2; comparison wanted for main variable. Default = datadist() defaults.
#' @param intVar Character string; name of interacting variable
#' @param intAdjust Character or numeric value; value to adjust intVar to. Default = datadist() default.
#' @param getRatios Indicator for whether to calculate ratios (exp(XB)) vs estimates on XB scale.
#' Defaults to TRUE if get(rmsObjName) is from cph() or lrm().
#'
#' @return data.frame containing reference, comparison, effect, lower and upper confidence limits,
#' variable name and indicator for whether row contains reference:reference comparison.
#'
#' @import rms
#' @importFrom broom tidy
#' @importFrom dplyr mutate select filter
#' @importFrom tidyr separate
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @seealso \code{\link[rms]{ols}}, \code{\link[rms]{lrm}}, \code{\link[rms]{cph}},
#' \code{\link[rms]{Gls}}, \code{\link[rms]{summary.rms}}, \code{\link[purrr]{map}},
#' \code{\link[purrr]{pmap}}.
#'
#' @examples
#' ## Fit logistic regression using lrm()
#' mymod <- lrm(Species ~ Sepal.Length * Sepal.Width, data = iris)
#'
#' ## Set datadist
#' dd.iris <- datadist(iris)
#' options(datadist = 'dd.iris')
#'
#' ## Continuous covariate, comparing all quantiles to median by default
#' rms_estimate_with_int('mymod', estVar = 'Sepal.Length', intVar = 'Sepal.Width')
#'

rms_estimate_with_int <- function(rmsObjName,
                                  estVar,
                                  estVals,
                                  intVar,
                                  intAdjust,
                                  getRatios){

  ## If no values for estVals are given, use datadist defaults
  if(missing(estVals)){
    dd <- get(getOption('datadist'))
    ## First choice is Low/High:effect values, used for numeric variables
    estVarLow <- dd$limits['Low:effect', estVar]
    estVarHigh <- dd$limits['High:effect', estVar]
    ## If variable is not numeric, these values are missing; use Low/High instead
    if(is.na(estVarLow)){
      estVarLow <- dd$limits['Low', estVar]
      estVarHigh <- dd$limits['High', estVar]
    }
    ## If est.vals is given, must be a vector of length 2
  } else if(length(estVals) != 2){
    stop('estVals must be a vector of length 2', call. = FALSE)
  } else{
    estVarLow <- estVals[1]
      estVarHigh <- estVals[2]
    }

  ## If variable is a character or factor, add extra quotes for inclusion in sprintf()
  if(!is.numeric(estVarLow)){
    estVarLow <- paste0("'", estVarLow, "'")
    estVarHigh <- paste0("'", estVarHigh, "'")
  }

  ## If no value for getRatios is given, use default based on model type:
  ## Logistic, Cox = yes; all others = no
  if(missing(getRatios)){
    getRatios <- inherits(get(rmsObjName), c('cph', 'lrm'))
  }

  ## If no value is given for intAdjust, use datadist default
  if(missing(intAdjust)){
    intAdjust <- get(getOption('datadist'))$limits['Adjust to', intVar]
  }

  ## If value for interaction term is character, need to add a set of quotes to use in sprintf()
  if(!is.numeric(intAdjust)){
    intAdjustSummary <- paste0("'", intAdjust, "'")
  } else{
    intAdjustSummary <- intAdjust
  }

  ## Use summary.rms() to calculate desired estimate
  summary.text <- sprintf("summary(%s, %s = c(%s, %s), %s = %s, est.all = FALSE, antilog = FALSE)",
                          rmsObjName, estVar, estVarLow, estVarHigh, intVar, intAdjustSummary)
  summary.obj <- eval(parse(text = summary.text))

  ## Use broom to make into a data frame, selecting only row with desired variable, and adding
  ## value that interacting variable is adjusted to
  summary.data <- summary.obj %>%
    broom::tidy() %>%
    tidyr::separate(.rownames,
                    into = c('var.main', 'comp.level', 'ref.level'),
                    sep = " - |:",
                    fill = 'right') %>%
    dplyr::filter(var.main == estVar) %>%
    dplyr::mutate(var.int = intVar,
                  adjust.int = intAdjust,
                  comp.level = ifelse(is.na(comp.level), High, comp.level),
                  ref.level = ifelse(is.na(ref.level), Low, ref.level)) %>%
    ## General cleanup: don't need SE, Type, numeric Low/High columns
    dplyr::select(var.main, ref.level, comp.level, var.int, adjust.int,
                  Effect, Lower.0.95, Upper.0.95)

  ## Change variable names to make easier to work with
  names(summary.data) <- tolower(gsub('\\.0\\.95', 'cl', names(summary.data)))

  ## If we want ratios, exponentiate effect and CLs
  if(getRatios){
    summary.data[,c('effect', 'lowercl', 'uppercl')] <-
      exp(summary.data[,c('effect', 'lowercl', 'uppercl')])
  }

  return(summary.data)

}
