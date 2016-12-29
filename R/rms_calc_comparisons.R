#' Calculate estimates (XB or hazard/odds ratios, depending on model type) for specified values of
#' a numeric model covariate, or all levels of a factor covariate, vs a reference level. Store in a
#' data frame.
#'
#' @param rmsObj Model fit object of class rms.
#' @param getRatios Indicator for whether to calculate ratios (exp(XB)) vs estimates on XB scale.
#' Defaults to TRUE if rmsObj is from cph() or lrm().
#' @param df Data frame from which to get variable class information and, if needed, reference
#' and comparison values.
#' @param vname String; name of variable in both names(df) and a covariate in rmsObj.
#' @param refVal Value or factor level to set reference level of df[,vname] to. Default is
#' median (numeric) or first level (factor).
#' @param compVals Numeric; value(s) to set comparison levels of vname to, if df[,vname] is numeric.
#' @param nUnique Integer; if compVals are not specified and number of unique values of variable
#' is < nUnique, will calculate estimates for every unique value vs reference.
#' @param compQuant Numeric vector; if compVals are not specified, gets these quantiles from df to
#' use as default comparison values. Defaults to c(0.1, 0.25, 0.5, 0.75, 0.9).
#' @param rndRC Integer; number of digits to round reference, comparison columns to for numeric
#' variables. Defaults to 2.
#'
#' @return data.frame containing reference, comparison, effect, lower and upper confidence limits,
#' variable name and indicator for whether row contains reference:reference comparison.
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#'
#' @seealso \code{\link[rms]{ols}}, \code{\link[rms]{lrm}}, \code{\link[rms]{cph}},
#' \code{\link[rms]{Gls}}, \code{\link[rms]{summary.rms}}.
#'
#' @examples
#' ## Fit linear regression using ols()
#' mymod <- ols(Sepal.Length ~ Species + Sepal.Width, data = iris)
#'
#' ## Set datadist
#' dd.iris <- datadist(iris)
#' options(datadist = 'dd.iris')
#'
#' ## Continuous covariate
#' rms_calc_comparisons(mymod, vname = Sepal.Width, df = iris)
#'
#' ## Categorical covariate
#' rms_calc_comparisons(mymod, vname = Species, df = iris)
#'

rms_calc_comparisons <- function(rmsObj,
                                 getRatios,
                                 df,
                                 vname,
                                 refVal,
                                 compVals,
                                 nUnique = 10,
                                 compQuant = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                 rndRC = 2){

  ## Argument type tests
  if(!(inherits(rmsObj, 'rms'))){
    stop("rmsObj must be an rms object")
  }
  if(!(vname %in% names(df))){
    stop("vname must be in names(df)")
  }

  ## Do we want ratios (cph or lrm fit) or differences (ols, Gls)?
  if(missing(getRatios)){
    getRatios <- inherits(rmsObj, c('cph', 'lrm'))
  }

  ## For continuous covariates:
  if(is.numeric(df[,vname])){
    ## summary.rms() takes two values for comparison
    summary.text <- "summary(rmsObj, %s = c(refVal, cv), est.all = FALSE, antilog = FALSE)"

    ## Comparison values, if not provided
    if(missing(compVals)){
      ## If variable is numeric but has fewer than nUnique values,
      ##  get estimates for all of them vs reference
      if(length(unique(df[,vname])) < nUnique){
        compVals <- sort(unique(df[,vname]))
      ## Otherwise, use compQuant to get comparison values from data
      } else{
        compVals <- unique(quantile(df[,vname], probs = compQuant, na.rm = TRUE))
      }
    }

    ## Set reference value to median, if not provided
    if(missing(refVal)){
      refVal <- median(df[,vname], na.rm = TRUE)
    }

    ## Calculate estimates on XB scale for all comparison values vs reference
    est.xb <-
      do.call(rbind,
              lapply(compVals,
                     FUN = function(cv){
                       eval(parse(text = sprintf(summary.text, vname)))
                     }))
    est.xb <- as.data.frame(est.xb)

    ## Replace column names to make easier to work with
    colnames(est.xb) <- c('ref', 'comp', 'diff', 'effect', 'se', 'lcl', 'ucl', 'type')

    ## Add variable name
    est.xb$variable <- vname

    ## Make reference, comparison columns character so they can be combined with data frames
    ## created for factors
    est.xb <- est.xb %>%
      mutate(ref.c = as.character(round(ref, rndRC)),
             comp.c = as.character(round(comp, rndRC)))

  ## For factor covariates:
  } else{
    ## summary.rms() takes one level to use as reference
    summary.text <- "summary(rmsObj, %s = refVal, est.all = FALSE, antilog = FALSE)"

    ## If reference value not provided, set to first factor level
    if(missing(refVal)){
      refVal <- levels(df[,vname])[1]
    }

    ## Get normal summary, plus row for reference level
    est.xb <- rbind(c(1, 1, NA, 0, 0, 0, 0, 1),
                    eval(parse(text = sprintf(summary.text, vname))))
    rownames(est.xb)[1] <- paste0(vname, ' - ', refVal, ':', refVal)

    ## Replace column names to make easier to work with
    colnames(est.xb) <- c('ref', 'comp', 'diff', 'effect', 'se', 'lcl', 'ucl', 'type')

    est.xb <- as.data.frame(est.xb)

    ## Replace numeric reference, comparison values with factor level labels originally stored
    ## in rownames
    est.xb$tmp <- rownames(est.xb)
    est.xb <- est.xb %>%
      separate(tmp, into = c('variable', 'comp.c', 'ref.c'), ' - |:')
  }

  estimates <- est.xb %>%
    mutate(is.ref = ref == comp) %>%
    select(-diff, -se, -type)

  if(getRatios){
    estimates <- estimates %>%
      mutate(effect = exp(effect),
             lcl = exp(lcl),
             ucl = exp(ucl))
  }

  return(estimates)

}

