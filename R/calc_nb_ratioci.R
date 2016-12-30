#' Calculate incidence rate ratio and 95% confidence limits from a negative binomial regression
#' model for a given comparison in a continuous predictor variable.
#'
#' Currently doesn't handle interaction terms or categorical predVars.
#'
#' @param nbObj Model fit of class glm.nb, or mice::mira object fit using glm.nb.
#' @param nbCoefs Vector of coefficients for a glm.nb model.
#' @param nbVcov Variance-covariance matrix for a glm.nb model.
#' @param predVar Character string; name of main predictor variable.
#' @param adjustTo Numeric vector of length 2; values to adjust predVar to. Defaults to c(25th,
#' 75th percentiles).
#' @param df Data frame used to get defaults for adjustTo and to calculate nonlinear terms if
#' applicable.
#' @param alpha Alpha level for confidence limits. Defaults to 0.05.
#'
#' @return data.frame with one row, containing columns for point estimate, confidence limits, and
#' reference and comparison values used.
#'
#' @importFrom MASS glm.nb
#' @import mice
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[mice]{is.mice}}, \code{\link[mice]{pool}};
#' \code{\link[Hmisc]{rcspline.eval}} for nonlinear terms.
#'
#' @examples
#'
#' ## Create data frame
#' df <- data.frame(y = round(rexp(n = 100, rate = 0.5)),
#'                  v1 = sample(1:5, size = 100, replace = TRUE),
#'                  v2 = rnorm(n = 100))
#'
#' ## Fit negative binomial model
#' mymod <- MASS::glm.nb(y ~ v1 * v2, data = df)
#'
#' calc_nb_ratioci(mymod, predVar = 'v1', df = df)
#'

calc_nb_ratioci <- function(nbObj,
                            nbCoefs,
                            nbVcov,
                            predVar,
                            adjustTo = NULL,
                            df,
                            alpha = 0.05){

  if(!(predVar %in% names(df))){
    stop('predVar must be a variable in df', call. = FALSE)
  }

  ## One of nbObj *or* (nbCoefs + nbVcov) must be supplied
  if(missing(nbObj) & (missing(nbCoefs) & missing(nbVcov))){
    stop('Must provide either nbObj or both nbCoefs and nbVcov', call. = FALSE)
  }

  ## If nbObj is given, make sure it is a glm.nb fit (or result of mice + glm.nb)
  if(!missing(nbObj)){
    is.mice <- inherits(nbObj, 'mira')

    ## Model object must be fit using glm.nb()
    if(!(inherits(nbObj, 'negbin') | (is.mice & inherits(nbObj$analyses[[1]], 'negbin')))){
      stop('nbObj must be of class negbin from glm.nb(), or a mice object using glm.nb to fit',
           call. = FALSE)
    }

    ## Get model coefficients, vcov matrix
    if(is.mice){
      modcoefs <- pool(nbObj)$qbar
      modvcov <- pool(nbObj)$t
    } else{
      modcoefs <- coef(nbObj)
      modvcov <- vcov(nbObj)
    }
  } else{
    ## Make sure coefficient, vcov matrix have equal dimensions
    if(!inherits(nbVcov, 'matrix') |
       nrow(nbVcov) != ncol(nbVcov) |
       nrow(nbVcov) != length(nbCoefs)){
      stop('nbVcov must be a square matrix with nrow, ncol = length(nbCoefs)')
    }

    modcoefs <- nbCoefs
    modvcov <- nbVcov
  }

  ## How many nonlinear terms does var have in the model?
  n.terms <- length(unique(grep(paste0(predVar, "'*$"), names(modcoefs), value = TRUE)))

  ## Determine what values to adjust var terms to (by default, 75th and 25th percentiles of df)
  if(is.null(adjustTo)){
    ## Get percentiles for main term
    adjustTo <- as.numeric(quantile(df[,predVar], probs = c(0.25, 0.75), na.rm = TRUE))
  } else if(length(!is.na(as.numeric(adjustTo))) != 2){
    stop('adjustTo must have exactly two non-missing numeric values')
  }

  ## If variable is linear, adjustTo = 2x1 matrix
  if(n.terms == 1){
    adjust.matrix <- matrix(adjustTo, ncol = 1)
    ## Otherwise, adjustTo = 2xn.terms matrix
  } else{
    use.knots <- Hmisc::rcspline.eval(df[,predVar], nk = n.terms + 1, knots.only = TRUE)
    adjust.matrix <- as.matrix(Hmisc::rcspline.eval(adjustTo, knots = use.knots, inclx = TRUE))
  }

  ## Calculate differences in adjustment values - this is what will be multiplied by coefficients
  adjust.diffs <- adjust.matrix[2,] - adjust.matrix[1,]
  names(adjust.diffs) <- paste0(predVar,
                                unlist(lapply(1:n.terms, FUN = function(k){
                                  paste(rep("'", k - 1), collapse = '') })))

  ## Which coefficients do we need to use?
  use.coefs <- grep(paste0(predVar, "'*$"), names(modcoefs))
  use.betas <- modcoefs[use.coefs]

  ## Get vcov matrix for involved components
  use.vcov <- modvcov[use.coefs, use.coefs]

  ## Calculate each component of linear predictor: beta * xvals
  beta.x <- unlist(lapply(1:length(use.betas), FUN = function(b){
    prod(c(use.betas[b], adjust.diffs[b]))
  }))

  irrvar.logor <- sum(beta.x)
  irrvar.or <- exp(irrvar.logor)

  ## Calculate SE
  irrvar.se <- sqrt(adjust.diffs %*% use.vcov %*% adjust.diffs)

  critval <- 1 - alpha / 2

  irrvar.lcl <- exp(irrvar.logor - qnorm(critval)*irrvar.se)
  irrvar.ucl <- exp(irrvar.logor + qnorm(critval)*irrvar.se)

  return(c('pointest' = irrvar.or, 'lcl' = irrvar.lcl, 'ucl' = irrvar.ucl,
           'ref.val' = adjustTo[1], 'comp.val' = adjustTo[2]))
}