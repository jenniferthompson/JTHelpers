#' Calculate incidence rate ratio and confidence limits from a negative binomial regression
#' model
#'
#' The most common way to express incidence rate ratios is simply exp(beta), representing a one-unit
#' increase in the covariate value; however, in the case of continuous covariates, this is often not
#' a practically meaningful difference (one year of age among adults, eg, or a one-unit change in
#' mean arterial pressure). This function allows you to specify a clinically meaningful comparison
#' in the case of continuous covariates, defaulting to the 75th vs the 25th percentiles of values in
#' the data frame specified.
#'
#' Currently does not handle interaction terms.
#'
#' @param nbObj Model fit of class glm.nb, or mice::mira object fit using glm.nb.
#' @param nbCoefs Vector of coefficients from a glm.nb.
#' @param nbVcov Variance-covariance matrix from a glm.nb fit.
#' @param predVar Character string; name of main predictor variable.
#' @param adjustTo Numeric vector of length 2; values to adjust predVar to if predVar is continuous.
#' Defaults to c(25th, 75th percentiles).
#' @param df Data frame used to determine type of predVar, get defaults for adjustTo and to
#' calculate nonlinear terms if applicable.
#' @param alpha Alpha level for confidence limits. Defaults to 0.05.
#'
#' @return matrix with one row per comparison, containing columns for point estimate, confidence
#' limits, and reference and comparison levels used.
#'
#' @importFrom MASS glm.nb
#' @import mice
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[mice]{mice}}, \code{\link[mice]{pool}};
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
#' @name calc_nb_ratioci

calc_nb_ratioci <- function(nbObj = NULL, ...){ UseMethod("calc_nb_ratioci", nbObj) }

#' @describeIn calc_nb_ratioci Method used when passing coefficients and variance-covariance matrix,
#' vs original model fit
#'
#' Default method uses given coefficients and variance-covariance matrix. Intended use cases include
#' bootstrapped negative binomial models.
#'
#' @importFrom MASS glm.nb
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[Hmisc]{rcspline.eval}} for nonlinear terms.
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
#' mycoefs <- coef(mymod)
#' myvcov <- vcov(mymod)
#'
#' calc_nb_ratioci(nbCoefs = mycoefs, nbVcov = myvcov, predVar = 'v1', df = df)
#'

## Default: use if nbObj is NULL; nbCoefs and nbVcov must be present -
calc_nb_ratioci.default <- function(nbObj = NULL,
                                    nbCoefs,
                                    nbVcov,
                                    predVar,
                                    adjustTo = NULL,
                                    df,
                                    alpha = 0.05){

  if(!(predVar %in% names(df))){
    stop('predVar must be a variable in df', call. = FALSE)
  }

  ## nbObj must be of class negbin or a mira object fit with glm.nb(); if it's either one of these,
  ## other methods will be used, so throw an error if nbObj is supplied
  if(!is.null(nbObj)){
    stop('nbObj must be either a glm.nb() fit, or a mira object fit with glm.nb().', call. = FALSE)
  }

  ## One of nbObj *or* (nbCoefs + nbVcov) must be supplied
  if(missing(nbCoefs) | missing(nbVcov)){
    stop('Must provide either nbObj or both nbCoefs and nbVcov', call. = FALSE)
  }

  ## Make sure coefficient, vcov matrix have equal dimensions
  if(!inherits(nbVcov, 'matrix') |
     nrow(nbVcov) != ncol(nbVcov) |
     nrow(nbVcov) != length(nbCoefs)){
    stop('nbVcov must be a square matrix with nrow, ncol = length(nbCoefs)')
  }

  ## Get model coefficients and variance-covariance matrix
  modcoefs <- nbCoefs
  modvcov <- nbVcov

  ## For continuous covariate
  if(inherits(df[,predVar], "factor")){
    ## Which coefficients reflect categories
    use.coefs <- grep(paste0('^', predVar), names(modcoefs))
    use.betas <- modcoefs[use.coefs]

    ## Since each category only uses one coefficient, only need the variance, not the covariance
    use.vcov <- diag(modvcov)[use.coefs]

    ## beta * xvals is just beta (each X is 1)
    irrvar.logor <- use.betas
    irrvar.or <- exp(irrvar.logor)

    ## Calculate SE
    irrvar.se <- sqrt(use.vcov)

    ## Calculate confidence limits
    critval <- 1 - alpha / 2

    irrvar.lcl <- exp(irrvar.logor - qnorm(critval)*irrvar.se)
    irrvar.ucl <- exp(irrvar.logor + qnorm(critval)*irrvar.se)

    ## Combine point estimate(s), confidence limits, reference (=1 for all) and comparison values
    ## into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl,
                         rep(1, length(irrvar.or)),
                         as.numeric(gsub(predVar, '', names(modcoefs[use.coefs]), fixed = TRUE)))
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  } else{
    ## How many nonlinear terms does var have in the model?
    n.nlterms <- length(unique(grep(paste0(predVar, "'*$"), names(modcoefs), value = TRUE)))

    ## Determine what values to adjust var terms to (by default, 75th and 25th percentiles of df)
    if(is.null(adjustTo)){
      ## Get percentiles for main term
      adjustTo <- as.numeric(quantile(df[,predVar], probs = c(0.25, 0.75), na.rm = TRUE))
    } else if(length(!is.na(as.numeric(adjustTo))) != 2){
      stop('adjustTo must have exactly two non-missing numeric values')
    }

    ## If variable is linear, adjustTo = 2x1 matrix
    if(n.nlterms == 1){
      adjust.matrix <- matrix(adjustTo, ncol = 1)
      ## Otherwise, adjustTo = 2xn.nlterms matrix
    } else{
      use.knots <- Hmisc::rcspline.eval(df[,predVar], nk = n.nlterms + 1, knots.only = TRUE)
      adjust.matrix <- as.matrix(Hmisc::rcspline.eval(adjustTo, knots = use.knots, inclx = TRUE))
    }

    ## Calculate differences in adjustment values - this is what will be multiplied by coefficients
    adjust.diffs <- adjust.matrix[2,] - adjust.matrix[1,]
    names(adjust.diffs) <- paste0(predVar,
                                  unlist(lapply(1:n.nlterms, FUN = function(k){
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

    ## Combine point estimate, confidence limits, reference and comparison values into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl, adjustTo[1], adjustTo[2])
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  }

  return(irr.results)
}

#' @describeIn calc_nb_ratioci Method for glm.nb models fit with mice objects.
#' @importFrom MASS glm.nb
#' @import mice
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[mice]{mice}} \code{\link[Hmisc]{rcspline.eval}}
#' for nonlinear terms.
#'
#' @examples
#'
#' ## Create data frame
#' df <- data.frame(y = round(rexp(n = 100, rate = 0.5)),
#'                  v1 = sample(c(1:5, NA), size = 100, replace = TRUE),
#'                  v2 = rnorm(n = 100))
#'
#' ## Impute missing data using mice
#' mids.df <- mice(df)
#'
#' ## Fit negative binomial model
#' mymod <- with(mids.df, MASS::glm.nb(y ~ v1 * v2))
#'
#' calc_nb_ratioci(mymod, predVar = 'v1', df = df)
#'

calc_nb_ratioci.mira <- function(nbObj,
                                 predVar,
                                 adjustTo = NULL,
                                 df,
                                 alpha = 0.05){

  if(!(predVar %in% names(df))){
    stop('predVar must be a variable in df', call. = FALSE)
  }

  ## Make sure model was fit using glm.nb
  if(!(inherits(nbObj$analyses[[1]], 'negbin'))){
    stop('nbObj must be of class negbin from glm.nb(), or a mice object using glm.nb to fit',
         call. = FALSE)
  }

  ## Get model coefficients and variance-covariance matrix
  modcoefs <- pool(nbObj)$qbar
  modvcov <- pool(nbObj)$t

  ## For continuous covariate
  if(inherits(df[,predVar], "factor")){
    ## Which coefficients reflect categories
    use.coefs <- grep(paste0('^', predVar), names(modcoefs))
    use.betas <- modcoefs[use.coefs]

    ## Since each category only uses one coefficient, only need the variance, not the covariance
    use.vcov <- diag(modvcov)[use.coefs]

    ## beta * xvals is just beta (each X is 1)
    irrvar.logor <- use.betas
    irrvar.or <- exp(irrvar.logor)

    ## Calculate SE
    irrvar.se <- sqrt(use.vcov)

    ## Calculate confidence limits
    critval <- 1 - alpha / 2

    irrvar.lcl <- exp(irrvar.logor - qnorm(critval)*irrvar.se)
    irrvar.ucl <- exp(irrvar.logor + qnorm(critval)*irrvar.se)

    ## Combine point estimate(s), confidence limits, reference (=1 for all) and comparison values
    ## into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl,
                         rep(1, length(irrvar.or)),
                         as.numeric(gsub(predVar, '', names(modcoefs[use.coefs]), fixed = TRUE)))
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  } else{
    ## How many nonlinear terms does var have in the model?
    n.nlterms <- length(unique(grep(paste0(predVar, "'*$"), names(modcoefs), value = TRUE)))

    ## Determine what values to adjust var terms to (by default, 75th and 25th percentiles of df)
    if(is.null(adjustTo)){
      ## Get percentiles for main term
      adjustTo <- as.numeric(quantile(df[,predVar], probs = c(0.25, 0.75), na.rm = TRUE))
    } else if(length(!is.na(as.numeric(adjustTo))) != 2){
      stop('adjustTo must have exactly two non-missing numeric values')
    }

    ## If variable is linear, adjustTo = 2x1 matrix
    if(n.nlterms == 1){
      adjust.matrix <- matrix(adjustTo, ncol = 1)
      ## Otherwise, adjustTo = 2xn.nlterms matrix
    } else{
      use.knots <- Hmisc::rcspline.eval(df[,predVar], nk = n.nlterms + 1, knots.only = TRUE)
      adjust.matrix <- as.matrix(Hmisc::rcspline.eval(adjustTo, knots = use.knots, inclx = TRUE))
    }

    ## Calculate differences in adjustment values - this is what will be multiplied by coefficients
    adjust.diffs <- adjust.matrix[2,] - adjust.matrix[1,]
    names(adjust.diffs) <- paste0(predVar,
                                  unlist(lapply(1:n.nlterms, FUN = function(k){
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

    ## Combine point estimate, confidence limits, reference and comparison values into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl, adjustTo[1], adjustTo[2])
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  }

  return(irr.results)
}

#' @describeIn calc_nb_ratioci Method for glm.nb models (no imputation).
#' @importFrom MASS glm.nb
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[Hmisc]{rcspline.eval}} for nonlinear terms.
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

calc_nb_ratioci.negbin <- function(nbObj,
                                   predVar,
                                   adjustTo = NULL,
                                   df,
                                   alpha = 0.05){

  if(!(predVar %in% names(df))){
    stop('predVar must be a variable in df', call. = FALSE)
  }

  ## Get model coefficients and variance-covariance matrix
  modcoefs <- coef(nbObj)
  modvcov <- vcov(nbObj)

  ## For continuous covariate
  if(inherits(df[,predVar], "factor")){
    ## Which coefficients reflect categories
    use.coefs <- grep(paste0('^', predVar), names(modcoefs))
    use.betas <- modcoefs[use.coefs]

    ## Since each category only uses one coefficient, only need the variance, not the covariance
    use.vcov <- diag(modvcov)[use.coefs]

    ## beta * xvals is just beta (each X is 1)
    irrvar.logor <- use.betas
    irrvar.or <- exp(irrvar.logor)

    ## Calculate SE
    irrvar.se <- sqrt(use.vcov)

    ## Calculate confidence limits
    critval <- 1 - alpha / 2

    irrvar.lcl <- exp(irrvar.logor - qnorm(critval)*irrvar.se)
    irrvar.ucl <- exp(irrvar.logor + qnorm(critval)*irrvar.se)

    ## Combine point estimate(s), confidence limits, reference (=1 for all) and comparison values
    ## into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl,
                         rep(1, length(irrvar.or)),
                         as.numeric(gsub(predVar, '', names(modcoefs[use.coefs]), fixed = TRUE)))
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  } else{
    ## How many nonlinear terms does var have in the model?
    n.nlterms <- length(unique(grep(paste0(predVar, "'*$"), names(modcoefs), value = TRUE)))

    ## Determine what values to adjust var terms to (by default, 75th and 25th percentiles of df)
    if(is.null(adjustTo)){
      ## Get percentiles for main term
      adjustTo <- as.numeric(quantile(df[,predVar], probs = c(0.25, 0.75), na.rm = TRUE))
    } else if(length(!is.na(as.numeric(adjustTo))) != 2){
      stop('adjustTo must have exactly two non-missing numeric values')
    }

    ## If variable is linear, adjustTo = 2x1 matrix
    if(n.nlterms == 1){
      adjust.matrix <- matrix(adjustTo, ncol = 1)
      ## Otherwise, adjustTo = 2xn.nlterms matrix
    } else{
      use.knots <- Hmisc::rcspline.eval(df[,predVar], nk = n.nlterms + 1, knots.only = TRUE)
      adjust.matrix <- as.matrix(Hmisc::rcspline.eval(adjustTo, knots = use.knots, inclx = TRUE))
    }

    ## Calculate differences in adjustment values - this is what will be multiplied by coefficients
    adjust.diffs <- adjust.matrix[2,] - adjust.matrix[1,]
    names(adjust.diffs) <- paste0(predVar,
                                  unlist(lapply(1:n.nlterms, FUN = function(k){
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

    ## Combine point estimate, confidence limits, reference and comparison values into matrix
    irr.results <- cbind(irrvar.or, irrvar.lcl, irrvar.ucl, adjustTo[1], adjustTo[2])
    rownames(irr.results) <- NULL
    colnames(irr.results) <- c('pointest', 'lcl', 'ucl', 'ref.val', 'comp.val')

  }

  return(irr.results)
}
