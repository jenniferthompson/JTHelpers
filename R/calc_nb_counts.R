#' Calculate predicted value and standard error of a negative binomial regression model for each
#' row in a design matrix.
#'
#' @param nbObj Model fit of class glm.nb, or mice::mira object fit using glm.nb.
#' @param designMatrix Matrix of covariate values. Number of columns = number of coefficients in
#' nbObj.
#' @param predVar Character string; name of main predictor variable. Defaults to NULL, in which case
#' predictor variable values will not be included in resulting data.frame.
#' @param intVar Character string (optional); name of interacting variable. Defaults to NULL. If
#' included, nothing changes except this column in designMatrix will also be included in the
#' returned data.frame.
#'
#' @return data.frame containing columns for main predictor variable; interacting variable (if
#' given); linear predictor and its SE; adjusted count and its lower and upper confidence limits.
#'
#' @importFrom MASS glm.nb
#' @import mice
#'
#' @export
#'
#' @seealso \code{\link[MASS]{glm.nb}}; \code{\link[mice]{is.mice}}, \code{\link[mice]{pool}}.
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
#' ## Create design matrix
#' mydmat <- matrix(c(rep(1, 5), 1:5, rep(median(df$v2), 5)), ncol = 3)
#' mydmat <- cbind(mydmat, mydmat[,2] * mydmat[,3])
#' colnames(mydmat) <- c('(Intercept)', 'v1', 'v2', 'v1:v2')
#'
#' calc_nb_counts(mymod, designMatrix = mydmat, predVar = 'v1', intVar = 'v2')
#'

calc_nb_counts <- function(nbObj,
                           designMatrix,
                           predVar = NULL,
                           intVar = NULL){

  is.mice <- inherits(nbObj, 'mira')

  ## Model object must be fit using glm.nb()
  if(!(inherits(nbObj, 'negbin') | (is.mice & inherits(nbObj$analyses[[1]], 'negbin')))){
    stop('nbObj must be of class negbin from glm.nb(), or a mice object using glm.nb to fit',
         call. = FALSE)
  }

  ## Get coefficients and vcov matrix, depending on whether object is imputed via mice or not
  if(is.mice){
    nb.coefs <- pool(nbObj)$qbar
    nb.vcov <- pool(nbObj)$t
  } else{
    nb.coefs <- coef(nbObj)
    nb.vcov <- vcov(nbObj)
  }
  coefnames <- names(nb.coefs)

  if(sum(!(coefnames %in% colnames(designMatrix))) > 0){
    stop("Variables in nbObj are not in designMatrix", call. = FALSE)
  }

  designMatrix <- designMatrix[,coefnames]

  ## Calculate linear predictors and their SEs
  lp <- apply(designMatrix, MARGIN = 1, FUN = function(x){ sum(nb.coefs * as.numeric(x)) })
  lp.se <- apply(designMatrix, MARGIN = 1, FUN = function(x){ sqrt(t(x) %*% nb.vcov %*% x) })

  ## Calculate LCL, UCLs for linear predictors
  lp.lcl <- lp - qnorm(0.975)*lp.se
  lp.ucl <- lp + qnorm(0.975)*lp.se

  ## Calculate predicted counts, CIs as exp(quantities)
  count.pe <- exp(lp)
  count.lcl <- exp(lp.lcl)
  count.ucl <- exp(lp.ucl)

  ## Bind all results into data frame for plotting
  if(!is.null(predVar)){
    if(!(predVar %in% colnames(designMatrix))){
      stop("Predictor variable name (predVar) not in column names of design matrix", call. = FALSE)
    }

    ## Pull out predictor values (only works for now if predictor is linear and continuous)
    xvalue <- designMatrix[,predVar]

    ## If intVar is supplied, do the same for it; also only works if intVar is linear & continuous
    if(!is.null(intVar)){
      intvalue <- designMatrix[,intVar]
      return(as.data.frame(cbind(xvalue, intvalue, lp, lp.se, count.pe, count.lcl, count.ucl)))
    } else{
      return(as.data.frame(cbind(xvalue, lp, lp.se, count.pe, count.lcl, count.ucl)))
    }
  } else{
    return(as.data.frame(cbind(lp, lp.se, count.pe, count.lcl, count.ucl)))
  }
}
