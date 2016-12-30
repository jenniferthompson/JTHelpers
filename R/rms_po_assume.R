#' Create figures to visually examine proportional odds assumption from an `lrm` model fit.
#'
#' Based on code and strategy outlined in Frank Harrell's Regression Modeling Strategies (see
#' citation below).
#'
#' @param rmsObj Model fit of class `lrm`.
#' @param cuts Numeric vector; sequence of points to cut outcome. Should not include lowest outcome
#' level.
#' @param plotVars Character vector; which variables to plot. Defaults to all.
#' @param mfrowAuto Logical; whether to determine par(mfrow = ...) automatically. Defaults to FALSE.
#' @param modelData Data.frame; data set used to fit original model. Used to fit logistic models
#' with outcome dichotomized at each cut point.
#'
#' @import rms
#'
#' @export
#'
#' @seealso Harrell FE. *Regression Modeling Strategies: with applications to linear models,
#' logistic regression, and survival analysis.* New York: Springer Science + Business Media, LLC,
#' 2001.
#'
#' @examples
#'
#' df <- data.frame(ptclass = sample(1:4, size = 20, replace = TRUE),
#'                  v1 = rnorm(n = 20),
#'                  v2 = rnorm(mean = 5, sd = 1, n = 20))
#'
#' mymod <- lrm(ptclass ~ v1 + v2, data = df)
#' rms_po_assume(mymod, cuts = 2:4, modelData = df)
#'

rms_po_assume <- function(lrmObj,          ## model.obj: model of class lrm()
                          cuts,               ## cuts: sequence of points to cut outcome
                          plotVars = NULL,   ## which variables to plot (defaults to all variables)
                          mfrowAuto = FALSE, ## whether to determine par(mfrow) automatically
                          modelData){        ## model.data: data used to fit model

  if(!(inherits(lrmObj, 'lrm'))){
    stop('lrmObj must be of class lrm', call. = FALSE)
  }

  ## Create data set for each coefficient in main model
  cof.names <- names(coef(lrmObj))
  all.rows <- 1:length(cof.names)
  int.rows <- grep('y>=', cof.names, fixed = TRUE)
  take.rows <- all.rows[all.rows %nin% int.rows]
  cof <- data.frame(var = cof.names[take.rows])

  ## Extract formula from model call
  comp.call <- as.character(formula(lrmObj))
  model.outcome <- comp.call[2]
  model.formula <- comp.call[3]

  for(k in 1:length(cuts)){
    cut.mod <-
      lrm(as.formula(paste0('as.numeric(', model.outcome, ' >= ', cuts[k], ') ~ ', model.formula)),
          data = modelData)

    cof.temp <- data.frame(var = names(coef(cut.mod)),
                           hold.place = coef(cut.mod))

    cof <- merge(cof, cof.temp, all.x = TRUE, all.y = FALSE)
    names(cof) <- gsub('hold.place', paste('coef.cut', cuts[k], sep = '.'), names(cof))
  }

  ## Subset in case some splines didn't make requested number of knots
  cof <- cof[rowSums(is.na(cof[,2:ncol(cof)])) == 0,]

  ## If plot.vars is not null, take only rows for variables of interest
  if(!is.null(plotVars)){
    cof <- cof[grep(plotVars, cof$var),]
  }

  if(mfrowAuto){
    ## Get number of rows/columns for plot (plot as close to square as possible)
    plot.rows <- ceiling(sqrt(nrow(cof)))
    par(mfrow = c(plot.rows, plot.rows), mar = c(2, 4, 1, 1))
  }

  for(k in 1:nrow(cof)){
    plot(cuts, cof[k, 2:ncol(cof)], type = 'l', ylab = '')
    title(main = model.outcome)
    title(ylab = cof[k, 'var'], line = 2.5)
    abline(h = 0, lty = 2, col = 'red')
  }
}
