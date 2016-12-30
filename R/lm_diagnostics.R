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

## -- Function to graphically examine proportional odds assumption ---------------------------------
po.assume <- function(model.obj,          ## model.obj: model of class lrm()
                      cuts,               ## cuts: sequence of points to cut outcome
                      plot.vars = NULL,   ## which variables to plot (defaults to all variables)
                      mfrow.auto = FALSE, ## whether to determine par(mfrow) automatically
                      model.data){        ## model.data: data used to fit model
  ## Create data set for each coefficient in main model
  cof.names <- names(coef(model.obj))
  all.rows <- 1:length(cof.names)
  int.rows <- grep('y>=', cof.names, fixed = TRUE)
  take.rows <- all.rows[all.rows %nin% int.rows]
  cof <- data.frame(var = cof.names[take.rows])

  ## Extract formula from model call
  comp.call <- as.character(formula(model.obj))
  model.outcome <- comp.call[2]
  model.formula <- comp.call[3]

  for(k in 1:length(cuts)){
    cut.mod <-
      lrm(as.formula(paste0('as.numeric(', model.outcome, ' >= ', cuts[k], ') ~ ', model.formula)),
          data = model.data)

    cof.temp <- data.frame(var = names(coef(cut.mod)),
                           hold.place = coef(cut.mod))

    cof <- merge(cof, cof.temp, all.x = TRUE, all.y = FALSE)
    names(cof) <- gsub('hold.place', paste('coef.cut', cuts[k], sep = '.'), names(cof))
  }

  ## Subset in case some splines didn't make requested number of knots
  cof <- cof[rowSums(is.na(cof[,2:ncol(cof)])) == 0,]

  ## If plot.vars is not null, take only rows for variables of interest
  if(!is.null(plot.vars)){
    cof <- cof[grep(plot.vars, cof$var),]
  }

  if(mfrow.auto){
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
