#' Combine information from summary.rms(), anova.rms(), and other rms object info to create a
#' single tidy data.frame of model results that could be printed in a report.
#'
#' @param rmsObj Object of class rms
#' @param labeldf Data frame from which to get variable labels. Defaults to NULL, in which case
#' variable names will be used.
#' @param shortLabs Named vector of variable labels to replace in interaction rows. Must be in
#' format c("variable name" = "shortened label").
#' @param rndDigits Number of digits to round reference, comparison, result and CI values to.
#' Defaults to 2.
#' @param rndTestStat Number of digits to round test statistic to. Defaults to 1.
#' @param rndModInfo Number of digits to round model LR, R2, etc to. Defaults to rndDigits.
#' @param printFormat How to format results for printing. Defaults to plain text; LaTeX or Markdown
#' available.
#' @param rmRows Which sets of rows to exclude from final model results. For example, it is more
#' appropriate to look at the test for total nonlinearity to determine whether to keep these terms
#' than to look at individual variables' tests for nonlinearity. Example 2: in cases like nonlinear
#' interactions, some rows of output can be redundant. Options: nl (removes rows labeled
#' Nonlinear [Interaction] (...) in anova.rms()); int (removes rows labeled All Interactions in
#' anova.rms()); nl.int (removes all of the above); none (keeps all rows). Defaults to nl.int.
#' @param whichInfo Vector of elements of rmsObj$stats to include. Defaults to model likelihood
#' ratio, R^2, C, Dxy, as available. (Model df are included in the main table.)
#'
#' @import rms
#' @importFrom Hmisc contents label
#'
#' @return data.frame containing quantity descriptions, reference and comparison values,
#' coefficient or ratio and 95% CI, test statistic, df, and p-value, along with observations and
#' model statistics.
#'
#' @seealso \code{\link{rms::ols}}, \code{\link{[rms]lrm}}, \code{\link{[rms]cph}},
#' \code{\link{[rms]Gls}}, \code{\link{[rms]summary.rms}}, \code{\link{[rms]anova.rms}}.
#'
#' @examples
#'
#' ## Fit linear regression model using ols()
#' mymod <- ols(mpg ~ cyl * wt, data = mtcars)
#'
#' ## Make sure to set datadist; summary.rms() requires this
#' dd.mtcars <- datadist(mtcars)
#' options(datadist = 'dd.mtcars')
#'
#' ## Defaults, using variable names
#' rms_model_results(mymod)
#'
#' ## Defaults, using variable labels
#' label(mtcars$cyl) <- 'Number of cylinders'
#' label(mtcars$wt) <- 'Weight (1000 lbs)'
#' rms_model_results(mymod, labeldf = mtcars)
#'
#' ## Use shortened variable labels for interaction row
#' rms_model_results(mymod, labeldf = mtcars, shortLabs = c('cyl' = 'Cylinders', 'wt' = 'Weight'))
#'
#' ## Use shortened variable labels, format for Markdown
#' rms_model_results(mymod, labeldf = mtcars, shortLabs = c('cyl' = 'Cylinders', 'wt' = 'Weight'), printFormat = 'markdown')
#'
rms_model_results <- function(rmsObj,
                              labeldf = NULL,
                              shortLabs = NULL,
                              rndDigits = 2,
                              rndTestStat = 1,
                              rndInfo = rndDigits,
                              printFormat = c('plain', 'markdown', 'latex'),
                              rmRows = c('nl.int', 'nl', 'int', 'none'),
                              addInfo){

  ## Argument type tests
  if(!(inherits(rmsObj, 'rms'))){
    stop('rmsObj must be an rms model fit')
  }

  if(!is.null(labeldf) & !(inherits(labeldf, 'data.frame'))){
    stop('labeldf must be a data.frame')
  }

  ## -- General setup ----------------------------------------------------------------------------
  ## Which rows (if any) to remove from final?
  rmRows <- match.arg(rmRows)

  ## Set printFormat if not set already; create indicator for Latex formatting to save typing
  printFormat <- match.arg(printFormat)
  formatLatex <- printFormat == 'latex'
  formatMd <- printFormat == 'markdown'

  ## Spacing: if not formatted for Latex, use three spaces; if formatted for latex, use ~~~
  space.str <- ifelse(formatLatex, '~~~',
               ifelse(formatMd, paste(rep('&nbsp;', 3), collapse = ''), '   '))

  ## These model classes will have ratios presented, not beta coefficients
  use.ratios <- c('lrm', 'cph')

  ## -- Variable labels ----------------------------------------------------------------------------
  ## If data.set is NULL or has no labels, use variable names
  label.data <- NULL
  if(!is.null(labeldf)){
    if('Labels' %in% names(contents(labeldf)$contents)){
      label.data <- data.frame(variable = names(labeldf),
                               varlabel = as.vector(contents(labeldf)$contents$Labels),
                               stringsAsFactors = FALSE)

      ## Create indicator for whether shortened labels will be used; if so, create new column with
      ## original labels, replaced with short versions as indicated by names of short.labels
      use.short <- FALSE

      ## Remove and warn of any elements of shortLabs that aren't in names(labeldf)
      if(length(setdiff(names(shortLabs), names(labeldf))) > 0){
        message(paste("Note: these elements of shortLabs do not appear in names(labeldf):",
                      paste(setdiff(names(shortLabs), names(labeldf)), collapse = ', ')))
      }

      shortLabs <- shortLabs[names(shortLabs) %in% names(labeldf)]
      if(length(shortLabs) > 0){
        use.short <- TRUE

        ## All short labels are the same as variable labels by default
        label.data$shortlabel <- label.data$varlabel

        ## Replace desired variable labels with short versions
        for(i in 1:length(shortLabs)){
          varnum <- match(names(shortLabs)[i], label.data$variable)
          if(!is.na(varnum)){
            label.data$shortlabel[varnum] <- shortLabs[i]
          }
        }
      }
    }
  }

  ## -- summary.rms() results --------------------------------------------------------------------
  mod.sum <- as.data.frame(summary(rmsObj)) ## Often won't print due to duplicate rownames

  ## Get variable name in every row
  mod.sum$quantity <- rownames(mod.sum)
  rownames(mod.sum) <- NULL

  ## Create variable column ##
  ## Models for which summary() produces both coefficients and ratios: Take only ratios, variable
  ## column = row above ratio row
  if(inherits(rmsObj, use.ratios)){
    mod.sum$var <- c(NA, mod.sum$quantity[1:(nrow(mod.sum) - 1)])
    mod.sum <- subset(mod.sum, Type == 2)
  } else{
    mod.sum$var <- mod.sum$quantity

    ## If model is a Poisson model, exponentiate all point estimates, CLs to get IRRs
    if(inherits(rmsObj, 'Glm')){
      if(model.obj$family$family == 'poisson'){
        mod.sum[,c('Effect', 'Lower 0.95', 'Upper 0.95')] <-
          exp(mod.sum[,c('Effect', 'Lower 0.95', 'Upper 0.95')])
      }
    }
  }

  ## For categorical covariates, var column currently is of format
  ##  "variable name - comparison category:reference category"; split out components
  var.split <- lapply(mod.sum$var, FUN = function(x){ strsplit(x, split = ' - |:')[[1]] })
  mod.sum$var <- unlist(lapply(var.split, FUN = function(x){ x[1] }))

  ## Create new low/high variables with categories (not factor numbers) or formatted numbers
  mod.sum$low.char <- unlist(lapply(var.split, FUN = function(x){ x[3] }))
  mod.sum$low.char <- with(mod.sum, ifelse(is.na(low.char),
                                           format(round(Low, rndDigits), nsmall = rndDigits),
                                           as.character(low.char)))
  mod.sum$high.char <- unlist(lapply(var.split, FUN = function(x){ x[2] }))
  mod.sum$high.char <- with(mod.sum, ifelse(is.na(high.char),
                                            format(round(High, rndDigits), nsmall = rndDigits),
                                            as.character(high.char)))

  ## -- anova.rms() results ----------------------------------------------------------------------
  ## Create data frame of anova() results, add column to indicate variable/quantity
  mod.anova <- as.data.frame(anova(rmsObj))
  mod.anova$line <- rownames(mod.anova); rownames(mod.anova) <- NULL
  mod.anova$var <- gsub('^ ', '', gsub(' +\\(.*\\)$| : .*$', '', mod.anova$line))

  ## Remove rows requested in rmRows
  if(rmRows != 'none'){
    remove.these <- NULL
    if(rmRows %in% c('nl', 'nl.int')){
      remove.these <- c(remove.these, '^ *Nonlinear', '^ *f\\(A,B\\)')
    }
    if(rmRows %in% c('int', 'nl.int')){
      remove.these <- c(remove.these, '^ *All Interactions', '^ *Nonlinear Interaction')
    }

    remove.these.rows <- unique(unlist(sapply(remove.these, grep, mod.anova[,'var'])))
    if(length(remove.these.rows) > 0){
      mod.anova <- mod.anova[-remove.these.rows,]
    }
  }

  ## For each row, create unique var value that is as descriptive as possible
  cur.var <- NULL
  for(i in 1:nrow(mod.anova)){
    if((!missing(labeldf) & mod.anova$var[i] %in% names(labeldf)) |
       length(grep('*', mod.anova$var[i], fixed = TRUE)) > 0){
      cur.var <- mod.anova$var[i]
    } else{
      if(length(grep('^TOTAL|ERROR', mod.anova$var[i])) == 0){
        mod.anova$var[i] <- ifelse(is.null(cur.var),
                                   gsub('^All ', '', mod.anova$var[i]),
                                   paste(cur.var, gsub('^All ', '', mod.anova$var[i])))
      }
    }
  }

  ## -- Combine anova, summary results and format for printing -----------------------------------
  ## Merge data frames
  mod.data <- merge(subset(mod.sum, select = -c(Low, High, Type, quantity)),
                    subset(mod.anova, select = -c(line)),
                    by = 'var', all = TRUE)

  ## Replace variable name for test statistic with something without punctuation
  names(mod.data) <- gsub('^Chi-Square|F$', 'test.stat', names(mod.data))

  ## Create count variable for each line for a given quantity
  ##  (eg, multiple levels of categorical covariate)
  mod.data$var.line <- unlist(lapply(unique(mod.data$var),
                                     FUN = function(x){1:nrow(subset(mod.data, var == x))}))

  ## Format columns which need to be rounded
  format.cols <- c('Effect', 'Lower 0.95', 'Upper 0.95')
  for(i in 1:length(format.cols)){
    colname <- format.cols[i]
    mod.data[,colname] <- ifelse(!is.na(as.numeric(as.character(mod.data[,colname]))),
                                 format(round(as.numeric(as.character(mod.data[,colname])),
                                              rnd.digits),
                                        nsmall = rnd.digits),
                                 as.character(mod.data[,colname]))
    mod.data[,colname] <- ifelse(is.na(mod.data[,colname]), '', mod.data[,colname])
  }

  mod.data$test.stat <- with(mod.data, {
    ifelse(!is.na(as.numeric(as.character(test.stat))),
           format(round(as.numeric(as.character(test.stat)), rndTestStat), nsmall = rndTestStat),
           as.character(test.stat)) })
  mod.data$test.stat <- ifelse(is.na(mod.data$test.stat), '', mod.data$test.stat)

  ## Format p-value to 3 places unless <0.001
  mod.data$pval <- with(mod.data, ifelse(is.na(P), '',
                                  ifelse(P < 0.001, '<0.001', format(round(P, 3), nsmall = 2))))

  ## Replace missing values for low, high columns with ''
  mod.data$low.char <- with(mod.data, ifelse(is.na(low.char), '', low.char))
  mod.data$high.char <- with(mod.data, ifelse(is.na(high.char), '', high.char))

  ## Replace column names with easier-to-work-with names
  mod.data <- subset(mod.data,
                     select = c('var', 'low.char', 'high.char', 'Effect', 'Lower 0.95',
                                'Upper 0.95', 'test.stat', 'd.f.', 'pval', 'var.line'))
  names(mod.data) <- c('var', 'low', 'high', 'effect', 'lcl', 'ucl', 'test.stat', 'df', 'pval',
                       'var.line')

  ## Create string for estimate (CI) results
  mod.data$est.ci <-
    with(mod.data, ifelse(effect == '', '',
                          paste0(effect, ' (',
                                 gsub(' ', '', lcl), ', ',
                                 gsub(' ', '', ucl), ')')))

  ## -- Format for printing: variable labels, sorting --------------------------------------------
  ## Create variable for sort order and resort data
  mod.data$sort <- rep(NA, nrow(mod.data))
  mod.data$sort[grep('^ERROR$', mod.data$var)] <- 7
  mod.data$sort[grep('^TOTAL$', mod.data$var)] <- 6
  mod.data$sort[grep('^TOTAL INTERACTION$', mod.data$var)] <- 5
  mod.data$sort[grep('^TOTAL NONLINEAR$', mod.data$var)] <- 4
  mod.data$sort[grep('^TOTAL NONLINEAR \\+ INTERACTION$', mod.data$var)] <- 3
  mod.data$sort[grep('*', mod.data$var, fixed = TRUE)] <- 2
  mod.data$sort[is.na(mod.data$sort)] <- 1

  mod.data <- mod.data[with(mod.data, order(sort, var, var.line)),]

  ## Function to determine label value given a string
  create.label <- function(v){
    ## For rows involving interactions, replace variable names with labels (if using) and keep
    ##  all other default descriptions (eg, Nonlinear Interactions, f(A,B) vs. Af(B) + Bg(A))
    int.label.str <- NULL

    if(length(grep('*', v, fixed = TRUE)) > 0){
      ## Separate all pieces of interaction description; replace variable names with labels
      int.terms <- unlist(lapply(strsplit(v, ' ')[[1]], FUN = function(x){
        if(x %in% names(labeldf)){
          if(use.short){
            labelcol <- 'shortlabel'
          } else{
            labelcol <- 'varlabel'
          }

          ifelse(x %in% label.data[,'variable'],
                 label.data[match(x, label.data$variable), labelcol],
                 x)
        } else{
          x
        }
      }))

      ## Recombine all elements of interaction row description
      int.label.str <- paste(int.terms, collapse = ' ')
    }

    ifelse(formatLatex & v == 'TOTAL', '\\emph{\\textbf{Overall Model}}',
    ifelse(formatMd & v == 'TOTAL', '***Overall Model***',
    ifelse(v == 'TOTAL', 'Overall Model',
    ifelse(formatLatex & v == 'TOTAL NONLINEAR', '\\emph{All Nonlinear Terms}',
    ifelse(formatMd & v == 'TOTAL NONLINEAR', '*All Nonlinear Terms*',
    ifelse(formatLatex & v == 'TOTAL NONLINEAR + INTERACTION',
           '\\emph{All Nonlinear \\& Interactions}',
    ifelse(formatMd & v == 'TOTAL NONLINEAR + INTERACTION', '*All Nonlinear & Interactions*',
    ifelse(formatLatex & v == 'TOTAL INTERACTION', '\\emph{All Interaction Terms}',
    ifelse(formatMd & v == 'TOTAL INTERACTION', '*All Interaction Terms*',
    ifelse(formatLatex & v == 'ERROR', '\\emph{Error}',
    ifelse(formatMd & v == 'ERROR', '*Error*',
    ifelse(v == 'ERROR', 'Error',
    ifelse(v == 'TOTAL NONLINEAR', 'All Nonlinear Terms',
    ifelse(v == 'TOTAL NONLINEAR + INTERACTION', 'All Nonlinear & Interaction Terms',
    ifelse(v == 'TOTAL INTERACTION', 'All Interaction Terms',
    ifelse(!is.null(int.label.str), int.label.str,
    ifelse(length(grep('interaction[s]*$', tolower(v))) > 0, paste0(space.str, 'Interactions'),
    ifelse(length(grep('Nonlinear|NONLINEAR$', v)) > 0, paste0(space.str, 'Nonlinear'),
    ifelse(!is.null(label.data) & v %in% label.data$variable,
           label.data$varlabel[match(v, label.data$variable)],
           v)))))))))))))))))))
    }

  mod.data$label <- unlist(lapply(1:nrow(mod.data), FUN = function(i){
    ifelse(mod.data$var.line[i] != 1, '', create.label(mod.data$var[i]))
  }))

  ## Remove test stat, df, p for lines for multiple categories of variable with >2 levels
  mod.data$test.stat <- with(mod.data, ifelse(var.line != 1, '', test.stat))
  mod.data$df <- with(mod.data, ifelse(var.line != 1, '', df))
  mod.data$pval <- with(mod.data, ifelse(is.na(pval) | var.line != 1, '', pval))

  ## Select only desired components of data frame in correct order
  mod.data <- subset(mod.data, select = c(label, low, high, est.ci, test.stat, df, pval))

  ## -- Add model stats of interest ----------------------------------------------------------------
  if(missing(addInfo)){
    addInfo <- c('Model L.R.', 'C', 'Dxy', 'R2')
  }
  use.stats <- intersect(addInfo, names(rmsObj$stats))
  mod.stats <- data.frame(label = use.stats, low = rmsObj$stats[use.stats],
                          stringsAsFactors = FALSE)

  ## Format statistics for printing
  mod.stats$low <- with(mod.stats, {
    ifelse(label == 'd.f.',
           format(round(low), nsmall = 0),
           format(round(low, digits = rndInfo), nsmall = rndInfo)) })
  if(formatLatex){
    mod.stats$label <- ifelse(mod.stats$label == 'R2',
                              '\\emph{$R^{2}$}',
                              paste0('\\emph{', mod.stats$label, '}'))
  } else if(formatMd){
    mod.stats$label <- ifelse(mod.stats$label == 'R2', '*R^2^*', paste0('*', mod.stats$label, '*'))
  }

  ## Add number of observations, and events for cph() (separately, in order to put at the end)
  n.rows <- data.frame(label = grep('^n|Obs|Events', names(rmsObj$stats), value = TRUE),
                       low = rmsObj$stats[grep('^n|Obs|Events', names(rmsObj$stats))])
  n.rows$label <- gsub('Obs', 'Observations', n.rows$label)

  if(formatLatex){
    n.rows$label <- paste0('\\emph{', n.rows$label, '}')
  } else if(formatMd){
    n.rows$label <- paste0('*', n.rows$label, '*')
  }

  mod.stats <- rbind(mod.stats, n.rows)

  ## Add frequencies of individual categories, if applicable
  if('freq' %in% names(rmsObj)){
    freq.data <- data.frame(label = paste0(space.str, names(rmsObj$freq)),
                            low = as.numeric(rmsObj$freq))
    mod.stats <- rbind(mod.stats, freq.data)
  }

  mod.stats$high <- mod.stats$est.ci <- mod.stats$test.stat <- mod.stats$df <- mod.stats$pval <-
    rep('', nrow(mod.stats))

  ## Add stats to main data frame
  mod.data <- rbind(mod.data, rep('', ncol(mod.data)), mod.stats)
  rownames(mod.data) <- NULL

  return(mod.data)
}
