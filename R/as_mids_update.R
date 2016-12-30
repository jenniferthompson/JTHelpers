#' Slightly edited version of as.mids() from mice package
#'
#' The original function relies on data being sorted by .imp to work properly. This is problematic
#' in data sets created manually - say, calculating summary statistics using original + imputed
#' data, where the final data set we want to combine does not have an .imp column. Original function
#' by Stef van Buuren, with detective work and changes by Cole Beck. See
#' \href{https://github.com/jenniferthompson/TBIMort/blob/master/tbi_mortality_riskfactors.Rmd}{this
#'  analysis} for a full example.
#'
#' @param data A multiply imputed data set in long format
#' @param .imp Mandatory column indicator for the multiple imputation stream, where 0 indicates the
#' incomplete data and 1 through m indicate the m multiple imputation streams. Default is 1.
#' @param .id Optional column indicator for the row numbers. Default is 2.
#'
#' @import mice
#'
#' @export
#'

as.mids.update <- function(data, .imp = 1, .id = 2){
  imps <- data[, .imp]
  if(is.factor(imps)){
    m = max(as.numeric(levels(imps))[imps])
  } else {
    m = max(imps)
  }

  ini <- mice(data[imps == 0, -c(.imp, .id)], m = m, maxit = 0)
  names  <- names(ini$imp)
  if (!is.null(.id)){
    rownames(ini$data) <- data[imps == 0, .id]
  }
  for (i in 1:length(names)){
    ## This section is tweaked from original
    if(!is.null(ini$imp[[i]])){
      missix <- is.na(data[imps == 0, names[i]])
      for(j in 1:m){
        indic <- which(imps == j)[missix]
        ini$imp[[names[i]]][j] <- data[indic, names[i]]
      }
    }
  }
  return(ini)
}
