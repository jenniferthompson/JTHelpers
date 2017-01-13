#' Calculate frequencies and proportions for a categorical or factor variable, optionally using
#' specified denominator(s), with or without overall totals.
#'
#' This function was inspired by \code{Hmisc::describe}, but adding the flexibility to use a
#' vector of denominators vs \code(length(non-missing values)) to calculate proportions. For
#' example, we might have a multicenter study and want the proportion of patients excluded out of
#' those enrolled per site.
#'
#' @param descVar Vector to be described.
#' @param allLevels Vector with all levels which should be present. Defaults to unique(descVar).
#' @param overall Logical; whether to provide overall count/%.
#' @param changeDenom Logical; whether to use a different denominator than length(descVar).
#' @param useDenom Numeric vector (length = length(allLevels)); denominators to be used if
#' changeDenom = TRUE.
#' @param pctString Logical; whether to include a second element of strings of form "N (Pct%)".
#'
#' @return List with 1 or 2 elements: always $freqNums, and $freqStrings if pctString = TRUE.
#'
#' @export
#'
#' @seealso \code{\link[Hmisc]{describe}}.
#'
#' @examples
#'
#' ## Create example vector
#' x <- sample(LETTERS[1:5], size = 50, replace = TRUE)
#'
#' ## Calculate frequencies and proportions using defaults
#' calc_cat_freqs(x, pctString = TRUE)
#'
#' ## But maybe we had possibilities for F and G levels too, they just don't show up
#' calc_cat_freqs(x, allLevels = LETTERS[1:7])
#'
#' ## Alternatively, maybe we don't care about E
#' calc_cat_freqs(x, allLevels = LETTERS[1:4])
#'
#' ## But above gives us % of total non-missing values of descVar; we want the total of just A:D
#' calc_cat_freqs(x, allLevels = LETTERS[1:4], changeDenom = TRUE, useDenom = sum(x %in% LETTERS[1:4]))
#'
#' ## Maybe A-E are sites and x represents a patient excluded from each; we want the denominator to
#' ## be number screened at each site
#' ## Let us pretend this is the number screened:
#' nScreened <- c(100, 110, 120, 130, 140)
#' calc_cat_freqs(x, changeDenom = TRUE, useDenom = nScreened)
#'

calc_cat_freqs <- function(descVar,
                           allLevels,
                           overall = TRUE,
                           changeDenom = FALSE,
                           useDenom,
                           pctString = FALSE){

  ## Test: descVar must be character or numeric
  if(!inherits(descVar, c('character', 'factor'))){
    stop('descVar must be a character or factor vector', call. = FALSE)
  }

  ## allLevels: Warning or message if allLevels and unique(descVar) are not the same; set default
  ## if allLevels is not provided
  if(!missing(allLevels)){
    if(!all(unique(descVar) %in% allLevels)){
      notInLevels <- setdiff(unique(descVar), allLevels)
      warning(paste('Not all values of descVar are in allLevels. The following values will be dropped:',
                    paste(sort(notInLevels), collapse = '; ')),
              call. = FALSE)
    }

    if(!all(allLevels %in% unique(descVar))){
      message(paste('Not all values of allLevels are in descVar. The following values will be added with N = 0:',
                    paste(sort(notInLevels), collapse = '; ')))
    }
  } else{
    allLevels <- unique(descVar)
  }

  ## useDenom: if changeDenom = TRUE, must be given and must be either length 1 or
  ## length = length(allLevels); otherwise, use length(descVar)
  if(missing(useDenom)){
    if(changeDenom){
      stop('If changeDenom = TRUE, useDenom must be given', call. = FALSE)
    }
    useDenom <- length(!is.na(descVar))
  } else{
    if(!(length(useDenom) == 1 | length(useDenom) == length(allLevels))){
      stop('length(useDenom) must be either 1 or length(allLevels)', call. = FALSE)
    }
  }

  ## Restrict descVar to only levels specified
  descVar <- descVar[descVar %in% allLevels]

  descTable <- table(descVar)
  levelNames <- names(descTable)
  levelNs <- as.numeric(descTable)

  ## If there are values in allLevels not represented in the data, add levels + 0s to levelNames,
  ## levelsNs, respectively
  if(!all(allLevels %in% levelNames)){
    notInTable <- setdiff(allLevels, levelNames)
    levelNames <- c(levelNames, notInTable)
    levelNs <- c(levelNs, rep(0, length(notInTable)))
  }

  ## Sort both vectors alphabetically
  levelNs <- levelNs[order(levelNames)]
  levelNames <- levelNames[order(levelNames)]

  ## Add Overall if requested
  if(overall){
    totalN <- sum(levelNs)
    levelNs <- c(levelNs, totalN)
    levelNames <- c(levelNames, 'Overall')
    useDenom <- c(useDenom, sum(useDenom))
  }

  ## Calculate proportions
  levelProps <- levelNs / useDenom

  levelPcts <- round(levelProps * 100)

  freqNums <- matrix(c(levelNs, levelPcts), nrow = 2, ncol = length(levelNs), byrow = TRUE)
  colnames(freqNums) <- levelNames

  if(pctString){
    freqStrings <- matrix(paste0(levelNs, ' (', levelPcts, '%)', sep = ''), nrow = 1)
    colnames(freqStrings) <- levelNames
    return(list("freqNums" = freqNums,
                "freqStrings" = freqStrings))
  } else{
    return(list("freqNums" = freqNums))
  }
}
