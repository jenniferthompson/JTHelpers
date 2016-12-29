#' Create count process data set, eg, for use in time-dependent Cox proportional hazards models.
#'
#' Note that Zhiguo (Alex) Zhao wrote the majority of this code, with later edits by Cole Beck and
#' minor ones by me. It is here for my convenience, since Alex shared it with me and to my knowledge
#' it is not available elsewhere, not in an attempt to take credit.
#'
#' @param org.data Original data frame with multiple records per patient
#' @param id.var Variable name for subject identifier (eg, ID)
#' @param record.var Variable name for record identifier (eg, study day)
#' @param time.var Variable indicating end of records (eg, day of death)
#' @param event.var Variable indicating whether event happened (eg, "died.yn")
#' @param event.string Level of event variable indicating whether event happened (eg, "yes")
#' @param data.set Data set (one record/pt) containing time/event variables
#' @param out.strings Labels for 0/1 values of outcome variable
#'
#' @return data.frame containing reference, comparison, effect, lower and upper confidence limits,
#' variable name and indicator for whether row contains reference:reference comparison.
#'
#' @import rms
#' @importFrom dplyr mutate select
#' @importFrom tidyr separate
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @seealso \href{https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf}{Terry
#' Therneau's vignette on time-dependent Cox models}, \code{\link[survival]{survival}}.
#'
#' @examples

## -- Function to create time-varying Cox data -- ##
## -- Original data set requires "id", "study.day" variables -- ##
create_countprocess_data <- function(org.data,
                                     id.var,
                                     record.var,
                                     time.var,
                                     event.var,
                                     event.string,
                                     data.set,
                                     out.strings = c('Alive through end of interval',
                                                     'Died at end of interval')){
  ## Get vector of unique IDs, create null list with one slot for each
  ids <- unique(org.data[, id.var])
  final.data <- vector("list", length(ids))
  cols2drop <- match(c(id.var, record.var), names(org.data))

  ## For each ID...
  for(i in seq_along(ids)) {
    cur.id <- ids[i]

    ## Get all relevant daily data for specific patient
    use.data <- org.data[org.data[,id.var] == cur.id, -cols2drop]

    ## find rows with missing data
    valid.rows <- apply(use.data, MARGIN=1, FUN=function(i) !any(is.na(i)))

    ## Create vector of study days with non-missing data and first day with non-missing data
    study.day.id <- org.data[org.data[,id.var] == cur.id, record.var][valid.rows]
    first.day <- min(study.day.id)

    ## Delete any rows with missing data
    use.data <- use.data[valid.rows,]

    ## For each ID, cycle through each record. If current row = last row, increment counter by one.
    nr <- nrow(use.data)
    if(nr > 1) {
      stop.day <- rep(NA, nrow(use.data))
      tmp <- apply(use.data, MARGIN=1, paste, collapse='')
      for(j in seq(nr-1)) {
        if(tmp[j] != tmp[j+1]) {
          stop.day[j] <- study.day.id[j]
        }
      }
    } else {
      stop.day <- NA
    }

    ## Cbind each stop day to daily data set, and add date of death/end of period to last row
    stop.day[length(stop.day)] <- data.set[data.set[,id.var] == cur.id, time.var]
    use.data <- cbind(use.data, stop.day)

    ## Subset to only non-repeated rows, create start date vector
    use.data <- use.data[!is.na(stop.day),]
    use.data$start.day <- c((first.day - 1), use.data$stop.day[-nrow(use.data)])
    nr <- nrow(use.data)

    ## If patient died or discharged on a day with assessments, or on the day of enrollment,
    ##  assign interval half a day
    use.data$stop.day <- ifelse(use.data$stop.day == use.data$start.day,
                                use.data$stop.day + 0.5,
                                use.data$stop.day)

    ## Create event vector: all 0 unless patient died, in which case last value = 1
    use.data$died <- 0
    use.data$died[nr] <- (data.set[data.set[,id.var] == cur.id, event.var] == event.string)*1

    final.data[[i]] <- cbind(cur.id, use.data)

    cat(sprintf("\rFinished %s", cur.id))
  }
  cat("\n")

  ## combine list of data.sets into one
  final.data <- do.call('rbind', final.data)

  final.data$died <- factor(final.data$died, levels = 0:1, labels = out.strings)
  names(final.data) <- gsub('cur.id', id.var, names(final.data), fixed = TRUE)
  return(final.data)
}
