#' Read parsed netCDF data for arbitrary variables into dataframes to input
#' into regression models.
#'
#' The variable data should be a pair of parsed netCDF files, one of which should
#' be globally averaged the other should have globally gridded values.
#' The data can have any time series (daily, monthly, annually, etc.) as long as
#' both parsed files have the same length.
#'
#' The output will be a dataframe with three columns:
#' \describe{
#'   \item{\strong{global_value}}{Globally averaged value for whichever variable was fed in}
#'   \item{\strong{grid_cell}}{A unique integer labeling each grid cell.}
#'   \item{\strong{cell_value}}{The value of the respective grid cell and the corresponding
#'   globally averaged value when the cell held this value}
#'
#' Each row corresponds to a time step, however these are not unique as the table is in long format.
#' For example if 1 year of monthly data was read in, then after every 10th line the
#' time value would loop back to the start. This restart would correspond to the change
#' in grid_cell ID.
#'
#' @param geo.ave_data The globally averaged data
#' @param gridded_data The globally gridded data
#' @return A long form \code{data.frame}.
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols mutate select %>%
#' @importFrom tidyr gather
#' @export
#'

gather_griddata <- function(griddata){
    assertthat::assert_that(class(griddata) == "griddata")

    griddata$vardata %>%
        data.frame() %>%
        tidyr::gather(grid_cell, value) %>%
        dplyr::mutate(grid_cell = as.integer(gsub("X", "", grid_cell)),
                      lat = rep(griddata$lat, length(griddata$lon), each = length(griddata$time)),
                      lon = rep(griddata$lon, each = length(griddata$lat)*length(griddata$time)))%>%
        return()
}

#Current code simply numbers the month
label_month <- function(monthly_df){
    monthly_df %>%
        mutate(row_num = as.integer(rownames(test_frame_monthly)),
               month = rep(seq(1,12), nrow(test_frame_monthly)/12)) %>%
        select(-row_num) %>%
        return()
}

#' Standardize data to zero mean and unit variance
#'
#' A single standardization is performed over the entire dataset, as
#' opposed to a bunch of group-wise standardizations (e.g., with the
#' groups being grid cells).  This preserves relationships between grid cells,
#' making it easier to interpret coefficients.
#'
#' @param df A data frame of values to standardize
#' @param stvars List of names of variables to standardize in the df
#' @return Data frame of standardized values.  The original means and standard
#' deviations will be added as attributes \code{means} and \code{stdevs}.
#' @export
std_normalize <- function(df, stvars = c('global_value', 'cell_value')) {

    means <- sapply(stvars, function(var) {mean(df[[var]])})
    stdevs <- sapply(stvars, function(var) {sd(df[[var]])})
    names(means) <- names(stdevs) <- stvars
    for(var in stvars) {
        df[[var]] <- (df[[var]] - means[[var]]) / stdevs[[var]]
    }
    attr(df, 'means') <- means
    attr(df, 'stdevs') <- stdevs
    df
}

