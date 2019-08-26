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
#' @importFrom dplyr mutate
#' @export
#'


#NOTE: I'm not sure why the dplyr pipes work. I assume it's because I have it in the environment from
# the run_r_proj.R script. I'll need to look up how to get pipes in the package without calling dplyr.

bind_geo.ave_grid <- function(geo.ave_data, gridded_data){
    assertthat::assert_that(all((class(geo.ave_data) == "griddata"), class(gridded_data) == "griddata"))
    assertthat::assert_that(nrow(geo.ave_data$vardata) == nrow(gridded_data$vardata))

    data.frame(global_value = geo.ave_data$vardata) %>%
        bind_cols(data.frame(gridded_data$vardata)) %>%
        gather(grid_cell,cell_value, -global_value) %>%
        mutate(grid_cell = as.integer(gsub("X","",grid_cell))) %>%
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

std_normalize <- function(df){
    df %>%
        mutate(global_value = (global_value - mean(global_value))/sd(global_value),
               cell_value = (cell_value - mean(cell_value))/sd(cell_value)) %>%
        return()
}

