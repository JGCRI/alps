## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
library(alps, quietly = TRUE)
library(rethinking, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)

tempo_aggregated_tas <- read.general(file.path("..", "data", "temperature", "tas_Amon_MIROC5_tempo_avg.nc"), varname = "tas")
geo_tempo_aggregated_tas <- read.geo.aggregate(file.path("..", "data", "temperature", "tas_Amon_MIROC5_geo_tempo_avg.nc"), varname = "tas")

first_row <- 1
last_row <- 100

assertthat::assert_that(all(last_row>first_row,
                            last_row<length(tempo_aggregated_tas$time)))

#Format the griddata into a table and bind global temperatures to downscaled data
long_tempo_aggregated_tas <- gather_griddata(tempo_aggregated_tas) %>% group_by(grid_cell) %>% slice(first_row:last_row) %>% ungroup()
long_tempo_aggregated_tas %>%
    mutate(global_value = geo_tempo_aggregated_tas$vardata[first_row:last_row,] %>%
               rep(length(unique(long_tempo_aggregated_tas$grid_cell)))) %>%
    rename(cell_value = value) ->
    temp_df_annual


read.area(file.path("..", "data", "land_fraction", "sftlf_fx_MIROC5_historical_r0i0p0.nc"), varname = "sftlf") %>%
    gather_griddata() %>%
    rename(land_fraction = value) ->
    land_sea_map

land_sea_temp_df_annual <- left_join(temp_df_annual, land_sea_map, by = c("grid_cell", "lat", "lon"))

#Test that there was no mismatch between grid/lat/lon for the temp df and land map
assertthat::assert_that(sum(is.na(land_sea_temp_df_annual$land_fraction))==0)



## ------------------------------------------------------------------------
all_lat <- unique(land_sea_temp_df_annual$lat)
all_lon <- unique(land_sea_temp_df_annual$lon)

lat_lon <- expand.grid(all_lat, all_lon) %>% rename(lat = Var1, lon = Var2)

#Group cells into n x m grids
x_len <- 2
y_len <- 2

assertthat::assert_that(all(length(all_lon)%%x_len == 0,
                            length(all_lat)%%y_len == 0))

lat_lon %>%
    mutate(lat_1 = as.integer(as.factor(lat)), lon_1 = as.integer(as.factor(lon))) %>%
    mutate(group = floor((lat_1+1)/y_len)+
               floor((lon_1-1)/x_len)*length(all_lon)/x_len ) %>%
    select(lat, lon, group) ->
    group_map

land_sea_temp_df_annual %>%
    left_join(group_map, by = c("lat", "lon")) ->
    land_sea_temp_df_annual

land_sea_temp_df_annual %>%
    filter(land_fraction > 66.67) ->
    land_temp_df_annual

land_sea_temp_df_annual %>%
    filter(land_fraction <= 66.67) ->
    sea_temp_df_annual




## ------------------------------------------------------------------------
#Select which parallel band of cells to look at
batch_size_list <- c(2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30, 40, 48, 60, 120, 240)
# parallel <-  51.1278666

land_temp_df_annual %>%
    # filter(between(lat, parallel-.1, parallel+.1)) %>%
    select(grid_cell, global_value, cell_value) %>%
    std_normalize() ->
    parallel_land_temp_df_annual

time_linear_fitting <- function(parallel_temp_df_annual, batch_size){
    grid_cells <- parallel_temp_df_annual$grid_cell %>%as.factor() %>% unique
    parallel_temp_df_annual_batch_n <- filter(parallel_temp_df_annual, grid_cell %in% grid_cells[1:batch_size]) %>% mutate(grid_cell = as.integer(as.factor(grid_cell)))

    assertthat::assert_that(batch_size<= length(grid_cells))

    time_fit <- system.time(
        a_temp_fit <- map2stan(
            alist(
                cell_value ~ dnorm(mu, sig),
                mu <- a[grid_cell] + b[grid_cell]*global_value,
                sig <- sigma[grid_cell],
                a ~ dnorm(0, 5),
                b ~ dnorm(0, 5),
                sigma ~ dcauchy(0,5)
            ),
            data=parallel_temp_df_annual_batch_n %>% as.data.frame(),
            start = list(a=rep(0, batch_size), b=rep(1, batch_size), sigma=rep(1, batch_size)),
            cores=4, chains=4, iter = 3000,
            debug = TRUE)
    )

    rethinking::precis(a_temp_fit, depth = 2)@output %>%
        tibble::rownames_to_column("parameter") %>%
        mutate(n = batch_size, time = time_fit[[3]], grid_cell = rep(as.integer(levels(grid_cells))[1:batch_size], 3)) ->
        annual_temp_fit_batch_n
    return(annual_temp_fit_batch_n)
}

land_size_time_fit <- tibble(parameter = character(), Mean = numeric(),
                        StdDev = numeric(), `lower 0.89`=numeric(),
                        `upper 0.89`=numeric(), n_eff=numeric(), Rhat=numeric(),
                        n = numeric(), time = numeric(), grid_cell = numeric())


## ---- echo=FALSE---------------------------------------------------------
for(i in 1:length(batch_size_list)){
    land_size_time_fit_holder <- time_linear_fitting(parallel_land_temp_df_annual, batch_size_list[[i]])
    land_size_time_fit <- bind_rows(land_size_time_fit, land_size_time_fit_holder)
}


## ------------------------------------------------------------------------
land_size_time_fit %>%
    select(n, time) %>%
    unique() %>%
    mutate(time_per_run = time/n) ->
    land_batch_efficiency


## ---- echo = FALSE-------------------------------------------------------

land_batch_efficiency_plot <- ggplot(land_batch_efficiency, aes(n, time_per_run)) +
    geom_point(size=3) +
    ylab("Run Time per Batch Unit (s)") +
    xlab("Batch Size") +
    ggtitle("Baseline Batch Efficiency") +
    theme(text = element_text(size=20))



## ------------------------------------------------------------------------
#Select which parallel band of cells to look at
# sea_parallel <-  -35.6

sea_temp_df_annual %>%
    # filter(between(lat, sea_parallel-1.5, sea_parallel+1.5)) %>%
    #This Parallel has 114 n=4, 3 n=3, 2 n=2, and 1 n=1 which *just* works with a batch size of 120
    select(group, global_value, cell_value) %>%
    std_normalize() ->
    parallel_sea_temp_df_annual

time_grouped_linear_fitting <- function(parallel_temp_df_annual, batch_size){
    groups <- parallel_temp_df_annual$group %>%as.factor() %>% unique
    parallel_temp_df_annual_batch_n <- filter(parallel_temp_df_annual, group %in% groups[1:batch_size]) %>% mutate(group = as.integer(as.factor(group)))

    assertthat::assert_that(batch_size<= length(unique(groups)))

    time_fit <- system.time(
        a_temp_fit <- map2stan(
            alist(
                cell_value ~ dnorm(mu, sig),
                mu <- a[group] + b[group]*global_value,
                sig <- sigma[group],
                a ~ dnorm(0, 5),
                b ~ dnorm(0, 5),
                sigma ~ dcauchy(0,5)
            ),
            data=parallel_temp_df_annual_batch_n %>% as.data.frame(),
            start = list(a=rep(0, batch_size), b=rep(1, batch_size), sigma=rep(1, batch_size)),
            cores=4, chains=4, iter = 3000,
            debug = TRUE)
    )

    rethinking::precis(a_temp_fit, depth = 2)@output %>%
        tibble::rownames_to_column("parameter") %>%
        mutate(n = batch_size, time = time_fit[[3]], group = rep(as.integer(levels(groups))[1:batch_size], 3)) ->
        annual_temp_fit_batch_n
    return(annual_temp_fit_batch_n)
}

sea_size_time_fit <- tibble(parameter = character(), Mean = numeric(),
                        StdDev = numeric(), `lower 0.89`=numeric(),
                        `upper 0.89`=numeric(), n_eff=numeric(), Rhat=numeric(),
                        n = numeric(), time = numeric(), group = numeric())


## ---- echo=FALSE---------------------------------------------------------
for(i in 1:length(batch_size_list)){
    sea_size_time_fit_holder <- time_grouped_linear_fitting(parallel_sea_temp_df_annual, batch_size_list[[i]])
    sea_size_time_fit <- bind_rows(sea_size_time_fit, sea_size_time_fit_holder)
}



## ------------------------------------------------------------------------
sea_size_time_fit %>%
    select(n, time) %>%
    unique() %>%
    mutate(time_per_run = time/n) ->
    sea_batch_efficiency


## ---- echo = FALSE-------------------------------------------------------

sea_batch_efficiency_plot <- ggplot(sea_batch_efficiency, aes(n, time_per_run)) +
    geom_point(size=3) +
    ylab("Run Time per Batch Unit (s)") +
    xlab("Batch Size") +
    ggtitle("Sea Grouped Batch Efficiency") +
    theme(text = element_text(size=20))

#write outputs
#==================================
write.csv(land_batch_efficiency, file.path("results", "land_batch_efficiency.csv"))
write.csv(sea_batch_efficiency, file.path("results", "sea_batch_efficiency.csv"))
ggsave(file.path("results", "land_batch_efficiency.png"), plot = land_batch_efficiency_plot,  dpi=600/2, width=6000/300, height=3000/300)
ggsave(file.path("results", "sea_batch_efficiency.png"), plot = sea_batch_efficiency_plot,  dpi=600/2, width=6000/300, height=3000/300)

