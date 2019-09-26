library(dplyr)
library(tidyr)
library(rethinking)
library(alps)



# #Read in Precipitation data (Not currently used)
# tempo_aggregated_pr <- read.general("data/pr_Amon_MIROC5_tempo_avg.nc", varname = "pr")
# geo_tempo_aggregated_pr <- read.geo.aggregate("data/pr_Amon_MIROC5_geo_tempo_avg.nc", varname = "pr")

#Read in temperature data.
tempo_aggregated_tas <- read.general("data/tas_Amon_MIROC5_tempo_avg.nc", varname = "tas")
geo_tempo_aggregated_tas <- read.geo.aggregate("data/tas_Amon_MIROC5_geo_tempo_avg.nc", varname = "tas")

# #Read in monthly temp data
# unaggregated_tas <- read.general("data/tas_Amon_MIROC5_cat.nc", varname = "tas")
# geo_aggregated_tas <- read.geo.aggregate("data/tas_Amon_MIROC5_geo_avg.nc", varname = "tas")

long_tempo_aggregated_tas <- gather_griddata(tempo_aggregated_tas)
long_geo_tempo_aggregated_tas <- data.frame(global_value = geo_tempo_aggregated_tas$vardata %>% rep(length(geo_tempo_aggregated_tas$time)))

long_tempo_aggregated_tas %>%
    mutate(global_value = geo_tempo_aggregated_tas$vardata %>%
                            rep(length(unique(long_tempo_aggregated_tas$grid_cell)))) %>%
    rename(cell_value = value) ->
    temp_df_annual
# temp_df_monthly <- bind_geo.ave_grid(geo_aggregated_tas, unaggregated_tas)

# ncell <- 2
# temp_df_annual %>%
#     filter(grid_cell<ncell) %>%
#     std_normalize() ->
#     test_frame_annual
#
# # Do we need to normalize before attaching the df?
# # The global temperature gets repeated a bunch before normalizing it.
# # test_frame_monthly <- temp_df_monthly %>% filter(grid_cell<ncell) %>% std_normalize()
# # test_frame_january <- test_frame_monthly %>% label_month() %>% filter(month == 1) %>% select(-month)
#
# a_temp_fit <- map2stan(
#     alist(
#         cell_value ~ dnorm(mu, sig),
#         mu <- a[grid_cell] + b[grid_cell]*global_value,
#         sig <- sigma[grid_cell],
#         a ~ dnorm(0, 5),
#         b ~ dnorm(0, 5),
#         sigma ~ dcauchy(0,5)
#     ),
#     data=test_frame_annual,
#     start = list(a=rep(0, ncell), b=rep(1, ncell), sigma=rep(1,ncell)),
#     cores=4, chains=4, iter = 2500,   # These parameters are only valid for map2stan
#     debug = TRUE)
#
# annual_temp_fit <- precis(a_temp_fit, depth = 2)


# m_temp_fit <- map(
#     alist(
#         cell_value ~ dnorm(mu, sigma),
#         mu <- a[grid_cell] + b[grid_cell]*global_value,
#         a[grid_cell] ~ dnorm(0, 5),
#         b[grid_cell] ~ dnorm(0, 5),
#         sigma ~ dnorm(0,10)
#     ),
#     data=test_frame_january)
#
# monthly_temp_fit <- precis(m_temp_fit, depth = 2)
#
# #Try plotting the data we're trying to fit to see if there are any odd values causing failures
# ggplot(test_frame_annual, aes(global_value, cell_value)) +
#     geom_point() +
#     facet_wrap(.~grid_cell)

#time how long running the model takes
#======================================================================================
parallel <-  0.7003838

temp_df_annual %>%
    filter(between(lat, parallel-.1, parallel+.1)) %>%
    std_normalize() ->
    parallel_temp_df_annual

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
            data=parallel_temp_df_annual_batch_n,
            start = list(a=rep(0, batch_size), b=rep(1, batch_size), sigma=rep(1, batch_size)),
            cores=4, chains=4, iter = 2500,   # These parameters are only valid for map2stan
            debug = TRUE)
    )

    rethinking::precis(a_temp_fit, depth = 2)@output %>%
        tibble::rownames_to_column("parameter") %>%
        mutate(n = batch_size, time = time_fit[[3]], grid_cell = rep(as.integer(levels(grid_cells))[1:batch_size], 3)) ->
    annual_temp_fit_batch_n
    return(annual_temp_fit_batch_n)
}
time_linear_fitting_strict <- function(batch_size, parallel_temp_df_annual, parameter_df){
    #Get a df of gridcells at the same lat
    parallel_temp_df_annual_batch_n <- filter(parallel_temp_df_annual, grid_cell %in% parallel_grid_cells[1:batch_size]) %>% mutate(grid_cell = as.factor(grid_cell))
    factor_grid_cell_mapping <- levels(parallel_temp_df_annual_batch_n$grid_cell)
    parallel_temp_df_annual_batch_n <- mutate(parallel_temp_df_annual_batch_n, grid_cell = as.integer(grid_cell))

    a_mean <- parameter_df %>% filter(n == batch_size, grepl("^a\\[", parameter)) %>% .$Mean %>% mean()
    a_std <- parameter_df %>% filter(n == batch_size, grepl("^a\\[", parameter)) %>% .$StdDev %>% mean()
    b_mean <- parameter_df %>% filter(n == batch_size, grepl("^b\\[", parameter)) %>% .$Mean %>% mean()
    b_std <- parameter_df %>% filter(n == batch_size, grepl("^b\\[", parameter)) %>% .$StdDev %>% mean()
    sigma_mean <- parameter_df %>% filter(n == batch_size, grepl("^sigma\\[", parameter)) %>% .$Mean %>% mean()
    sigma_std <- parameter_df %>% filter(n == batch_size, grepl("^sigma\\[", parameter)) %>% .$StdDev %>% mean()


    time_fit <- system.time(
        a_temp_fit <- map2stan(
            alist(
                cell_value ~ dnorm(mu, sig),
                mu <- a[grid_cell] + b[grid_cell]*global_value,
                sig <- sigma[grid_cell],
                a ~ dnorm(0.220, 0.002205),
                b ~ dnorm(0.321, 0.002205),
                sigma ~ dcauchy(0.175, 0.001565)
            ),
            data=parallel_temp_df_annual_batch_n,
            start = list(a=rep(0, batch_size), b=rep(1, batch_size), sigma=rep(1, batch_size)),
            cores=4, chains=4, iter = 2500,   # These parameters are only valid for map2stan
            debug = TRUE)
    )

    rethinking::precis(a_temp_fit, depth = 2)@output %>%
        tibble::rownames_to_column("parameter") %>%
        mutate(n = batch_size, time = time_fit[[3]], grid_cell = rep(parallel_grid_cells[1:batch_size], 3)) ->
        annual_temp_fit_batch_n
    return(annual_temp_fit_batch_n)
}
time_linear_fitting_mid <- function(batch_size, parallel_temp_df_annual, parameter_df){
    #Get a df of gridcells at the same lat
    parallel_temp_df_annual_batch_n <- filter(parallel_temp_df_annual, grid_cell %in% parallel_grid_cells[1:batch_size]) %>% mutate(grid_cell = as.factor(grid_cell))
    factor_grid_cell_mapping <- levels(parallel_temp_df_annual_batch_n$grid_cell)
    parallel_temp_df_annual_batch_n <- mutate(parallel_temp_df_annual_batch_n, grid_cell = as.integer(grid_cell))

    a_mean <- parameter_df %>% filter(n == batch_size, grepl("^a\\[", parameter)) %>% .$Mean %>% mean()
    a_std <- parameter_df %>% filter(n == batch_size, grepl("^a\\[", parameter)) %>% .$StdDev %>% mean()
    b_mean <- parameter_df %>% filter(n == batch_size, grepl("^b\\[", parameter)) %>% .$Mean %>% mean()
    b_std <- parameter_df %>% filter(n == batch_size, grepl("^b\\[", parameter)) %>% .$StdDev %>% mean()
    sigma_mean <- parameter_df %>% filter(n == batch_size, grepl("^sigma\\[", parameter)) %>% .$Mean %>% mean()
    sigma_std <- parameter_df %>% filter(n == batch_size, grepl("^sigma\\[", parameter)) %>% .$StdDev %>% mean()


    time_fit <- system.time(
        a_temp_fit <- map2stan(
            alist(
                cell_value ~ dnorm(mu, sig),
                mu <- a[grid_cell] + b[grid_cell]*global_value,
                sig <- sigma[grid_cell],
                a ~ dnorm(0.220, 0.005),
                b ~ dnorm(0.321, 0.005),
                sigma ~ dcauchy(0.175, 0.0032)
            ),
            data=parallel_temp_df_annual_batch_n,
            start = list(a=rep(0, batch_size), b=rep(1, batch_size), sigma=rep(1, batch_size)),
            cores=4, chains=4, iter = 2500,   # These parameters are only valid for map2stan
            debug = TRUE)
    )

    rethinking::precis(a_temp_fit, depth = 2)@output %>%
        tibble::rownames_to_column("parameter") %>%
        mutate(n = batch_size, time = time_fit[[3]], grid_cell = rep(parallel_grid_cells[1:batch_size], 3)) ->
        annual_temp_fit_batch_n
    return(annual_temp_fit_batch_n)
}

batch_size_list <- c(2, 3, 4, 5)#, 6, 8, 10, 12, 15, 20, 24, 30, 40, 60, 120)

#Baseline
#===================================================================================
size_time_fit <- tibble(parameter = character(), Mean = numeric(),
                        StdDev = numeric(), `lower 0.89`=numeric(),
                        `upper 0.89`=numeric(), n_eff=numeric(), Rhat=numeric(),
                        n = numeric(), time = numeric(), grid_cell = numeric())
for(i in 1:length(batch_size_list)){
    size_time_fit_holder <- time_linear_fitting(parallel_temp_df_annual, batch_size_list[[i]])
    size_time_fit <- bind_rows(size_time_fit, size_time_fit_holder)
}
#===================================================================================

#Strict
#===================================================================================
size_time_fit_strict <- tibble(parameter = character(), Mean = numeric(),
                        StdDev = numeric(), `lower 0.89`=numeric(),
                        `upper 0.89`=numeric(), n_eff=numeric(), Rhat=numeric(),
                        n = numeric(), time = numeric(), grid_cell = numeric())
for(i in 1:length(batch_size_list)){
    size_time_fit_strict_holder <- time_linear_fitting_strict(batch_size_list[[i]], parallel_temp_df_annual, size_time_fit)
    size_time_fit_strict <- bind_rows(size_time_fit_strict, size_time_fit_strict_holder)
}
#===================================================================================

#Mid
#===================================================================================
size_time_fit_mid <- tibble(parameter = character(), Mean = numeric(),
                               StdDev = numeric(), `lower 0.89`=numeric(),
                               `upper 0.89`=numeric(), n_eff=numeric(), Rhat=numeric(),
                               n = numeric(), time = numeric(), grid_cell = numeric())
for(i in 1:length(batch_size_list)){
    size_time_fit_mid_holder <- time_linear_fitting_mid(batch_size_list[[i]], parallel_temp_df_annual, size_time_fit)
    size_time_fit_mid <- bind_rows(size_time_fit_strict, size_time_fit_mid_holder)
}
#===================================================================================

#Calculate Efficiency
#===================================================================================
size_time_fit %>%
    select(n, time) %>%
    unique() %>%
    mutate(time_per_run = time/n) ->
    batch_efficiency

size_time_fit_strict %>%
    select(n, time) %>%
    unique() %>%
    mutate(time_per_run = time/n) ->
    batch_efficiency_strict

size_time_fit_mid %>%
    select(n, time) %>%
    unique() %>%
    mutate(time_per_run = time/n) ->
    batch_efficiency_mid
#===================================================================================


#Plotting
#===================================================================================

ggplot(batch_efficiency, aes(n, time_per_run)) +
    geom_point(size=3) +
    ylab("Run Time per Batch Unit (s)") +
    xlab("Batch Size") +
    ylim(10, 40) +
    ggtitle("Baseline Batch Efficiency") +
    theme(text = element_text(size=20)) ->
    baseline_batch_efficiency_plot

ggplot(batch_efficiency_strict, aes(n, time_per_run)) +
    geom_point(size=3) +
    ylab("Run Time per Batch Unit (s)") +
    xlab("Batch Size") +
    ylim(10, 40) +
    ggtitle("Stricter Batch Efficiency") +
    theme(text = element_text(size=20)) ->
    strict_batch_efficiency_plot

ggsave("scratch_running/batch_efficiency.png", plot = baseline_batch_efficiency_plot,  dpi=600/2, width=6000/300, height=3000/300)
ggsave("scratch_running/batch_efficiency_strict.png",  plot = strict_batch_efficiency_plot, dpi=600/2, width=6000/300, height=3000/300)

#===================================================================================
