library(dplyr)
library(tidyr)
library(rethinking)
library(alps)

#Read in Precipitation data (Not currently used)
tempo_aggregated_pr <- read.general("data/pr_Amon_MIROC5_tempo_avg.nc", varname = "pr")
geo_tempo_aggregated_pr <- read.geo.aggregate("data/pr_Amon_MIROC5_geo_tempo_avg.nc", varname = "pr")

#Read in temperature data.
tempo_aggregated_tas <- read.general("data/tas_Amon_MIROC5_tempo_avg.nc", varname = "tas")
geo_tempo_aggregated_tas <- read.geo.aggregate("data/tas_Amon_MIROC5_geo_tempo_avg.nc", varname = "tas")

#Read in monthly temp data
unaggregated_tas <- read.general("data/tas_Amon_MIROC5_cat.nc", varname = "tas")
geo_aggregated_tas <- read.geo.aggregate("data/tas_Amon_MIROC5_geo_avg.nc", varname = "tas")

temp_df_annual <- bind_geo.ave_grid(geo_tempo_aggregated_tas, tempo_aggregated_tas)
temp_df_monthly <- bind_geo.ave_grid(geo_aggregated_tas, unaggregated_tas)

ncell <- 25
temp_df_annual %>%
    filter(grid_cell<ncell) %>%
    std_normalize() ->
    test_frame_annual

# Do we need to normalize before attaching the df?
# The global temperature gets repeated a bunch before normalizing it.
test_frame_monthly <- temp_df_monthly %>% filter(grid_cell<50) %>% std_normalize()
test_frame_january <- test_frame_monthly %>% label_month() %>% filter(month == 1) %>% select(-month)

a_temp_fit <- map2stan(
    alist(
        cell_value ~ dnorm(mu, sig),
        mu <- a[grid_cell] + b[grid_cell]*global_value,
        sig <- sigma[grid_cell],
        a ~ dnorm(0, 5),
        b ~ dnorm(0, 5),
        sigma ~ dcauchy(0,5)#This parameter seems to be causing the failure.
        # As we add more grid cells the slopes become smaller until they are pushed to 0
    ),
    data=test_frame_annual,
    start = list(a=rep(0, ncell), b=rep(1, ncell), sigma=rep(1,ncell)),
    cores=4, chains=4, iter = 2500,   # These parameters are only valid for map2stan
    debug = TRUE)

annual_temp_fit <- precis(a_temp_fit, depth = 2)


m_temp_fit <- map(
    alist(
        cell_value ~ dnorm(mu, sigma),
        mu <- a[grid_cell] + b[grid_cell]*global_value,
        a[grid_cell] ~ dnorm(0, 5),
        b[grid_cell] ~ dnorm(0, 5),
        sigma ~ dnorm(0,10)
    ),
    data=test_frame_january)

monthly_temp_fit <- precis(m_temp_fit, depth = 2)

#Try plotting the data we're trying to fit to see if there are any odd values causing failures
ggplot(test_frame_annual, aes(global_value, cell_value)) +
    geom_point() +
    facet_wrap(.~grid_cell)
