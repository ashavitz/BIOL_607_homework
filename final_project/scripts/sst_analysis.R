#'--------------------------
#' SST impact on eelgrass cover - Analysis
#' SST data from NOAA CoastWatch
#' Eelgrass cover data from seagrassnet.org
#'--------------------------

# ---- Load Libraries ----
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)


# ---- Load and Organize SST Data ----
# ---
# SST data previously downloaded from NOAA CoastWatch ERDDAP
# Subsetted data stored on local drive for speed of access locally
# Original SST data source:
# Sea Surface Temperature, Multi-scale Ultra-high Resolution (MUR JPL),
# Daily 1km East Coast EEZ
# June 2002 - August 2022
# https://coastwatch.noaa.gov/erddap/info/noaacwecnMURdaily/index.html
# ---

# Load each SST df, adding column for site name
sst_data_dh <- read_csv("data/sst_data_dh.csv") |> mutate(site = "duck_harbor")
sst_data_fg <- read_csv("data/sst_data_fg.csv") |> mutate(site = "fort_getty")
sst_data_fi <- read_csv("data/sst_data_fi.csv") |> mutate(site = "fire_island")
sst_data_gb_proxy <- read_csv("data/sst_data_gb_proxy.csv") |> mutate(site = "great_bay_proxy")
sst_data_pb <- read_csv("data/sst_data_pb.csv") |> mutate(site = "pleasant_bay")
sst_data_pi <- read_csv("data/sst_data_pi.csv") |> mutate(site = "prudence_island")
sst_data_ti <- read_csv("data/sst_data_ti.csv") |> mutate(site = "tingles_island")
sst_data_wb <- read_csv("data/sst_data_wb.csv") |> mutate(site = "west_beach")

# Bind data vertically
sst_data <-
  bind_rows(
    sst_data_dh,
    sst_data_fg, 
    sst_data_fi, 
    sst_data_gb_proxy, 
    sst_data_pb, 
    sst_data_pi, 
    sst_data_ti, 
    sst_data_wb
  ) |> 
  # Remove unnecessary columns (long, lat, level)
  select(-c(latitude, longitude, level)) |> 
  # Add columns for year and month
  mutate(year = year(time),
         month = month(time))

# Remove unnecessary dfs from environment
remove(list = c("sst_data_dh", "sst_data_fg", "sst_data_fi", "sst_data_gb_proxy",
                "sst_data_pb", "sst_data_pi", "sst_data_ti", "sst_data_wb"))

# Viz all data
ggplot(sst_data, aes(x = time, y = sst)) +
  geom_point() +
  facet_wrap(~site)

# Calculate summer mean values (June - August)
sst_summer_means <- sst_data |> 
  filter(month %in% c(6,7,8)) |> 
  group_by(year, site) |> 
  summarize(mean_sst = mean(sst), .groups = "drop")

# Calculate average summer mean by site across entire time period
sst_summer_average_means <- sst_data |> 
  filter(month %in% c(6,7,8)) |> 
  group_by(site) |> 
  summarize(average_mean_sst = mean(sst), .groups = "drop")

# Create df with centered mean summer temperatures
# (Calculate difference between each site summer mean and site average summer mean)
sst_summer <- sst_summer_means |> 
  left_join(sst_summer_average_means, by = c("site")) |> 
  mutate(sst_centered = mean_sst - average_mean_sst)

# Visualize
ggplot(sst_summer,
       aes(x = year, y = sst_centered, color = site)) +
  geom_line()


# ---- Load and Organize SeagrassNet Data ----

# Load each SeagrassNet df, adding column for site name
dh_data <- read_csv("data/MA20.1.csv") |> mutate(site = "duck_harbor")
fg_data <- read_csv("data/RN31.1.csv") |> mutate(site = "fort_getty")
fi_data <- read_csv("data/NY41.1.csv") |> mutate(site = "fire_island")
gb_data <- read_csv("data/NH9.2.csv") |> mutate(site = "great_bay_proxy")
pb_data <- read_csv("data/MA20.2.csv") |> mutate(site = "pleasant_bay")
pi_data <- read_csv("data/RN31.2.csv") |> mutate(site = "prudence_island")
ti_data <- read_csv("data/MD12.2.csv") |> mutate(site = "tingles_island")
wb_data <- read_csv("data/MS45.1.csv") |> mutate(site = "west_beach")


# Bind data vertically
eelgrass_data <-
  bind_rows(
    dh_data,
    fg_data, 
    fi_data, 
    gb_data, 
    pb_data, 
    pi_data, 
    ti_data,
    wb_data
  ) |> 
  # Only include Zostera marina data, and only include shallow transect (transect A)
  filter(scientific_name == "Zostera marina", transect_name == "A") |>
  mutate(year = year(observation_date)) |> 
  # Reorganize and remove unnecessary variables
  select(c(year, site, transect_name, quadrat_number, observation_date, percent_cover)) |> 
  relocate(site) |> 
  # Remove duplicated rows
  distinct()
  
# Remove unnecessary dfs from environment
remove(list = c("dh_data", "fg_data", "fi_data", "gb_data",
                "pb_data", "pi_data", "ti_data", "wb_data"))


# I am not interested in seasonal variation, only annual, and my predictor is annual (temp),
# so aggregate (average) across quadrats within year site each year
transect_a_aggregated <- eelgrass_data |> 
  group_by(site, year, quadrat_number) |> 
  summarize(percent_cover = mean(percent_cover), .groups = "drop") |> 
  # average at the site year level |> 
  group_by(site, year) |> 
  summarize(percent_cover = mean(percent_cover), .groups = "drop") |> 
  # Join sst and transect A aggregated percent cover data
  left_join(sst_summer, by = c("site", "year"))



