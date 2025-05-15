#'--------------------------
#' Load SST from NOAA CoastWatch ERDDAP
#' Sea Surface Temperature, Multi-scale Ultra-high Resolution (MUR JPL),
#' Daily 1km East Coast EEZ
#' June 2002 - August 2022
#' https://coastwatch.noaa.gov/erddap/info/noaacwecnMURdaily/index.html
#'--------------------------


# ---- Load Libraries ----
library(dplyr)
library(lubridate)
library(purrr) # for map
library(readr)
library(rerddap) # for accessing ERDDAP servers

# ---- Define Functions ----

# Attempting to pull too much data at once from the NOAA CoastWatch ERDDAP server caused failure
# so it was necessary to request subsets of the data. The following function is designed to work
# with the rerddap::griddap function to break up requests into smaller time chunks, and then 
# bind then back together.

# Define a function to download data in chunks (necessary to prevent crashing/errors)
pull_sst_chunks <- function(data_set_info, 
                            fields, 
                            latitude, 
                            longitude, 
                            start_date, 
                            end_date,
                            chunk_months = 6) {
  
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Create sequence of chunk start dates
  chunk_starts <- seq.Date(from = start_date, to = end_date, by = paste0(chunk_months, " months"))
  
  # Make chunk end dates (6 months later or up to end_date)
  # Note - %m+% lets you add a month to a date without exceeding the last day of the new month
  # https://search.r-project.org/CRAN/refmans/lubridate/html/mplus.html
  # Note - this line basically take the starting date, adds 6 months, and then for each section
  # (pmin for each combo instead of just 1 min) takes that 6 month later value, or the set end_date
  # if end_date is earlier that the calculated end date
  chunk_ends <- pmin(chunk_starts %m+% months(chunk_months) - days(1), end_date)
  
  # Pull data for each chunk
  sst_list <- map2(chunk_starts, chunk_ends, function(chunk_start, chunk_end) {
    message(paste("Pulling chunk:", chunk_start, "to", chunk_end))
    
    # Much of the structure below using tryCatch is inspired by:
    # https://stackoverflow.com/questions/12193779/how-to-use-the-trycatch-function
    tryCatch({
      griddap(data_set_info,
              fields = fields,
              latitude = latitude,
              longitude = longitude,
              time = c(as.character(chunk_start), as.character(chunk_end))
      )$data |> 
        
        # Above returns some columns as arrays - need to flatten
        mutate(across(where(is.array), as.vector))
    }, error = function(e) {
      message(paste("Failed chunk:", chunk_start, "-", chunk_end))
      return(NULL)
    })
  })
  
  sst_combined <- bind_rows(sst_list)
  
  return(sst_combined)
}


# ---- Load Data ----
# NOTE - The earliest available day is 2002-06-02

# Get information about the data set
data_set_info <- rerddap::info("noaacwecnMURdaily",
                               url = "coastwatch.noaa.gov/erddap")


# Duck Harbor, Cqpe Cod, MA, 2002-06-02 through 2022-08-28
sst_data_dh <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(41.94, 41.94),
  longitude = c(-70.09, -70.09),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)

# Examine data
str(sst_data_dh)
summary(sst_data_dh)
visdat::vis_dat(sst_data_dh)
View(sst_data_dh)

# Export data to locally saved csv (to avoid having to do pull from ERDDAP each time)
# write_csv(sst_data_dh, "data/sst_data_dh.csv")

# # One at a time
# # Duck Harbor
# sst_data_<- griddap(data_set_info,
#                     fields = c("time", "latitude", "longitude", "sst"),
#                     latitude = c(41.93,41.95),
#                     longitude = c(-70.09, -70.08),
#                     time = c("2002-06-02", "2002-12-31")
#                     )$data



# Tinges Island, MD, 2002-06-02 through 2022-08-28
sst_data_ti <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(38.18, 38.18),
  longitude = c(-75.18, -75.18),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_ti, "data/sst_data_ti.csv")


# Fire Island, NY, 2002-06-02 through 2022-08-28
sst_data_fi <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(40.72, 40.72),
  longitude = c(-72.94, -72.94),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 3
)
# write_csv(sst_data_fi, "data/sst_data_fi.csv")


# Fort Getty, Narragansett Bay, RI, 2002-06-02 through 2022-08-28 
sst_data_fg <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(41.51, 41.51),
  longitude = c(-71.39, -71.39),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_fg, "data/sst_data_fg.csv")


# Prudence Island (T-Dock), Narragansett Bay, RI, 2002-06-02 through 2022-08-28 
sst_data_pi <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(41.58, 41.58),
  longitude = c(-71.32, -71.32),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_pi, "data/sst_data_pi.csv")


# Pleasant Bay, Cape Cod, MA, 2002-06-02 through 2022-08-28 
sst_data_pb <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(41.75, 41.75),
  longitude = c(-69.96, -69.96),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_pb, "data/sst_data_pb.csv")


# West Beach, Beverly, MA, 2002-06-02 through 2022-08-28 
sst_data_wb <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(42.55, 42.55),
  longitude = c(-70.81, -70.81),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_wb, "data/sst_data_wb.csv")


# Greay Bat, NH, 2002-06-02 through 2022-08-28
# The Sea Surface Temperature, Multi-scale Ultra-high Resolution (MUR JPL), Daily 1km East Coast EEZ
# data set does not have SST for Great Bay (GB is a tidal estuary). SST from the closest point,
# the mouth of the Piscataqua River, is used as a proxy (may be excluded in analysis).  
sst_data_gb_proxy <- pull_sst_chunks(
  data_set_info = data_set_info,
  fields = c("time", "latitude", "longitude", "sst"),
  latitude = c(43.08, 43.08),
  longitude = c(-70.71, -70.71),
  start_date = "2002-06-02",
  end_date = "2022-08-28",
  chunk_months = 6
)
# write_csv(sst_data_gb_proxy, "data/sst_data_gb_proxy.csv")


##################################
# Process for collecting SST data
# Determine specific lat/long coordinates for site, at hundredth decimal places
# Determine full time range (earlier initial sampling date to latest last sampling date)
# - Eelgrass cover data public range: 2003 - 2019 on seagrassnet. 
#   2003-2021 is available from Plaisted 2022 paper. Would need 2002-06-01 through 2022-08-31
#   BUT sst temp data from the dataset I'm using is only through 2022-08-28. That's OK.
# 
# NOT NEEDED, DATA IS COMPLETE
# Pull data for each the specific site by lat long - if there are data gaps
# Two step data interpolation: 
# 1) Use nearest neighbor OR average of all values in a small radius. Likely radius (e.g. 3-5km)
#    because nearest neighbor can be an issue for coastlines (source?)
# 2) Where there are still missing values, use simple linear interpolation* - Make sure there
#    aren't long time gaps when this method is used. 








