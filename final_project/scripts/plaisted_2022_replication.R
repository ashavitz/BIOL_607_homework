#'--------------------------
#' Attempt to replicate data analysis from Plaisted et al. 2022:
#' https://doi.org/10.3389/fmars.2022.920699
#' Data download:
#' https://datadryad.org/landing/show?id=doi%3A10.5061%2Fdryad.z34tmpggg#citations
#'--------------------------

# ---- Load Libraries ----
library(dplyr)
# library(broom)
# library(broom.mixed)
library(car) #Companion to Applied Regression - for Anova()
# library(DHARMa)
library(ggplot2)
# library(glmmTMB)
# library(lme4)
library(lubridate)
library(mgcv)
# library(modelbased)
# library(performance)
library(readr)


# Load data from Plaisted et al. 2022
eelgrass_data <- read_csv("data/doi_10_5061_dryad_z34tmpggg__v20220628/seagrass_for_modeling.csv")

eelgrass_data <- eelgrass_data |> 
  mutate(
    # create a presence absence variable
    presence = ifelse(percent_cover > 0, 1, 0),
    
    # convert all instances of 100% cover to 99% cover, and convert to proportion (between 0 and 1)
    percent_cover_prop = ifelse(percent_cover == 100, 0.99, percent_cover / 100),
    
    # create centered temperature variable:
    temp_centered = prior_yr_mean_daily_mean_temp - period_avg_mean_daily_mean_temp,
    
    # created nested ID for site_quadrat for random effects inclusion
    site_quadrat = interaction(location, quadrat, drop = TRUE),
    
    # make location a factor for modeling
    location = as.factor(location)
  )


# Fit the presence model: 
presence_model <- gam(
  presence ~ s(year) + s(temp_centered) + s(location, bs = "re") + s(site_quadrat, bs = "re"),
  data = eelgrass_data,
  family = binomial(link = "logit"),
  method = "REML"
)

# Fit the abundance model: 
abundance_data <- eelgrass_data |>  filter(presence == 1)

abundance_model <- gam(
  percent_cover_prop ~ s(year) + s(temp_centered) + s(location, bs = "re") + s(site_quadrat, bs = "re"),
  data = abundance_data,
  family = betar(link = "logit"),
  method = "REML"
)

# Model diagnostics
library(gratia)

draw(presence_model)
draw(abundance_model)

gam.check(presence_model)
gam.check(abundance_model)


library(visreg)
visreg(presence_model, "temp_centered", scale = "response")
visreg(abundance_model, "temp_centered", scale = "response")

