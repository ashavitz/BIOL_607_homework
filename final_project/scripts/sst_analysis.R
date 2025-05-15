# ---- Header ----
#'--------------------------
#' SST impact on eelgrass cover - Analysis
#' SST data from NOAA CoastWatch
#' Eelgrass cover data from Plaisted et al. 2022: https://doi.org/10.3389/fmars.2022.920699
#' Data download:
#' https://datadryad.org/landing/show?id=doi%3A10.5061%2Fdryad.z34tmpggg#citations
#' Partial data is available at seagrassnet.org
#'--------------------------

# ---- Load Libraries ----
library(dplyr)
library(broom)
library(broom.mixed)
library(car) # Companion to Applied Regression - for Anova()
library(DHARMa)
library(emmeans)
library(ggplot2)
library(glmmTMB)
library(gratia) # For examining GAMs made with mgcv
library(lme4)
library(lubridate)
library(mgcv) # For GAMs and thin plate regression spline smooths in modeling, per Plaisted 2022
library(modelbased)
library(performance)
library(readr)
library(tidyr)
library(visreg)


# ---- Set Global ggplot Themes ----

# Set ggplot theme to minimal, rotate x axis labels, and center plot titles
theme_set(
  theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)


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
sst_data_dh <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_dh.csv")) |> mutate(site = "Duck Harbor Beach")
sst_data_fg <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_fg.csv")) |> mutate(site = "Fort Getty")
sst_data_fi <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_fi.csv")) |> mutate(site = "Fire Island")
sst_data_gb_proxy <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_gb_proxy.csv")) |> mutate(site = "Great Bay")
sst_data_pb <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_pb.csv")) |> mutate(site = "Pleasant Bay")
sst_data_pi <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_pi.csv")) |> mutate(site = "T-Dock")
sst_data_ti <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_ti.csv")) |> mutate(site = "Tingles Island")
sst_data_wb <- read_csv(here::here("data/noaa_coastwatch_sst/sst_data_wb.csv")) |> mutate(site = "Salem Sound")

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

# # Viz all data
# ggplot(sst_data, aes(x = time, y = sst)) +
#   geom_point() +
#   facet_wrap(~site)

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
  mutate(sst_centered = mean_sst - average_mean_sst) |> 
  # create column for previous summer's centered sst
  arrange(site, year) |> 
  group_by(site) |> 
  mutate(
    mean_sst_prev = lag(mean_sst),
    sst_centered_prev = lag(sst_centered)) |> 
  ungroup()

# Visualize
ggplot(sst_summer,
       aes(x = year, y = sst_centered, color = site)) +
  geom_line()


# ---- Load and Organize Data from Plaisted et al. 2022 ----
#'--------------------------
#' Eelgrass cover data from Plaisted et al. 2022: https://doi.org/10.3389/fmars.2022.920699
#' Data download:
#' https://datadryad.org/landing/show?id=doi%3A10.5061%2Fdryad.z34tmpggg#citations
#'--------------------------

# NOTE - This data only includes the shallowest transects, to standardize depth.
# This is transect "A" in all cases, except for Fire Island it is "B"
eelgrass_data <- read_csv(here::here("data/doi_10_5061_dryad_z34tmpggg__v20220628/seagrass_for_modeling.csv")) |> 
  rename(site = location)

eelgrass_data <- eelgrass_data |> 
  mutate(
    # create centered temperature variable
    temp_centered_prev = prior_yr_mean_daily_mean_temp - period_avg_mean_daily_mean_temp
  )


# ---- Compare CoastWatch sst to Plaisted et al. 2022 Temperature Data ----

# Create a joined data frame
temp_compare <- left_join(eelgrass_data, sst_summer, by = c("site", "year")) |> 
  select(c(year, site, prior_yr_mean_daily_mean_temp, 
           mean_sst_prev, sst_centered_prev, temp_centered_prev)) |> 
  distinct() |> 
  rename(plaisted_mean_summer_sst_prev = prior_yr_mean_daily_mean_temp) |> 
  rename(coastwatch_mean_summer_sst_prev = mean_sst_prev) |> 
  rename(plaisted_centered_sst_prev = temp_centered_prev) |> 
  rename(coastwatch_centered_sst_prev = sst_centered_prev)

# Pivot longer for visualization
temp_compare_long <- temp_compare |> 
  select(-c(plaisted_centered_sst_prev, coastwatch_centered_sst_prev)) |> 
  pivot_longer(
    cols = c(plaisted_mean_summer_sst_prev, coastwatch_mean_summer_sst_prev),
    names_to = "source",
    values_to = "mean_summer_sst_prev"
  )

centered_temp_compare_long <- temp_compare |> 
  select(-c(plaisted_mean_summer_sst_prev, coastwatch_mean_summer_sst_prev)) |> 
  pivot_longer(
    cols = c(plaisted_centered_sst_prev, coastwatch_centered_sst_prev),
    names_to = "source",
    values_to = "centered_mean_summer_sst_prev"
  )

# ---- Plot Compare CoastWatch SST to Plaisted Temperature ----
# Previous Summer mean sst
ggplot(temp_compare_long, aes(x = year, y = mean_summer_sst_prev, color = source)) +
  geom_line() +
  facet_wrap(~site) +
  labs(x = "Year",
       y = "Previous Summer Mean SST (°C)",
       color = "Data Source"
       ) + 
  scale_color_discrete(
    labels = c("coastwatch_mean_summer_sst_prev" = "NOAA Coastwatch JPL",
               "plaisted_mean_summer_sst_prev" = "Plaisted at al. 2022"))

# ---- Plot Compare Anomalies CoastWatch SST to Plaisted Temperature ----
# Previous summer centered mean sst
ggplot(centered_temp_compare_long, aes(x = year, y = centered_mean_summer_sst_prev, color = source)) +
  geom_line() +
  facet_wrap(~site) +
  labs(x = "Year",
       y = "Centered Previous Summer Mean SST (°C)",
       color = "Data Source"
  ) + 
  scale_color_discrete(
    labels = c("coastwatch_centered_sst_prev" = "NOAA Coastwatch JPL",
               "plaisted_centered_sst_prev" = "Plaisted at al. 2022"))

# ---- Remove Salem Sound Data from Analysis ----
# After observing that the satellite SSTs were consistently lower than in situ measurements
# except at Salem Sound, I reported this discrepency to Mass Bays, who informed me that the 
# shallowest transect depth at Salem South is significantly deeper than at other sites. 
# Thus, I feel it is appropriate to remove Salem Sound from the analysis...

# ---- Correlation Between CoastWatch SST to Plaisted Temperature ----
# Build simple linear model to examine correlation between 
# CoastWatch sst and Plaisted at el. 2022 temps
temp_compare_lm <- lm(plaisted_mean_summer_sst_prev ~ coastwatch_mean_summer_sst_prev, 
                      data = temp_compare)

summary(temp_compare_lm)
r2(temp_compare_lm)
visreg(temp_compare_lm, "coastwatch_mean_summer_sst_prev")

# Include site as a random effect
temp_compare_glmm <-
  glmmTMB(plaisted_mean_summer_sst_prev ~ coastwatch_mean_summer_sst_prev + (1 | site),
          data = temp_compare)

summary(temp_compare_glmm)
r2(temp_compare_glmm)
# The conditional R2 is near 1, so the entire model explains the variance well,
# but the conditional R2 is very small, so CoastWatch mean summer sst doesn't 
# explain the variance of the mean summer sst values used in the Plaisted paper well.
# This suggests to me that we may get different results from our model compared to theirs.
visreg(temp_compare_glmm, "coastwatch_mean_summer_sst_prev")


# ---- Correlation Between Centered Temperatures ----
# Build simple linear model to examine correlation between centered
# CoastWatch sst and centered Plaisted at el. 2022 temps
centered_temp_compare_lm <- lm(plaisted_centered_sst_prev ~ coastwatch_centered_sst_prev, 
                               data = temp_compare)

summary(centered_temp_compare_lm)
r2(centered_temp_compare_lm)
visreg(centered_temp_compare_lm, "coastwatch_centered_sst_prev")


# ---- Prepare Data for Causal Modeling ----

# Reorganize and remove unnecessary variables
eelgrass_data <- eelgrass_data |> 
  select(c(year, site, quadrat, percent_cover)) |> 
  relocate(site)

## Prepare data for Aggregated Percent Cover Model
# Aggregate (average) at the site level in order to avoid requiring inclusion of 
# quadrat random effects, and to reduce 0 and 1 values for binomial regression.
# Also, My predictor is annual site level (temp), so my outcome (percent cover) will be the same.
eelgrass_data_aggregated <- eelgrass_data |> 
  group_by(site, year) |> 
  summarize(percent_cover = mean(percent_cover) / 100, .groups = "drop") |> 
  
  # Join sst and aggregated percent cover data
  left_join(sst_summer, by = c("site", "year")) |> 

  # If the average percent cover is 100%, convert to 99% because there is no meaningful ecological
  # difference, and it allows beta regression to work
  mutate(percent_cover = ifelse(percent_cover == 1, 0.99, percent_cover))


## Prepare data for Presence Absence Model
# For the presence absence model, I will not aggregate at the site level, and will include 
# quadrat level random effects

# Join sst data to eelgrass_data
eelgrass_data <- eelgrass_data |> 
  left_join(sst_summer, by = c("site", "year")) |> 
  mutate(
    # create column for presence absence (1 = presence, 0 = absence)
    presence = ifelse(percent_cover == 0, 0, 1),
    # create nested ID for site_quadrat for random effects inclusion
    site_quadrat = paste0(site, "_", quadrat)
    ) |> 
  relocate(presence, .after = percent_cover) 


# ## Prepare data for Non-aggregated Percent Cover Model
# # Convert 1 to 0.99
# # 0% cover values are removed on the premise that this model is given coverage
# # Convert from 0-100 scale to 0-1, for modeling
# eelgrass_data_non_ag <- eelgrass_data |> 
#   mutate(
#     percent_cover = ifelse(percent_cover == 100, 99, percent_cover) / 100
#   ) |> 
#   # filter out zero percent cover values
#   filter(percent_cover != 0)


# ---- Visualize Data Prior to Modeling ----

## Visualization for Aggregated Percent Cover Model
ggplot(eelgrass_data_aggregated,
       aes(x = sst_centered_prev, y = percent_cover)) +
  geom_point() +
  facet_wrap(~site)

## Visualization for Presence-Absence Model

# Create df with centered previous sst and propotion of plots where eelgrass is present
props_present <- eelgrass_data |> 
  group_by(sst_centered_prev) |> 
  summarize(prop_present = sum(presence) / n(), .groups = "drop")

ggplot(props_present,
       aes(x = sst_centered_prev, y = prop_present)) +
  geom_point()

# ## Visualization for Non-aggregated Percent Cover Model
# ## (Quadrat Percent Cover *Given Coverage* Model)
# ggplot(eelgrass_data_non_ag,
#        aes(x = sst_centered_prev, y = percent_cover)) +
#   geom_point() +
#   facet_wrap(~site)


# ---- Aggregated Percent Cover Model ----

# Build Model to predict aggregated percent cover:
# - Outcome is percent cover 
# - Predictor is the centered mean sst from the previous year
# - Beta error distribution because % cover data
# - Logit link because working with percentages bounded between 0 & 1
# - There are zero values, need to account for zero inflation
# - Random effects include year and site, as these both will cause group-level variability
eelgrass_pc_glmm <- glmmTMB(
  percent_cover ~ sst_centered_prev + (1 | year) + (1 + sst_centered_prev | site),
  family = beta_family(link = "logit"),
  ziformula = ~1,
  data = eelgrass_data_aggregated)

# # Include year as predictor version
# eelgrass_pc_glmm <- glmmTMB(
#   percent_cover ~ sst_centered_prev + year + (1 + sst_centered_prev | site),
#   family = beta_family(link = "logit"),
#   ziformula = ~1,
#   data = eelgrass_data_aggregated)

# # Mundlak version
# eelgrass_pc_glmm <- glmmTMB(
#   percent_cover ~ sst_centered_prev + mean_sst_prev + (1 + sst_centered_prev | site),
#   family = beta_family(link = "logit"),
#   ziformula = ~1,
#   data = eelgrass_data_aggregated)

# eelgrass_pc_glmm_spline <- glmmTMB(
#   percent_cover ~ splines::ns(sst_centered_prev, df = 2) + (1 | year) + (1 + sst_centered_prev | site),
#   family = beta_family(link = "logit"),
#   ziformula = ~1,
#   data = eelgrass_data_aggregated)


# AIC(eelgrass_pc_glmm, eelgrass_pc_glmm_spline)


# ---- Check Assumptions - Aggregated Percent Cover Model ----
check_model(eelgrass_pc_glmm)

# Validity - Do X and Y Reflect Concepts I'm interested In
# - Conceptual check
#   - Yes, we are interested in water temperature, and sst anomalies can serve as a proxy for
#     bottom water anomalies in shallow areas where eelgrass grows.

# Representativeness: Does Your Data Represent the Population?
# - Conceptual check 
#   - Our sample we are working with only represents the population of the most shallow growing 
#     eelgrass on the Eastern US. Variety of different mean temps.

# Model captures features in the data
# - Posterior predictive check - Looks good!

# Linearity
# - Assuming linearity between the predictor and the log-odds (for binomial) 
#   or logit-transformed mean (for beta regression)
# - Simulate residual
simulateResiduals(eelgrass_pc_glmm, plot = TRUE)
# - No issues identified

# Additivity
# - Assumption that fixed effects impact is additive - they don't affect each other
#   - All set here, because I only have one fixed effect

# Independence of errors
# - Are all replicates truly independent
#   - Look at fitted values vs residuals
simulateResiduals(eelgrass_pc_glmm, plot = TRUE)
testResiduals(eelgrass_pc_glmm)
#   - This suggest OK
#   - However, logically samples will be correlated year to year, within site, 
#     and between neighboring quadrats. Aggregating at the site level to deal with quadrat
#     autocorrelation takes care of that. 
#   - Inclusion of random effects for year and site attempts to handle how different years 
#     and sites behave, but I suspect temporal autocorrelation is still an issue.

# Homoscedasticity of the errors (constant variance AKA homogeneity of variance)
# - Plot fitted values vs sqrt of sd of residuals
# - But this is more for guassian dist. This is relaxed for non-gaussian dist.
# - Look at quantile residuals
check_model(eelgrass_pc_glmm, check = c("homogeneity"))
# Clear curvature in fitted values vs sqrt of sd of residuals plot. But for non-gaussian 
# we should actually check residuals QQ plot
simulateResiduals(eelgrass_pc_glmm, plot = TRUE)
check_model(eelgrass_pc_glmm, check = c("qq"))
# looks good

# Overdispersion
# - Is observed variance higher than the variance of the theoretical model?
simulateResiduals(eelgrass_pc_glmm, plot = TRUE)
testDispersion(eelgrass_pc_glmm)

# Endogeneity - Predictor is not correlated with the error term
# - For glmm, this is predominantly argued logically, based on theory/DAG
#   - This may not pass, because SST anomaly could be correlated with storm events 
#     (e.g. mixing impacting SST and eelgrass), turbidity from chlorophyll, and other error. 
#     But generally, we assume the majority of factors influencing eelgrass cover would not
#     be correlated with SST anomaly, such as boating/fishing activity...? Including year as 
#     random effect also helps account for storm events. 

# Normality of the error distribution
check_model(eelgrass_pc_glmm)
# looks good

# Normality of random effects
check_model(eelgrass_pc_glmm, check = c("reqq"))
# looks good

# Minimal outlier influence
simulateResiduals(eelgrass_pc_glmm, plot = TRUE)
testOutliers(eelgrass_pc_glmm)
# no clear outliers


# ---- Evaluate Aggregated Percent Cover Model ----
tidy(eelgrass_pc_glmm)
# r2(eelgrass_pc_glmm)
# summary(eelgrass_pc_glmm)
# Anova(eelgrass_pc_glmm)


# ---- Visualize Aggregated Percent Cover Model - 1 ----
estimate_relation(eelgrass_pc_glmm) |> plot(show_data = TRUE)
visreg(eelgrass_pc_glmm, "sst_centered_prev", scale = "response")

# ---- Make Data Frame of Predicted Values for Model Plotting - 1 ----
 
# Create a prediction df following guidance from bbolker (see glmmTMB): 
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#glmmtmb 

# Create equence of SST anomalies 
sst_seq <- seq(min(eelgrass_data$sst_centered_prev),
               max(eelgrass_data$sst_centered_prev),
               length.out = 100)

# List of all sites
sites <- unique(eelgrass_data_aggregated$site)
# years <- unique(eelgrass_data_aggregated$year)

# Create a prediction grid
pred_grid_ag <- expand.grid(
  sst_centered_prev = sst_seq,
  site = NA,
  year = NA
)

# Get predictions with standard errors on response scale
preds_ag <- predict(
  eelgrass_pc_glmm,
  newdata = pred_grid_ag,
  type = "response",
  se.fit = TRUE,
)

# Attach predicted values and intervals to grid
pred_grid_ag$predicted <- preds_ag$fit
pred_grid_ag$se <- preds_ag$se.fit
pred_grid_ag$lower <- pred_grid_ag$predicted - 2 * pred_grid_ag$se
pred_grid_ag$upper <- pred_grid_ag$predicted + 2 * pred_grid_ag$se


# ---- Make Data Frame of Predicted Values for Model Plotting - 2 ----

# Create a prediction grid by site
pred_grid_ag_site <- expand.grid(
  sst_centered_prev = sst_seq,
  site = sites,
  year = NA
)

# # Create a prediction grid by year
# pred_grid_ag_year <- expand.grid(
#   sst_centered_prev = sst_seq,
#   site = NA,
#   year = years
# )

# Get predictions with standard errors on response scale - site
preds_ag <- predict(
  eelgrass_pc_glmm,
  newdata = pred_grid_ag_site,
  type = "response",
  se.fit = TRUE,
)

# # Get predictions with standard errors on response scale - year
# preds_ag <- predict(
#   eelgrass_pc_glmm,
#   newdata = pred_grid_ag_year,
#   type = "response",
#   se.fit = TRUE,
# )

# Attach predicted values and intervals to grid - site
pred_grid_ag_site$predicted <- preds_ag$fit
pred_grid_ag_site$se <- preds_ag$se.fit
pred_grid_ag_site$lower <- pred_grid_ag_site$predicted - 2 * pred_grid_ag_site$se
pred_grid_ag_site$upper <- pred_grid_ag_site$predicted + 2 * pred_grid_ag_site$se

# # Attach predicted values and intervals to grid - year
# pred_grid_ag_year$predicted <- preds_ag$fit
# pred_grid_ag_year$se <- preds_ag$se.fit
# pred_grid_ag_year$lower <- pred_grid_ag_year$predicted - 2 * pred_grid_ag_year$se
# pred_grid_ag_year$upper <- pred_grid_ag_year$predicted + 2 * pred_grid_ag_year$se


# ---- Visualize Aggregated Percent Cover Model - 2 ----
# Plot model across all data, not separating by site
ggplot(eelgrass_data_aggregated, aes(x = sst_centered_prev, y = percent_cover)) +
  geom_point(alpha = 0.3) +
  geom_line(
    data = pred_grid_ag,
    mapping = aes(y = predicted),
    linewidth = 1.5) +
  geom_ribbon(
    data = pred_grid_ag,
    mapping = aes(x = sst_centered_prev, y = predicted,
                  ymin = lower, ymax = upper),
    alpha = 0.2) +
  geom_line(
    data = pred_grid_ag_site, 
    mapping = aes(y = predicted, color = site),
    alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Eelgrass Percent Cover Predicted by SST Anomaly",
    x = "SST Anomaly (°C)",
    y = "Percent Cover",
    color = "Site",
    caption = "Ribbon = 95% CI"
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme(plot.caption.position = "plot", 
        legend.position = "none")

# ---- Visualize Aggregated Percent Cover Model - 3 ----
# Plot model across all data, faceting by site
ggplot(eelgrass_data_aggregated, 
       aes(x = sst_centered_prev, y = percent_cover, color = site)) +
  geom_point(alpha = 0.6) +
  geom_line(
    data = pred_grid_ag_site, 
    mapping = aes(y = predicted, group = site)) +
  geom_ribbon(
    data = pred_grid_ag_site,
    aes(x = sst_centered_prev, y = predicted,
        ymin = lower, ymax = upper, 
        fill = site, group = site),
    linetype = 0, alpha = 0.2,
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Eelgrass Percent Cover Predicted by SST Anomaly",
    x = "SST Anomaly (°C)",
    y = "Percent Cover",
    color = "Site",
    fill = "Site",
    caption = "Ribbon = 95% CI"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(plot.caption.position = "plot", 
        legend.position = "none") +
  scale_y_continuous(labels = scales::label_percent()) +
  facet_wrap(~site)

# # ---- Visualize Aggregated Percent Cover Model - 4----
# # Plot model across all data, faceting by year
# ggplot(eelgrass_data_aggregated, 
#        aes(x = sst_centered_prev, y = percent_cover, color = year)) +
#   geom_point(alpha = 0.6) +
#   geom_line(
#     data = pred_grid_ag_year, 
#     mapping = aes(y = predicted, group = year)) +
#   geom_ribbon(
#     data = pred_grid_ag_year,
#     aes(x = sst_centered_prev, y = predicted,
#         ymin = lower, ymax = upper, 
#         fill = year, group = year),
#     linetype = 0, alpha = 0.2,
#   ) +
#   labs(
#     title = "Predicted vs Observed Eelgrass Percent Cover with CI",
#     x = "Centered SST Anomaly (Previous Summer)",
#     y = "Percent Cover",
#     color = "Year",
#     fill = "Year",
#     caption = "Ribbons depict 95% CI"
#   ) +
#   theme(plot.caption.position = "plot") +
#   facet_wrap(~year)


# ---- Attempt at GAM - Aggregated Percent Cover ----

# eelgrass_data_aggregated <- eelgrass_data_aggregated |>
#   mutate(site = as.factor(site))
# 
# eelgrass_pc_gam <- gam(
#   # percent_cover ~ s(sst_centered_prev, by = site, bs = "tp") + s(year, bs = "re") + s(site, bs = "re"),
#   percent_cover ~ s(sst_centered_prev) + s(year, bs = "re") + s(site, bs = "re"),
#   data = eelgrass_data_aggregated,
#   family = betar(link = "logit"),
#   method = "REML"
# )
# # Warning that saturated likelihood may be inaccurate does not affect model fitting or inference
# # Rather it is only a concern if comparing models using AIC or related metrics
# # where deviance comparisons may be inaccurate
# 
# draw(eelgrass_pc_gam)
# gam.check(eelgrass_pc_gam)
# summary(eelgrass_pc_gam)
# r2(eelgrass_pc_gam)


# ---- Quadrat Presence Absence Model ----

# Build Model to predict presence absence:
# - Outcome is presence or absence (1 or 0)
# - Predictor is the centered mean sst from the previous year
# - Included a spline
# - Binomial error distribution because outcome is 1 or 0 (presence absence)
# - Logit link because working with percentages bounded between 0 & 1 (probability of presence)
# - Random effects include year, site, and site_quadrat, as these are expected to cause 
#   group-level variability. Year and site_quadrat are only expected to impact intercept, not
#   the influence of sst on presence
eelgrass_pa_glmm <- glmmTMB(
  presence ~ sst_centered_prev + (1 | year) + (1 + sst_centered_prev | site) + (1 | site_quadrat),
  family = binomial(link = "logit"),
  data = eelgrass_data)

# # Include year as predictor version - doesn't converge?
# eelgrass_pa_glmm <- glmmTMB(
#   presence ~ sst_centered_prev + year + (1 | site) + (1 | site_quadrat),
#   family = binomial(link = "logit"),
#   data = eelgrass_data)

# # Mundlak version
# eelgrass_pa_glmm <- glmmTMB(
#   presence ~ sst_centered_prev + mean_sst_prev + year + (1 | year) + (1 | site_quadrat),
#   family = binomial(link = "logit"),
#   data = eelgrass_data)

# # Spline version
# eelgrass_pa_glmm_spline <- glmmTMB(
#   presence ~ splines::ns(sst_centered_prev, df = 2) + (1 | year) + (1 + sst_centered_prev | site) + (1 | site_quadrat),
#   family = binomial(link = "logit"),
#   data = eelgrass_data)
# 
# AIC(eelgrass_pa_glmm, eelgrass_pa_glmm_spline)


# ---- Check Assumptions - Presence Absence Model ----
check_model(eelgrass_pa_glmm)

# Validity - Do X and Y Reflect Concepts I'm interested In
# - Conceptual check
#   - Yes, we are interested in water temperature, and sst anomalies can serve as a proxy for
#     bottom water anomalies in shallow areas where eelgrass grows.

# Representativeness: Does Your Data Represent the Population?
# - Conceptual check 
#   - Our sample we are working with only represents the population of the most shallow growing 
#     eelgrass on the Eastern US. Variety of different mean temps.

# Model captures features in the data
# - Posterior predictive check
check_predictions(eelgrass_pa_glmm)
# - Looks good!

# Linearity
# - Assuming linearity between the predictor and the log-odds (for binomial) 
#   or logit-transformed mean (for beta regression)
# - Simulate residual
simulateResiduals(eelgrass_pa_glmm, plot = TRUE)
check_model(eelgrass_pa_glmm, check = c("binned_residuals"))
# - Residuals are inconsistent, so linearity likely fails.
# - I tried adding splines. Spline with df = 2 seems better, but any higher df was bad
#   and AIC was worse

# Additivity
# - Assumption that fixed effects impact is additive - they don't affect each other
#   - All set here, because I only have one fixed effect

# Independence of errors
# - Are all replicates truly independent
#   - Samples will be correlated year to year within sites, and between neighboring quadrats. 
#   - Butinclusion of random effects for year, site, and quadrat takes care of this.
#   - Look at fitted values vs residuals

# Independence of errors
# - Are all replicates truly independent
#   - Look at fitted values vs residuals
simulateResiduals(eelgrass_pa_glmm, plot = TRUE)
testResiduals(eelgrass_pa_glmm)
#   - This suggest OK
#   - However, logically samples will be correlated year to year, within site, 
#     and between neighboring quadrats.
#   - Inclusion of random effects for year, site, and site_quadrat attempts to handle how
#     different years, sites, and quadrats behave, but I suspect temporal and within-site 
#     autocorrelation is still an issue.

# Homoscedasticity of the errors (constant variance AKA homogeneity of variance)
# - Plot fitted values vs sqrt of sd of residuals
check_model(eelgrass_pa_glmm, check = c("homogeneity"))
# - This is the least "flat and horizontal" line ever. But for non-gaussian this is
#   relaxed. 
# - We should actually check residuals QQ plot
# - Look at quantile residuals
simulateResiduals(eelgrass_pa_glmm, plot = TRUE)
check_model(eelgrass_pa_glmm, check = c("qq"))
# looks good

# Overdispersion
# - Is observed variance higher than the variance of the theoretical model?
simulateResiduals(eelgrass_pa_glmm, plot = TRUE)
testDispersion(eelgrass_pa_glmm)
# - I think this looks okay

# Endogeneity - Predictor is not correlated with the error term
# - For glmm, this is predominantly argued logically, based on theory/DAG
#   - This may not pass, because SST anomaly could be correlated with storm events 
#     (e.g. mixing impacting SST and eelgrass), turbidity from chlorophyll, and other error. 
#     But generally, we assume the majority of factors influencing eelgrass cover would not
#     be correlated with SST anomaly, such as boating/fishing activity...? Including year as 
#     random effect also helps account for storm events. 

# Normality of the error distribution
check_model(eelgrass_pa_glmm)
# - looks good

# Normality of random effects
check_model(eelgrass_pa_glmm, check = c("reqq"))
# looks good

# Minimal outlier influence
simulateResiduals(eelgrass_pa_glmm, plot = TRUE)
testOutliers(eelgrass_pa_glmm)
# no clear outliers


# ---- Evaluate Aggregated Percent Cover Model ----
tidy(eelgrass_pa_glmm)
# summary(eelgrass_pa_glmm)

# ---- Visualize Presence Absence Model - 1 ----
estimate_relation(eelgrass_pa_glmm) |> plot(show_data = TRUE)
visreg(eelgrass_pa_glmm, "sst_centered_prev", scale = "response")

# ---- Make Data Frame of Predicted Values for Model Plotting - 3 ----

# Create a prediction df for presence absence
pred_grid_pa <- expand.grid(
  sst_centered_prev = sst_seq,
  site = NA,
  year = NA,
  site_quadrat = NA
)

# Get predictions with standard errors on response scale
preds_pa <- predict(
  eelgrass_pa_glmm,
  newdata = pred_grid_pa,
  type = "response",
  se.fit = TRUE,
)

# Attach predicted values and intervals to grid
pred_grid_pa$predicted <- preds_pa$fit
pred_grid_pa$se <- preds_pa$se.fit
pred_grid_pa$lower <- pred_grid_pa$predicted - (2 * pred_grid_pa$se)
pred_grid_pa$upper <- pred_grid_pa$predicted + (2 * pred_grid_pa$se)


# ---- Make Binned Raw Data Frame for Plotting Proportions Present - 1 ----

# Create binned data frame to plot proportion present within temperature anomaly bins
binned_eg_data <- eelgrass_data |> 
  # using cut_number to put the about same number of sites in each bin
  mutate(sst_bin = cut_number(sst_centered_prev, n = 20)) |> 
  group_by(sst_bin) |> 
  summarize(
    sst_anomaly = mean(sst_centered_prev, na.rm = TRUE),
    prop_present = mean(presence),
    n = n(),
    .groups = "drop"
  )

# ---- Make Data Frame of Predicted Values for Model Plotting - 4 ----

# Create a site level prediction df
pred_grid_pa_site <- expand.grid(
  sst_centered_prev = sst_seq,
  site = sites,
  year = NA,
  site_quadrat = NA
)

# Get predictions with standard errors on response scale
preds_pa <- predict(
  eelgrass_pa_glmm,
  newdata = pred_grid_pa_site,
  type = "response",
  se.fit = TRUE,
)

# Attach predicted values and intervals to grid
pred_grid_pa_site$predicted <- preds_pa$fit
pred_grid_pa_site$se <- preds_pa$se.fit
pred_grid_pa_site$lower <- pred_grid_pa_site$predicted - (2 * pred_grid_pa_site$se)
pred_grid_pa_site$upper <- pred_grid_pa_site$predicted + (2 * pred_grid_pa_site$se)

# ---- Make Binned Raw Data Frame for Plotting Proportions Present - 2 ---- 
# Bins within site level
binned_site_eg_data <- eelgrass_data |> 
  # using cut to maintain consistent temperature ranges per bin
  mutate(sst_bin = cut(sst_centered_prev, breaks = 10)) |> 
  group_by(site, sst_bin) |> 
  summarize(
    sst_anomaly = mean(sst_centered_prev, na.rm = TRUE),
    prop_present = mean(presence),
    n = n(),
    .groups = "drop"
  )

# ---- Visualize Presence Absence Model - 2 ----

# Plot modeled probability of presence with proportion of presence binned by sst anomaly
ggplot(binned_eg_data, aes(x = sst_anomaly, y = prop_present)) +
  geom_point(shape = 12) +
  geom_line(pred_grid_pa,
            mapping = aes(x = sst_centered_prev,
                          y = predicted),
            linewidth = 1.5) +
  geom_ribbon(pred_grid_pa,
              mapping = aes(x = sst_centered_prev, y = predicted,
                            ymin = lower, ymax = upper),
              linetype = 0, alpha = 0.2) +
  geom_line(
    data = pred_grid_pa_site, 
    mapping = aes(x = sst_centered_prev, y = predicted, color = site),
    alpha = 0.4) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(title = "Modeled Probability of Presence Absence",
       subtitle = "[grids represent bins of site-years with similar SST anomalies]",
       y = "Modeled Probability (line) & Actual Proportion (grids) of Presence",
       x = "SST Anomaly (°C)",
       caption = "Ribbon = 95% CI"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(plot.caption.position = "plot",
        legend.position = "none")


# ---- Visualize Presence Absence Model - 3 ----

# Plot modeled probability of presence with proportion of sites with presence binned by sst anomaly
ggplot(binned_site_eg_data, aes(x = sst_anomaly, y = prop_present, color = site)) +
  geom_point(shape = 12) +
  geom_line(pred_grid_pa_site,
            mapping = aes(x = sst_centered_prev, y = predicted, color = site)) +
  geom_ribbon(pred_grid_pa_site,
              mapping = aes(x = sst_centered_prev, y = predicted,
                            ymin = lower, ymax = upper,
                            fill = site),
              linetype = 0, alpha = 0.2) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~site) +
  labs(y = "Modeled Probability (line) & Actual Proportion (grids) of Presence",
       x = "SST Anomaly (°C)",
       caption = "Ribbon = 95% CI",
       color = "Site",
       fill = "Site"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme(plot.caption.position = "plot")



# ---- Visualize Presence Absence Model - 4 ----

# # Plot showing actual presence absence (1 or 0) data behind model probability predictions
# # NEEDS REVISION
# ggplot(eelgrass_data, aes(x = sst_centered_prev, y = presence, color = site)) +
#   geom_point(alpha = 0.6) +
#   geom_line(data = pred_grid_pa, aes(y = predicted, group = site)) +
#   geom_ribbon(
#     data = pred_grid_pa,
#     aes(x = sst_centered_prev, y = predicted, ymin = lower, ymax = upper, fill = site, group = site),
#     alpha = 0.2,
#   ) +
#   labs(
#     title = "Predicted vs Observed Eelgrass Presence Absence with CI",
#     x = "Centered SST Anomaly (Previous Summer)",
#     y = "Presence (1 = Presence, 0 = Absence)",
#     color = "Site",
#     fill = "Site"
#   ) +
#   coord_cartesian(ylim = c(0, 1)) +
#   facet_wrap(~site) 


# ---- Attempt at GAM - Presence Absence ----
### Attempt at GAM
# 
# eelgrass_data <- eelgrass_data |>
#   mutate(site = as.factor(site), 
#          site_quadrat = as.factor(site_quadrat))
# 
# eelgrass_pa_gam <- gam(
#   # presence ~ s(sst_centered_prev, by = site, bs = "tp") + 
#   presence ~ s(sst_centered_prev) + 
#     s(year, bs = "re") + s(site, bs = "re") + s(site_quadrat, bs = "re"),
#   data = eelgrass_data,
#   family = binomial(link = "logit"),
#   method = "REML"
# )
# 
# draw(eelgrass_pa_gam)
# gam.check(eelgrass_pa_gam)
# summary(eelgrass_pa_gam)
# r2(eelgrass_pa_gam)


# ---- Quadrat Percent Cover *Given Coverage* Model ----
# 
# # Build another model to predict abundance, at the quadrat level, given any coverage:
# # - Outcome is percent cover 
# # - Predictor is the centered mean sst from the previous year
# # - Beta error distribution because % cover data
# # - Logit link because working with percentages bounded between 0 & 1
# # - REMOVE zero values, as this is given coverage
# # - Random effects include year, site, and site_quadrat, as these are expected to cause 
# #   group-level variability. Year and site_quadrat are only expected to impact intercept, not
# #   the influence of sst on presence. 
# #   Removing random effect impact on slope because model won't converge otherwise...
# 
# 
# eelgrass_pc_glmm_quadrat <- glmmTMB(
#   # percent_cover ~ sst_centered_prev + (1 | year) + (1 + sst_centered_prev | site) + (1 | site_quadrat),
#   percent_cover ~ splines::ns(sst_centered_prev, df = 2) + (1 | year) + (1 + sst_centered_prev | site) + (1 | site_quadrat),
#   family = beta_family(link = "logit"),
#   # ziformula = ~1,
#   data = eelgrass_data_non_ag
# )
# 
# # Check assumptions
# check_model(eelgrass_pc_glmm_quadrat)
# 
# # Simulate residuals
# simulateResiduals(fittedModel = eelgrass_pc_glmm_quadrat, plot = TRUE)
# 
# 
# # Evaluate model
# tidy(eelgrass_pc_glmm_quadrat)
# # summary(eelgrass_pc_glmm_quadrat)
# 
# # Visualize model
# estimate_relation(eelgrass_pc_glmm_quadrat) |> plot(show_data = TRUE)
# visreg(eelgrass_pc_glmm_quadrat, "sst_centered_prev", scale = "response")
# 
# 
# ### Attempt at GAM
# 
# # eelgrass_data_non_ag <- eelgrass_data_non_ag |>
# #   mutate(site = as.factor(site), 
# #          site_quadrat = as.factor(site_quadrat))
# # 
# # 
# # eelgrass_pc_gam_quadrat <- gam(
# #   percent_cover ~ s(sst_centered_prev, by = site, bs = "tp") +
# #   # percent_cover ~ s(sst_centered_prev) + 
# #     s(year, bs = "re") + s(site, bs = "re") + s(site_quadrat, bs = "re"),
# #   data = eelgrass_data_non_ag,
# #   family = betar(link = "logit"),
# #   method = "REML"
# # )
# # 
# # draw(eelgrass_pc_gam_quadrat)
# # gam.check(eelgrass_pc_gam_quadrat)
# # summary(eelgrass_pc_gam_quadrat)
# # r2(eelgrass_pc_gam_quadrat)

