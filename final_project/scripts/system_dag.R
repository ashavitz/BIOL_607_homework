#'--------------------------
#' Directed Acyclic Graphs (DAG) of Eelgrass System
#'--------------------------

# Load Libraries for DAGs
library(dagitty)
library(ggdag)
library(ggplot2)

# Load other Libraries
library(dplyr)

# ---- DAG Set Up ----

# Build eeglrass DAG, where all effects other than SST are included in 
# unmeasured random site effect, year effects, or quadrat effects
# Include site_effects impact SST, which would cause endogeneity problems for my models
eelgrass_dag_initial <- dagify(
  eelgrass_cover ~ sst + site_effects + year + quadrat_effects,
  quadrat_effects ~ site_effects, 
  sst ~ site_effects
)

# Initial look at DAG
eelgrass_dag_initial |> plot()


# Check conditional independence relationships
impliedConditionalIndependencies(eelgrass_dag_initial)

# SST -> eelgrass cover (direct effects)
adjustmentSets(eelgrass_dag_initial,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "direct",
               type = "minimal")

# SST -> eelgrass cover (total effects)
adjustmentSets(eelgrass_dag_initial,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "total",
               type = "all")


# If we calculate the SST anomaly, that is the annual SST relative to the average SST
# at each site. This breaks the causal connection from site_effects to SST (anomaly), because
# the site_effect would theoretically be consistent (e.g. cold current)

# DAG with SST Anomaly
eelgrass_dag_anomaly <- dagify(
  eelgrass_cover ~ sst_anomaly + site_effects + year + quadrat_effects,
  quadrat_effects ~ site_effects
)

# Look at DAG
eelgrass_dag_anomaly |> plot()


## Beautify DAG
coordinates(eelgrass_dag_anomaly) <- list(
  x = c(sst_anomaly = 0,
        site_effects = 3,
        quadrat_effects = 6,
        year = 6,
        eelgrass_cover = 6),
  
  y = c(sst_anomaly = 2,
        site_effects = 4,
        quadrat_effects = 4,
        year = 0,
        eelgrass_cover = 2)
)

# Tidy and relabel
eelgrass_dag_anomaly_tidy <- tidy_dagitty(eelgrass_dag_anomaly)
eelgrass_dag_anomaly_tidy$data$name <- dplyr::recode(
  eelgrass_dag_anomaly_tidy$data$name,
  eelgrass_cover = "Eelgrass Cover",
  sst_anomaly = "SST Anomaly",
  site_effects = "Site Effects",
  quadrat_effects = "Quadrat Effects",
  year = "Year"
)

# Adjust arrows to reduce overlap
eelgrass_dag_anomaly_tidy$data <- eelgrass_dag_anomaly_tidy$data |> 
  mutate(
    xstart = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 3.11,
      name == "Site Effects" & to == "quadrat_effects" ~ 3.25,
      name == "SST Anomaly" ~ 0.3,
      name == "Year" ~ 6
    ),
    ystart = case_when(
      name == "Quadrat Effects" ~ 3.85,
      name == "Site Effects" & to == "eelgrass_cover" ~ 4,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "SST Anomaly" ~ 2,
      name == "Year" ~ 0.15
    ),
    xend = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 5.9,
      name == "Site Effects" & to == "quadrat_effects" ~ 5.7,
      name == "SST Anomaly" ~ 5.7,
      name == "Year" ~ 6
    ),
    yend = case_when(
      name == "Quadrat Effects" ~ 2.2,
      name == "Site Effects" & to == "eelgrass_cover" ~ 2.1,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "SST Anomaly" ~ 2,
      name == "Year" ~ 1.85
    )
  )

# ---- Final Anomaly DAG ----

# Plot the DAG
ggplot(eelgrass_dag_anomaly_tidy,
       aes(x = x, y = y,
           xend = xend, yend = yend)) +
  geom_dag_node(shape = 5, size = 25, color = "blue") +
  geom_dag_text(size = 3, color = "black") +
  geom_dag_edges(
    aes(x = xstart, y = ystart),
    edge_width = 1.2
  ) +
  xlim(c(-1, 10)) +
  ylim(c(-0.5, 5)) +
  theme_dag_blank()


# ---- Confounded DAG ----

## Beautify DAG
coordinates(eelgrass_dag_initial) <- list(
  x = c(sst = 0,
        site_effects = 3,
        quadrat_effects = 6,
        year = 6,
        eelgrass_cover = 6),
  
  y = c(sst = 2,
        site_effects = 4,
        quadrat_effects = 4,
        year = 0,
        eelgrass_cover = 2)
)

# Tidy and relabel
eelgrass_dag_initial_tidy <- tidy_dagitty(eelgrass_dag_initial)
eelgrass_dag_initial_tidy$data$name <- dplyr::recode(
  eelgrass_dag_initial_tidy$data$name,
  eelgrass_cover = "Eelgrass Cover",
  sst = "SST",
  site_effects = "Site Effects",
  quadrat_effects = "Quadrat Effects",
  year = "Year"
)

# Adjust arrows to reduce overlap
eelgrass_dag_initial_tidy$data <- eelgrass_dag_initial_tidy$data |> 
  mutate(
    xstart = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 3.11,
      name == "Site Effects" & to == "quadrat_effects" ~ 3.25,
      name == "Site Effects" & to == "sst" ~ 2.95,
      name == "SST" ~ 0.3,
      name == "Year" ~ 6
    ),
    ystart = case_when(
      name == "Quadrat Effects" ~ 3.85,
      name == "Site Effects" & to == "eelgrass_cover" ~ 4,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "Site Effects" & to == "sst" ~ 3.95,
      name == "SST" ~ 2,
      name == "Year" ~ 0.15
    ),
    xend = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 5.9,
      name == "Site Effects" & to == "quadrat_effects" ~ 5.7,
      name == "Site Effects" & to == "sst" ~ 0.1,
      name == "SST" ~ 5.7,
      name == "Year" ~ 6
    ),
    yend = case_when(
      name == "Quadrat Effects" ~ 2.2,
      name == "Site Effects" & to == "eelgrass_cover" ~ 2.1,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "Site Effects" & to == "sst" ~ 2,
      name == "SST" ~ 2,
      name == "Year" ~ 1.85
    )
  )

# ---- Confounded DAG Final ----

# Plot the DAG
ggplot(eelgrass_dag_initial_tidy,
       aes(x = x, y = y,
           xend = xend, yend = yend)) +
  geom_dag_node(shape = 5, size = 25, color = "blue") +
  geom_dag_text(size = 3, color = "black") +
  geom_dag_edges(
    aes(x = xstart, y = ystart),
    edge_width = 1.2
  ) +
  xlim(c(-1, 10)) +
  ylim(c(-0.5, 5)) +
  theme_dag_blank()




# ---- Full Final DAG - Setup ----

# Build eeglrass DAG, where all effects other than SST are included in 
# unmeasured random site effect, year effects, or quadrat effects
# Include site_effects impact SST, but SST anomaly is endogenous
eelgrass_dag_final <- dagify(
  eelgrass_cover ~ sst + sst_anomaly + site_effects + year + quadrat_effects,
  quadrat_effects ~ site_effects, 
  sst ~ site_effects
)

# Initial look at DAG
eelgrass_dag_final |> plot()


# SST -> eelgrass cover (direct effects)
adjustmentSets(eelgrass_dag_final,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "direct",
               type = "minimal")

# SST -> eelgrass cover (total effects)
adjustmentSets(eelgrass_dag_final,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "total",
               type = "all")


# If we calculate the SST anomaly, that is the annual SST relative to the average SST
# at each site. Including average mean SST breaks the causal connection from site_effects to 
# SST (anomaly), allowing for exogenous link from SST anomaly to eelgrass cover


## Beautify DAG
coordinates(eelgrass_dag_final) <- list(
  x = c(sst = 0,
        sst_anomaly = 0,
        site_effects = 3,
        quadrat_effects = 6,
        year = 6,
        eelgrass_cover = 6),
  
  y = c(sst = 2,
        sst_anomaly = 0,
        site_effects = 4,
        quadrat_effects = 4,
        year = 0,
        eelgrass_cover = 2)
)

# Tidy and relabel
eelgrass_dag_final_tidy <- tidy_dagitty(eelgrass_dag_final)
eelgrass_dag_final_tidy$data$name <- dplyr::recode(
  eelgrass_dag_final_tidy$data$name,
  eelgrass_cover = "Eelgrass Cover",
  sst = "Avg-Mean SST",
  sst_anomaly = "SST Anomaly",
  site_effects = "Site Effects",
  quadrat_effects = "Quadrat Effects",
  year = "Year"
)

# Adjust arrows to reduce overlap
eelgrass_dag_final_tidy$data <- eelgrass_dag_final_tidy$data |> 
  mutate(
    xstart = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 3.11,
      name == "Site Effects" & to == "quadrat_effects" ~ 3.25,
      name == "Site Effects" & to == "sst" ~ 2.95,
      name == "Avg-Mean SST" ~ 0.3,
      name == "SST Anomaly" ~ 0.1,
      name == "Year" ~ 6
    ),
    ystart = case_when(
      name == "Quadrat Effects" ~ 3.85,
      name == "Site Effects" & to == "eelgrass_cover" ~ 4,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "Site Effects" & to == "sst" ~ 3.95,
      name == "Avg-Mean SST" ~ 2,
      name == "SST Anomaly" ~ 0.1,
      name == "Year" ~ 0.15
    ),
    xend = case_when(
      name == "Quadrat Effects" ~ 6,
      name == "Site Effects" & to == "eelgrass_cover" ~ 5.9,
      name == "Site Effects" & to == "quadrat_effects" ~ 5.7,
      name == "Site Effects" & to == "sst" ~ 0.1,
      name == "Avg-Mean SST" ~ 5.7,
      name == "SST Anomaly" ~ 6,
      name == "Year" ~ 6
    ),
    yend = case_when(
      name == "Quadrat Effects" ~ 2.2,
      name == "Site Effects" & to == "eelgrass_cover" ~ 2.1,
      name == "Site Effects" & to == "quadrat_effects" ~ 4,
      name == "Site Effects" & to == "sst" ~ 2,
      name == "Avg-Mean SST" ~ 2,
      name == "SST Anomaly" ~ 1.85,
      name == "Year" ~ 1.85
    )
  )

# ---- Full Final DAG ----

# Plot the DAG
ggplot(eelgrass_dag_final_tidy,
       aes(x = x, y = y,
           xend = xend, yend = yend)) +
  geom_dag_node(shape = 5, size = 25, color = "blue") +
  geom_dag_text(size = 3, color = "black") +
  geom_dag_edges(
    aes(x = xstart, y = ystart),
    edge_width = 1.2
  ) +
  xlim(c(-1, 10)) +
  ylim(c(-0.5, 5)) +
  theme_dag_blank()






# ---- Previous Versions ----
#### PREVIOUS VERSIONS ##############


# Build simplified eeglrass DAG, where all effects other than SST 
# are included in unmeasured "site effects"
eelgrass_dag_initial <- dagify(
  eelgrass_cover ~ sst + latitude_effects + site_effects,
  sst ~ latitude_effects + site_effects,
  site_effects ~ latitude_effects
  )

# Initial look at DAG
eelgrass_dag_initial |> plot()


# Check conditional independence relationships
impliedConditionalIndependencies(eelgrass_dag_initial)

# SST -> eelgrass cover (direct effects)
adjustmentSets(eelgrass_dag_initial,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "direct",
               type = "minimal")

# SST -> eelgrass cover (total effects)
adjustmentSets(eelgrass_dag_initial,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "total",
               type = "all")


# Simplified DAG with just eelgrass cover, SST, and site_effects controlled by SST
eelgrass_dag_simple <- dagify(
  eelgrass_cover ~ sst + site_effects,
  site_effects ~ sst
)
# Look at DAG
eelgrass_dag_simple |> plot()


# TODO - Beautify DAGs

coordinates(eelgrass_dag_simple) <- list(
  x = c(eelgrass_cover = 10, sst = 0,  site_effects = 5),

  y = c(eelgrass_cover = 0, sst = 0, site_effects = 1)
)

# Make a df of DAG
eelgrass_dag_simple_tidy <- tidy_dagitty(eelgrass_dag_simple)
eelgrass_dag_simple_tidy$data$name <- dplyr::recode(
  eelgrass_dag_simple_tidy$data$name,
  eelgrass_cover = "Eelgrass Cover",
  sst = "SST",
  site_effects = "Site Effects"
)

# Shorted arrows to avoid overlap
eelgrass_dag_simple_tidy$data <- eelgrass_dag_simple_tidy$data |> 
  mutate(
    xstart = x + (0.025*xend),
    ystart = 0.95*y,
    xend = 0.95*xend,
    yend = 0.95*yend
  )
# Good enough...

# ---- Plot fancy DAG ----
ggplot(eelgrass_dag_simple_tidy,
       aes(x = x, y = y,
           xend = xend, yend = yend)) +
  geom_dag_node(shape = 5, size = 25, color = "blue") +
  geom_dag_text(size = 3, color = "black") +
  geom_dag_edges(
    aes(x = xstart, y = ystart),
    edge_width = 1.2,
  ) +
  xlim(c(-0.5, 10.5)) +
  ylim(c(-0.5, 2)) +
  theme_dag_blank()



# THE LOGIC BELOW IS WRONG. ARCHIVING HERE AT BOTTOM.

# Since in reality, SST would impact site effects, we might have SST <-> site_effects.
# But since this is a DAG, we can't have that feedback, so we can separate out the site_effects
# which are most likely to impact SST (), such as oceanography and weather, and other site
# effects which may be impacted by SST 
eelgrass_dag_expanded <- dagify(
  eelgrass_cover ~ sst + site_effects + oceanography + weather + latitude_effects,
  sst ~ oceanography + weather + latitude_effects,
  site_effects ~ sst + oceanography + weather,
  oceanography ~ latitude_effects,
  weather ~ latitude_effects
)
# Look at DAG
eelgrass_dag_expanded |> plot()

# Check conditional independence relationships
impliedConditionalIndependencies(eelgrass_dag_expanded)

# SST -> eelgrass cover (direct effects)
adjustmentSets(eelgrass_dag_expanded,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "direct",
               type = "minimal")
# If we wanted to know the direct of sst on eelgrass cover, we would have to control for 
# any site effect controlled by sst, which we can't do. 

# SST -> eelgrass cover (total effects)
adjustmentSets(eelgrass_dag_expanded,
               exposure = "sst",
               outcome = "eelgrass_cover",
               effect = "total",
               type = "all")
# But we can look at the total effect of sst on eelgrass cover, because we are able to control 
# for site effects not controlled by sst, such as weather, oceanography, and latitude by using
# longitudinal data (time series) at the same sites.

