#'--------------------------
#' Directed Acyclic Graphs (DAG) of Eelgrass System
#'--------------------------

# Load Libraries for DAGs
library(dagitty)
library(ggdag)
library(ggplot2)

# Load other Libraries
library(dplyr)

# Build simplified eeglrass DAG, where all effects other than SST 
# are included in unmeasured "site effects"
eelgrass_dag_initial <- dagify(
  eelgrass_cover ~ sst + latitude_effects + site_effects,
  sst ~ latitude_effects + site_effects,
  site_effects ~ latitude_effects
  )

# Initial look at DAG
eelgrass_dag_initial |> plot()

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

# Plot fancy DAG
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