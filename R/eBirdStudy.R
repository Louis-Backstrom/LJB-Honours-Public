# Libraries and Functions ----

library(tidyverse)
library(auk)
library(sf)
library(lubridate)
library(ggpubr)

source("R/CustomFunctions.R")

# Set-up ----

obs <- read_rds("Data/eBird Filtered Dataset/Observation.rds")                  # read in eBird observations
sam <- read_rds("Data/eBird Filtered Dataset/Sampling.rds")                     # read in eBird checklist information
hab <- read_sf("Data/NewSpatial/NewHabPA.gpkg")                                 # read in landscape spatial data

all_habitats <- hab %>%                                                         # get list of habitat types
  pull(habitat) %>% 
  unique()

hab_proportions <- hab %>%                                                      # get relative proportions of each habitat in landscape
  filter(brisbane == TRUE) %>% 
  st_make_valid() %>%
  mutate(area = st_area(.)) %>% 
  st_drop_geometry() %>% 
  group_by(habitat) %>% 
  summarise(area = as.numeric(sum(area))) %>% 
  ungroup() %>% 
  mutate(relative_proportion = area / sum(area)) %>% 
  select(-area)

loc <- sam %>%                                                                  # get habitat type for each location
  select(locality_id, latitude, longitude) %>%
  unique() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 3112) %>%
  st_join(., hab) %>% 
  st_drop_geometry() %>% 
  select(locality_id, habitat)

sp <- obs %>%                                                                   # get list of species to model (>1000 total observations)
  group_by(scientific_name) %>% 
  summarise(n_obs = n()) %>% 
  filter(n_obs >= 1000) %>% 
  pull(scientific_name)

zf <- auk_zerofill(x = obs, sampling_events = sam, collapse = TRUE) %>%         # get zerofilled presence-absence dataframe of all model species
  left_join(., loc) %>% 
  filter(scientific_name %in% sp) %>% 
  select(checklist_id, locality_id, latitude, longitude, observation_date,
         time_observations_started, observer_id, protocol_type, 
         duration_minutes, effort_distance_km, number_observers,
         scientific_name, observation_count, species_observed, habitat) %>% 
  mutate(year = year(observation_date),
         yday = yday(observation_date),
         aday = as.integer(as_date(observation_date) - as_date("2011-12-31")))

# Models ----

species_model_changes <- tibble(
  species = character(),
  habitat = character(),
  nullinformed = character(),
  semiinformed = character(),
  fullinformed = character(),
)

for (i in c(65, 106)) {                                                       # run each species
  species <- sp[i]
  
  zf_data <- zf %>%                                                             # species-specific dataframe
    filter(scientific_name == species) %>% 
    mutate(species_observed = as.integer(species_observed))
  
  nullinformed_glm <- glm(data = zf_data,                                       # run a basic, uninformed GLM
                          formula = species_observed ~ aday,
                          family = "binomial") 
  semiinformed_glm <- glm(data = zf_data,                                       # run a simple, semi-informed GLM
                          formula = species_observed ~ aday + habitat,
                          family = "binomial")
  fullinformed_glm <- glm(data = zf_data,                                       # run a complex, fully-informed GLM
                          formula = species_observed ~ aday * habitat,
                          family = "binomial")
  
  zf_data_plots <- expand_grid(                                                 # get data into a plottable format
    date = seq(as_date("2012-01-01"), as_date("2021-12-31"), by = "1 day"),
    habitat = all_habitats) %>% 
    mutate(year = year(date),
           yday = yday(date),
           aday = as.integer(as_date(date) - as_date("2011-12-31")))
  
  zf_data_plots$nullinformed_pred <- nullinformed_glm %>%                       # run uninformed model predictions on plot data
    predict(., zf_data_plots) %>% 
    plogis()
  zf_data_plots$semiinformed_pred <- semiinformed_glm %>%                       # run semi-informed model predictions on plot data
    predict(., zf_data_plots) %>% 
    plogis()
  zf_data_plots$fullinformed_pred <- fullinformed_glm %>%                       # run fully-informed model predictions on plot data
    predict(., zf_data_plots) %>% 
    plogis()
  
  zf_data_tests <- expand_grid(                                                 # get data into a plottable format
    date = c(as_date("2012-01-01"), as_date("2021-12-31")),
    habitat = all_habitats) %>% 
    mutate(year = year(date),
           yday = yday(date),
           aday = as.integer(as_date(date) - as_date("2011-12-31")))
  
  zf_data_tests$nullinformed_pred <- predict.glm(                               # run uninformed model predictions on plot data      
    nullinformed_glm, zf_data_tests, se.fit = TRUE)$fit 
  zf_data_tests$nullinformed_se <- predict.glm(
    nullinformed_glm, zf_data_tests, se.fit = TRUE)$se.fit
  
  zf_data_tests$semiinformed_pred <- predict.glm(                               # run semi-informed model predictions on plot data
    semiinformed_glm, zf_data_tests, se.fit = TRUE)$fit 
  zf_data_tests$semiinformed_se <- predict.glm(
    semiinformed_glm, zf_data_tests, se.fit = TRUE)$se.fit
  
  zf_data_tests$fullinformed_pred <- predict.glm(                               # run fully-informed model predictions on plot data
    fullinformed_glm, zf_data_tests, se.fit = TRUE)$fit 
  zf_data_tests$fullinformed_se <- predict.glm(
    fullinformed_glm, zf_data_tests, se.fit = TRUE)$se.fit
  
  plot_base <- ggplot(mapping = aes(x = date)) +                                # set up the basic (shared) elements of the plots
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    lims(y = c(0, 1)) +
    labs(x = "", y = "")
  
  plot_nullinformed <- plot_base +                                              # plot un-informed model
    geom_line(
      data = zf_data_plots, size = 1, alpha = 0.5,
      aes(y = nullinformed_pred, colour = habitat)) +
    geom_point(
      data = zf_data_tests, size = 1, alpha = 0.5,
      aes(y = plogis(nullinformed_pred), colour = habitat)) +
    geom_errorbar(
      data = zf_data_tests, size = 1, alpha = 0.5, width = 100,
      aes(ymin = plogis(nullinformed_pred - qnorm(0.975) * nullinformed_se),
          ymax = plogis(nullinformed_pred + qnorm(0.975) * nullinformed_se),
          colour = habitat)) +
    labs(x = "Year", y = "Reporting Rate", subtitle = "Null Model",
         colour = "Habitat Type")
  
  plot_semiinformed <- plot_base +                                              # plot semi-informed model
    geom_line(
      data = zf_data_plots, size = 1, alpha = 0.5,
      aes(y = semiinformed_pred, colour = habitat)) +
    geom_point(
      data = zf_data_tests, size = 1, alpha = 0.5,
      aes(y = plogis(semiinformed_pred), colour = habitat)) +
    geom_errorbar(
      data = zf_data_tests, size = 1, alpha = 0.5, width = 100,
      aes(ymin = plogis(semiinformed_pred - qnorm(0.975) * semiinformed_se),
          ymax = plogis(semiinformed_pred + qnorm(0.975) * semiinformed_se),
          colour = habitat)) +
    labs(x = "Year", y = "Reporting Rate", subtitle = "Additive Model")
  
  plot_fullinformed <- plot_base +                                              # plot fully-informed model
    geom_line(
      data = zf_data_plots, size = 1, alpha = 0.5,
      aes(y = fullinformed_pred, colour = habitat)) +
    geom_point(
      data = zf_data_tests, size = 1, alpha = 0.5, 
      aes(y = plogis(fullinformed_pred), colour = habitat)) +
    geom_errorbar(
      data = zf_data_tests, size = 1, alpha = 0.5, width = 100,
      aes(ymin = plogis(fullinformed_pred - qnorm(0.975) * fullinformed_se),
          ymax = plogis(fullinformed_pred + qnorm(0.975) * fullinformed_se),
          colour = habitat)) +
    labs(x = "Year", y = "Reporting Rate", subtitle = "Interactive Model")
  
  plot_combined <- ggarrange(                                                   # arrange plots
    ncol = 3, plot_nullinformed, plot_semiinformed, plot_fullinformed,
    legend = "bottom", common.legend = TRUE)
  
  plot_combined <- annotate_figure(p = plot_combined, top = text_grob(
    species, face = "bold.italic", size = 16))
  
  ggsave(filename = paste0("Figures/eBirdStudy/", species, ".svg"),             # save full plot
         plot = plot_combined, width = 8, height = 4)
  ggsave(filename = paste0("Figures/eBirdStudy/", species, ".png"),             
         plot = plot_combined, width = 8, height = 4)

  zf_data_tests <- zf_data_tests %>%                                            # determine overall change from start to finish
    group_by(habitat) %>% 
    summarise(
      nullinformed_difference = abs(diff(nullinformed_pred)),
      nullinformed_differencese = sqrt(sum(nullinformed_se^2)),
      nullinformed_sign = sign_chr(diff(nullinformed_pred)),
      semiinformed_difference = abs(diff(semiinformed_pred)),
      semiinformed_differencese = sqrt(sum(semiinformed_se^2)),
      semiinformed_sign = sign_chr(diff(semiinformed_pred)),
      fullinformed_difference = abs(diff(fullinformed_pred)),
      fullinformed_differencese = sqrt(sum(fullinformed_se^2)),
      fullinformed_sign = sign_chr(diff(fullinformed_pred))) %>% 
    mutate(
      nullinformed_sign = if_else(
        nullinformed_difference > qnorm(0.975) * nullinformed_differencese,
        nullinformed_sign, "n.s."),
      semiinformed_sign = if_else(
        semiinformed_difference > qnorm(0.975) * semiinformed_differencese,
        semiinformed_sign, "n.s."),
      fullinformed_sign = if_else(
        fullinformed_difference > qnorm(0.975) * fullinformed_differencese,
        fullinformed_sign, "n.s."),
      species = species) %>% 
    select(species, habitat, nullinformed_sign, 
           semiinformed_sign, fullinformed_sign)
  
  species_model_changes <- rbind(species_model_changes, zf_data_tests)
  
  print(species)
}

# Overall Changes ----

for (habitats in hab_proportions$habitat) {
  print(habitats)
  
  table <- species_model_changes %>%                                            # Table Null vs Semi
    mutate(across(.cols = 3:5, 
                  .fns = ~factor(.x, levels = c("-", "n.s.", "+")))) %>%
    filter(habitat == habitats) %>%
    select(nullinformed_sign, semiinformed_sign) %>% 
    table()
  
  print(table)
  
  table <- species_model_changes %>%                                            # Table Semi vs Full
    mutate(across(.cols = 3:5, 
                  .fns = ~factor(.x, levels = c("-", "n.s.", "+")))) %>%
    filter(habitat == habitats) %>%
    select(semiinformed_sign, fullinformed_sign) %>% 
    table()
  
  print(table)
  
  table <- species_model_changes %>%                                            # Table Full vs Null
    mutate(across(.cols = 3:5, 
                  .fns = ~factor(.x, levels = c("-", "n.s.", "+")))) %>%
    filter(habitat == habitats) %>%
    select(fullinformed_sign, nullinformed_sign) %>% 
    table()
  
  print(table)
}
