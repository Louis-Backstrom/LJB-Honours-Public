# Libraries and Functions ----

library(tidyverse) 
library(lubridate)
library(scales)
library(hms)
library(sf)
library(exactextractr)
library(terra)
library(raster)
library(ggpubr)

source("R/CustomFunctions.R")

extract <- raster::extract
select <- dplyr::select

# Set-up ----

eBirdFilter <- askYesNo(msg = "Do you want to re-run the eBird Filter?")        # for re-running the eBird filter - shouldn't need to do this unless making big changes

if (eBirdFilter == TRUE) {
  source("R/eBirdFilter.R")  
}

preplots <- askYesNo(msg = "Do you want to reload from the saved state?")       # for reloading saved state (saves lots of time!) - leave if running fresh

if (preplots == TRUE) {
  load("pre-plots.RData")
} else {
  day.abb <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")                 # set up handy text vectors for later use
  day.init <- c("S", "M", "T", "W", "T", "F", "S")
  month.init <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
  # month.abb already exists as a built-in R object
  
  hab.abb <- c("Built Up", "Dry Forest", "Estuary", 
               "Non-Remnant", "Wet Forest", "Wetland")
  hab.init <- c("BU", "DF", "Es", "NR", "WF", "We")
  
  std_buffer <- 200                                                             # standard buffer width (in metres) - can change if desired
  
  sam <- read_rds("Data/eBird Filtered Dataset/Sampling.rds")                   # read in the sampling dataset
  
  sam <- sam %>%                                                                # create new columns with useful temporal information
    mutate(
      year = year(observation_date),                # y
      month = month(observation_date),              # m
      week = week(observation_date),                # w
      yday = yday(observation_date),                # d
      wday = wday(observation_date),                # d
      time = as_hms(time_observations_started),     # t
      hour = hour(time),                            # h
      minute = minute(time),                        # m
      effort_distance_km = round(replace_na(effort_distance_km, 0), digits = 2)
    )
  
  sam_sf <- sam %>%                                                             # convert sampling data to spatial format
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = 3112)
  
  bne_poly_sf <- read_sf("Data/NewSpatial/NewBrisbane.gpkg") %>%                # city boundary
    filter(brisbane == TRUE)
  
  dem_ra <- raster("Data/NewSpatial/NewBrisbaneDEM.tif") %>%                    # elevation
    projectRaster(crs = "EPSG:3112") 
  
  bne_hab_pas_sf <- read_sf("Data/NewSpatial/NewHabPA.gpkg")                    # habitat and protected areas
  
  bne_hab_pas_stat <- bne_hab_pas_sf %>%                                        # calculate areas of each
    mutate(area = as.numeric(st_area(.))) %>% 
    st_drop_geometry()
  
  sam_stdbuf <- sam_sf %>%                                                      # append elevation data to buffered sampling data
    st_buffer(dist = std_buffer) %>% 
    mutate(stdbuf_elevation = exact_extract(dem_ra, ., fun = "mean"),
           stdbuf_elevation_round = flexi_round(stdbuf_elevation, 10),
           area = as.numeric(st_area(.))) %>%
    st_drop_geometry() %>% 
    select(-area)
  
  sam_stdbuf_poly <- sam_sf %>%                                                 # append habitat and protected area data to buffered sampling data
    st_buffer(dist = std_buffer) %>% 
    st_intersection(., bne_hab_pas_sf) %>% 
    mutate(area = as.numeric(st_area(.))) %>% 
    st_drop_geometry() %>% 
    group_by(checklist_id) %>% 
    mutate(area_proportion = area / sum(area))
  
  bne_dem_vec <- extract(dem_ra, bne_poly_sf, cell = TRUE)[[1]] %>%             # get DEM for just inside Brisbane (not all SEQ)
    as_tibble() %>% 
    rename(elevation = value) %>% 
    mutate(elevation_round = flexi_round(elevation, 10)) %>% 
    group_by(elevation_round) %>% 
    summarise(p_cells = n()/nrow(.))
  
  sam_stdbuf_interaction <- left_join(sam_stdbuf_poly, sam_stdbuf)              # join two spatial for spatial interactions
  
  bne_hab_pas_stat_interaction <- bne_hab_pas_stat %>%                          # area of each habitat & PA class - inside Brisbane
    filter(brisbane == TRUE) %>% 
    mutate(id = row_number())
  
  bne_dem_interaction <- terra::extract(                                        # spatial interaction with DEM
    dem_ra, filter(bne_hab_pas_sf, brisbane == TRUE), cells = TRUE, list = TRUE)
  
  bne_dem_interaction <- bne_dem_interaction %>% 
    tibble() %>% 
    rename("values" = ".") %>% 
    mutate(id = row_number()) %>%  
    unnest(., cols = values)
  
  bne_dem_interaction <- bne_dem_interaction %>% 
    left_join(., bne_hab_pas_stat_interaction, by = "id") %>% 
    select(-area, -brisbane, -id) %>% 
    rename(elevation = values) %>% 
    mutate(elevation_round = flexi_round(elevation, 10)) %>% 
    group_by(elevation_round, habitat, protected) %>% 
    summarise(p_cells = n()/nrow(.))
  
  temporal_weights <- seq.Date(
    as.Date("2012-01-01"), as.Date("2021-12-31"), by = "days") %>% 
    expand_grid(date = ., hour = 0:23) %>% 
    mutate(year = year(date), 
           month = month(date), 
           wday = wday(date),
           weight = nrow(sam) / nrow(.)) %>% 
    select(date, year, month, wday, hour, weight)
  
  temporal_bne_dem_interaction <- expand_grid(
    bne_dem_interaction, temporal_weights) %>% 
    mutate(st_weight = p_cells * weight)
  
  ## Save State ----
  
  save.image("pre-plots.RData")                                                 # useful to save time, as the above process is very slow!
}

# Plots ----

ggs <- str_subset(ls(), "_plot")                                                # remove old plots
rm(list = ggs, ggs)

## Temporal (Single) ----

y_of_a_checklists_plot <- ggplot() +
  geom_bar(data = sam, mapping = aes(x = year), stat = "count",
           fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_weights, mapping = aes(x = year, weight = weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Year", y = "Number of Checklists") +
  scale_x_continuous(breaks = 2012:2021)

m_of_y_checklists_plot <- ggplot() +
  geom_bar(data = sam, mapping = aes(x = month), stat = "count",
           fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_weights, mapping = aes(x = month, weight = weight),
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Month of Year", y = "Number of Checklists") +
  scale_x_continuous(breaks = 1:12, labels = month.abb)

d_of_w_checklists_plot <- ggplot() +
  geom_bar(data = sam, mapping = aes(x = wday), stat = "count",
           fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_weights, mapping = aes(x = wday, weight = weight),
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Day of Week", y = "Number of Checklists") +
  scale_x_continuous(breaks = 1:7, labels = day.abb)

h_of_d_checklists_plot <- ggplot() +
  geom_bar(data = sam, mapping = aes(x = hour), stat = "count",
           fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_weights, mapping = aes(x = hour, weight = weight),
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Hour of Day", y = "Number of Checklists") +
  scale_x_continuous(breaks = seq(0, 24, 3), 
                     labels = c(paste0(seq(0, 21, 3), ":00"), "0:00"))

### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w_of_\\w_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w_of_\\w_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0("Figures/Temporal/", gg_names[[i]], ".svg"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  ggsave(filename = paste0("Figures/Temporal/", gg_names[[i]], ".png"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  
  x_vals <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(x, y)
  
  distribution_vals <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(x, y) %>% 
    rename(distribution = y)
  
  vals <- left_join(distribution_vals, x_vals, by = "x") %>% 
    mutate(y = replace_na(y, 0))
  
  label_text <- vals %>% 
    summarise(hoover_index = hoover(x = y, distribution = distribution)) %>% 
    pull() %>% 
    round(digits = 2) %>%
    paste0("Bias: ", ., " (n = ", nrow(vals), ")")
  
  plot_bias <- gg_plots[[i]] +
    labs(subtitle = label_text)
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0("Figures/Temporal/", gg_names[[i]], "_hoover.svg"), 
         plot = plot_bias, width = 5, height = 5)
  ggsave(filename = paste0("Figures/Temporal/", gg_names[[i]], "_hoover.png"), 
         plot = plot_bias, width = 5, height = 5)
}

rm(label_text, i, x_vals, distribution_vals, vals, plot_bias)

temporal_plot <- ggarrange(
  plotlist = gg_bias[c(4, 3, 1, 2)], align = "hv", ncol = 2, nrow = 2)
ggsave(filename = "Figures/Temporal/TemporalBias.svg", 
       plot = temporal_plot, width = 9, height = 7)
ggsave(filename = "Figures/Temporal/TemporalBias.png", 
       plot = temporal_plot, width = 9, height = 7)

## Temporal by Temporal ----

### Temporal by Year ----

m_of_y_by_y_checklists_plot <- m_of_y_checklists_plot +
  scale_x_continuous(breaks = 1:12, labels = month.init) +
  facet_wrap(~ year, ncol = 5)

d_of_w_by_y_checklists_plot <- d_of_w_checklists_plot +
  scale_x_continuous(breaks = 1:7, labels = day.init) +
  facet_wrap(~ year, ncol = 5)

h_of_d_by_y_checklists_plot <- h_of_d_checklists_plot +
  scale_x_continuous(breaks = seq(0, 24, 6), 
                     labels = c(paste0(seq(0, 21, 6), ":00"), "0:00")) +
  facet_wrap(~ year, ncol = 5)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w_of_\\w_by_y_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w_of_\\w_by_y_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Year/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 10, height = 5)
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Year/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 10, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(year = 2011 + as.integer(PANEL)) %>%
    select(year, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(year, "\n", "Bias: ", hoover_index, " (n = ", n, ")\n", 
                         "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- pull(label_text_df, year)
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ year, ncol = 5, labeller = labeller(year = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0("Figures/Temporal-Temporal/by Year/", 
                           gg_names[[i]], "_hoover.svg"),
         plot = plot_bias, width = 10, height = 5)
  ggsave(filename = paste0("Figures/Temporal-Temporal/by Year/", 
                           gg_names[[i]], "_hoover.png"),
         plot = plot_bias, width = 10, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

temporaltemporalyear_plot <- ggarrange(
  plotlist = gg_bias[c(3, 1, 2)], align = "hv", ncol = 1, nrow = 3)
ggsave(
  filename = "Figures/Temporal-Temporal/by Year/TemporalTemporalYearBias.svg", 
  plot = temporaltemporalyear_plot, width = 9, height = 12)
ggsave(
  filename = "Figures/Temporal-Temporal/by Year/TemporalTemporalYearBias.png", 
  plot = temporaltemporalyear_plot, width = 9, height = 12)

### Temporal by Month of Year ----

month.labs <- month.abb
names(month.labs) <- 1:12

d_of_w_by_m_checklists_plot <- d_of_w_checklists_plot +
  scale_x_continuous(breaks = 1:7, labels = day.init) +
  facet_wrap(~ month, ncol = 4, labeller = labeller(month = month.labs))

h_of_d_by_m_checklists_plot <- h_of_d_checklists_plot +
  scale_x_continuous(breaks = seq(0, 24, 6), 
                     labels = c(paste0(seq(0, 21, 6), ":00"), "0:00")) +
  facet_wrap(~ month, ncol = 4, labeller = labeller(month = month.labs))

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w_of_\\w_by_m_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w_of_\\w_by_m_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Month/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 9, height = 5)
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Month/", gg_names[[i]], ".png"),
    plot = gg_plots[[i]], width = 9, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(month = month.abb[as.integer(PANEL)]) %>%
    select(month, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(month, "\n", "Bias: ", hoover_index, 
                         " (n = ", n, ")\n", "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- 1:12
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ month, ncol = 4, labeller = labeller(month = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Month/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 9, height = 5)
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Month/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 9, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution,
   month.labs)

temporaltemporalmonth_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 2)
ggsave(
  filename = "Figures/Temporal-Temporal/by Month/TemporalTemporalMonthBias.svg", 
  plot = temporaltemporalmonth_plot, width = 9, height = 12)
ggsave(
  filename = "Figures/Temporal-Temporal/by Month/TemporalTemporalMonthBias.png", 
  plot = temporaltemporalmonth_plot, width = 9, height = 12)

### Temporal by Day of Week ----

day.labs <- day.abb
names(day.labs) <- 1:7

h_of_d_by_d_checklists_plot <- h_of_d_checklists_plot +
  scale_x_continuous(breaks = seq(0, 24, 8), 
                     labels = c(paste0(seq(0, 16, 8), ":00"), "0:00")) +
  facet_wrap(~ wday, ncol = 7, labeller = labeller(wday = day.labs))

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w_of_\\w_by_d_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w_of_\\w_by_d_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Day/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 14, height = 5)
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Day/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 14, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(wday = day.abb[as.integer(PANEL)]) %>%
    select(wday, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(wday, "\n", "Bias: ", hoover_index, " (n = ", n, ")\n", 
                         "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- 1:7
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ wday, ncol = 7, labeller = labeller(wday = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Day/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 14, height = 5)
  ggsave(filename = paste0(
    "Figures/Temporal-Temporal/by Day/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 14, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_distribution, 
   plot_bias, day.labs)

temporaltemporalday_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 1)
ggsave(
  filename = "Figures/Temporal-Temporal/by Day/TemporalTemporalDayBias.svg", 
  plot = temporaltemporalday_plot, width = 9, height = 4)
ggsave(
  filename = "Figures/Temporal-Temporal/by Day/TemporalTemporalDayBias.png", 
  plot = temporaltemporalday_plot, width = 9, height = 4)

## Spatial (Single) ----

gg_bias <- list()

### Spatial (Single) Elevation ----

elevation10_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf, mapping = aes(x = stdbuf_elevation_round), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = bne_dem_vec, fill = "red", alpha = 0.5, mapping = aes(
    x = elevation_round, y = nrow(sam_stdbuf) * p_cells), stat = "identity") +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "elevation10_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "elevation10_of_\\w+_\\w+_plot"))

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".svg"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".png"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  
  x_vals <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(x, y)
  
  distribution_vals <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(x, y) %>% 
    rename(distribution = y)
  
  vals <- left_join(distribution_vals, x_vals, by = "x") %>% 
    mutate(y = replace_na(y, 0))
  
  label_text <- vals %>% 
    summarise(hoover_index = hoover(x = y, distribution = distribution)) %>% 
    pull() %>% 
    round(digits = 2) %>%
    paste0("Bias: ", ., " (n = ", nrow(vals), ")")
  
  plot_bias <- gg_plots[[i]] +
    labs(subtitle = label_text)
  
  gg_bias <- append(gg_bias, list(plot_bias))
  
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.svg"), 
         plot = plot_bias, width = 5, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.png"), 
         plot = plot_bias, width = 5, height = 5)
}

rm(label_text, i, x_vals, distribution_vals, vals, plot_bias)

### Spatial (Single) Protected Areas ----

protectedareas_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = filter(bne_hab_pas_stat, brisbane == TRUE), mapping = aes(
    x = protected, weight = nrow(sam_sf) * area / sum(area)), 
    stat = "count", fill = "red", alpha = 0.5, ) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists")

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "protectedareas_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "protectedareas_of_\\w+_\\w+_plot"))

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".svg"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".png"), 
         plot = gg_plots[[i]], width = 5, height = 5)
  
  x_vals <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>%
    pull(y)
  
  distribution_vals <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    pull(y)
  
  label_text <- hoover(x = x_vals, distribution = distribution_vals) %>% 
    round(digits = 2) %>%
    paste0("Bias: ", ., " (n = ", length(x_vals), ")")
  
  plot_bias <- gg_plots[[i]] +
    labs(subtitle = label_text)
  
  gg_bias <- append(gg_bias, list(plot_bias))
  
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.svg"), 
         plot = plot_bias, width = 5, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.png"), 
         plot = plot_bias, width = 5, height = 5)
}

rm(label_text, i, x_vals, distribution_vals, plot_bias)

### Spatial (Single) Habitats ----

habitats_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = filter(bne_hab_pas_stat, brisbane == TRUE), mapping = aes(
    x = habitat, weight = nrow(sam_sf) * area / sum(area)), 
    stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists")

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "habitats_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "habitats_of_\\w+_\\w+_plot"))

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".svg"), 
         plot = gg_plots[[i]], width = 8, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], ".png"), 
         plot = gg_plots[[i]], width = 8, height = 5)
  
  x_vals <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>%
    pull(y)
  
  distribution_vals <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    pull(y)
  
  label_text <- hoover(x = x_vals, distribution = distribution_vals) %>% 
    round(digits = 2) %>%
    paste0("Bias: ", ., " (n = ", length(x_vals), ")")
  
  plot_bias <- gg_plots[[i]] +
    labs(subtitle = label_text)
  
  gg_bias <- append(gg_bias, list(plot_bias))
  
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.svg"), 
         plot = plot_bias, width = 8, height = 5)
  ggsave(filename = paste0("Figures/Spatial/", gg_names[[i]], "_hoover.png"), 
         plot = plot_bias, width = 8, height = 5)
}

rm(label_text, i, x_vals, distribution_vals, plot_bias)

### Save ----

spatial_plot <- ggarrange(
  ggarrange(plotlist = gg_bias[c(1, 2)], align = "hv", ncol = 2, nrow = 1),
  ggarrange(plotlist = gg_bias[c(3)], ncol = 1, nrow = 1),
  ncol = 1, nrow = 2)
ggsave(filename = "Figures/Spatial/SpatialBias.svg", 
       plot = spatial_plot, width = 9, height = 9)
ggsave(filename = "Figures/Spatial/SpatialBias.png", 
       plot = spatial_plot, width = 9, height = 9)

## Spatial by Spatial ----

### Spatial by Habitats ----

elevation10_by_habitats_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_interaction, 
           mapping = aes(x = stdbuf_elevation_round, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = bne_dem_interaction, 
           mapping = aes(x = elevation_round, 
                         weight = nrow(sam_stdbuf) * p_cells), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ habitat, ncol = 3)

protectedareas_by_habitats_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_interaction, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = bne_hab_pas_stat_interaction, 
           mapping = aes(x = protected, 
                         weight = nrow(sam_sf) * area / sum(area)), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists") +
  facet_wrap(~ habitat, ncol = 3)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w+_by_habitats_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w+_by_habitats_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Habitat/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 8, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Habitat/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 8, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(habitat = hab.abb[as.integer(PANEL)]) %>%
    select(habitat, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(habitat, "\n", "Bias: ", hoover_index, 
                         " (n = ", n, ")\n", "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- hab.abb
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ habitat, ncol = 3, labeller = labeller(habitat = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Habitat/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 8, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Habitat/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 8, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialspatialhabitat_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 2)
ggsave(filename = 
         "Figures/Spatial-Spatial/by Habitat/SpatialSpatialHabitatBias.svg", 
       plot = spatialspatialhabitat_plot, width = 9, height = 12)
ggsave(filename = 
         "Figures/Spatial-Spatial/by Habitat/SpatialSpatialHabitatBias.png", 
       plot = spatialspatialhabitat_plot, width = 9, height = 12)

### Spatial by Protected Areas ----

elevation10_by_protectedareas_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_interaction, 
           mapping = aes(x = stdbuf_elevation_round, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = bne_dem_interaction, 
           mapping = aes(x = elevation_round, y = nrow(sam_stdbuf) * p_cells), 
           stat = "identity", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ protected, ncol = 2)

habitats_by_protectedareas_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_interaction, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = bne_hab_pas_stat_interaction, 
           mapping = aes(x = habitat, weight = nrow(sam_sf) * area / sum(area)), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists") +
  scale_x_discrete(labels = hab.init) +
  facet_wrap(~ protected, ncol = 2)

#### Save ----

gg_plots <- as.list(mget(str_subset(
  ls(), "\\w+_by_protectedareas_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(
  ls(), "\\w+_by_protectedareas_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Protected Area/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 8, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Protected Area/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 8, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(protected = c(FALSE, TRUE)[as.integer(PANEL)]) %>%
    select(protected, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(protected, "\n", "Bias: ", hoover_index, 
                         " (n = ", n, ")\n", "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- c(FALSE, TRUE)
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ protected, ncol = 4, 
               labeller = labeller(protected = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Protected Area/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 8, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Spatial/by Protected Area/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 8, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialspatialprotectedarea_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 2)
ggsave(filename = paste0("Figures/Spatial-Spatial/by Protected Area/", 
                         "SpatialSpatialProtectedAreaBias.svg"), 
       plot = spatialspatialprotectedarea_plot, width = 9, height = 7)
ggsave(filename = paste0("Figures/Spatial-Spatial/by Protected Area/", 
                         "SpatialSpatialProtectedAreaBias.png"),
       plot = spatialspatialprotectedarea_plot, width = 9, height = 7)

## Spatial by Temporal ----

### Spatial by Year ----

elevation10_by_y_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf, 
           mapping = aes(x = stdbuf_elevation_round), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = elevation_round, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ year, ncol = 5)

protectedareas_by_y_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = protected, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists") +
  facet_wrap(~ year, ncol = 5)

habitats_by_y_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = habitat, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists") +
  scale_x_discrete(labels = hab.init) +
  facet_wrap(~ year, ncol = 5)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w+_by_y_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w+_by_y_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Year/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 12, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Year/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 12, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(year = 2011 + as.integer(PANEL)) %>%
    select(year, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(year, "\n", "Bias: ", hoover_index, " (n = ", n, ")\n", 
                         "Adjustment: ", hoover_adjust))
  
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- pull(label_text_df, year)
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ year, ncol = 5, labeller = labeller(year = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Year/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 12, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Year/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 12, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialtemporalyear_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 3)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Year/SpatialTemporalYearBias.svg", 
       plot = spatialtemporalyear_plot, width = 9, height = 12)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Year/SpatialTemporalYearBias.png", 
       plot = spatialtemporalyear_plot, width = 9, height = 12)

### Spatial by Month of Year ----

elevation10_by_m_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf, 
           mapping = aes(x = stdbuf_elevation_round), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = elevation_round, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ month, ncol = 4)

protectedareas_by_m_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = protected, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists") +
  facet_wrap(~ month, ncol = 4)

habitats_by_m_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = habitat, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists") +
  scale_x_discrete(labels = hab.init) +
  facet_wrap(~ month, ncol = 4)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w+_by_m_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w+_by_m_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Month/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 10, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Month/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 10, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(month = month.abb[as.integer(PANEL)]) %>%
    select(month, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(month, "\n", "Bias: ", hoover_index, 
                         " (n = ", n, ")\n", "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- 1:12
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ month, ncol = 4, labeller = labeller(month = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Month/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 10, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Month/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 10, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialtemporalmonth_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 3)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Month/SpatialTemporalMonthBias.svg", 
       plot = spatialtemporalmonth_plot, width = 9, height = 12)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Month/SpatialTemporalMonthBias.png", 
       plot = spatialtemporalmonth_plot, width = 9, height = 12)

### Spatial by Day of Week ----

elevation10_by_d_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf, 
           mapping = aes(x = stdbuf_elevation_round), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = elevation_round, weight = st_weight),
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ wday, ncol = 7)

protectedareas_by_d_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = protected, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists") +
  facet_wrap(~ wday, ncol = 7)

habitats_by_d_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = habitat, weight = st_weight),
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists") +
  scale_x_discrete(labels = hab.init) +
  facet_wrap(~ wday, ncol = 7)

#### Save ----

gg_plots <- as.list(mget(str_subset(ls(), "\\w+_by_d_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w+_by_d_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Day/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 14, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Day/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 14, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(wday = day.abb[as.integer(PANEL)]) %>%
    select(wday, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(wday, "\n", "Bias: ", hoover_index, " (n = ", n, ")\n", 
                         "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- 1:7
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ wday, ncol = 7, labeller = labeller(wday = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Day/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 14, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Day/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 14, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialtemporalday_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 3)
ggsave(filename = "Figures/Spatial-Temporal/by Day/SpatialTemporalDayBias.svg", 
       plot = spatialtemporalday_plot, width = 9, height = 12)
ggsave(filename = "Figures/Spatial-Temporal/by Day/SpatialTemporalDayBias.png", 
       plot = spatialtemporalday_plot, width = 9, height = 12)

### Spatial by Hour of Day ----

elevation10_by_h_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf, 
           mapping = aes(x = stdbuf_elevation_round), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = elevation_round, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Number of Checklists") + 
  scale_x_continuous(limits = c(-10, 710), oob = function(x, limits) x) +
  facet_wrap(~ hour, ncol = 6)

protectedareas_by_h_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = protected, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = protected, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Protected Area Status", y = "Number of Checklists") +
  facet_wrap(~ hour, ncol = 6)

habitats_by_h_of_stdbuf_checklists_plot <- ggplot() +
  geom_bar(data = sam_stdbuf_poly, 
           mapping = aes(x = habitat, weight = area_proportion), 
           stat = "count", fill = "blue", alpha = 0.5) +
  geom_bar(data = temporal_bne_dem_interaction, 
           mapping = aes(x = habitat, weight = st_weight), 
           stat = "count", fill = "red", alpha = 0.5) +
  theme_classic() +
  labs(x = "Habitat Type", y = "Number of Checklists") +
  scale_x_discrete(labels = hab.init) +
  facet_wrap(~ hour, ncol = 6)

#### Save ----

hour.labs <- c(paste0(seq(0, 23, 1), ":00"))
names(hour.labs) <- 0:23

gg_plots <- as.list(mget(str_subset(ls(), "\\w+_by_h_of_\\w+_\\w+_plot")))
gg_names <- as.list(str_subset(ls(), "\\w+_by_h_of_\\w+_\\w+_plot"))

gg_bias <- list()

for (i in 1:length(gg_plots)) {
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Hour/", gg_names[[i]], ".svg"), 
    plot = gg_plots[[i]], width = 12, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Hour/", gg_names[[i]], ".png"), 
    plot = gg_plots[[i]], width = 12, height = 5)
  
  plot_x <- gg_plots[[i]] %>% 
    layer_data(i = 1) %>% 
    select(PANEL, x, y)
  
  plot_distribution <- gg_plots[[i]] %>% 
    layer_data(i = 2) %>% 
    select(PANEL, x, y) %>% 
    rename(distribution = y)
  
  plot_all <- full_join(plot_distribution, plot_x, by = c("PANEL", "x")) %>% 
    mutate(y = replace_na(y, 0),
           distribution = replace_na(distribution, 0),
           PANEL = as.integer(PANEL))
  
  label_text_df <- plot_all %>% 
    group_by(PANEL) %>%
    summarise(hoover_index = hoover(y, distribution)[1],
              hoover_adjust = hoover(y, distribution)[2],
              n = n_distinct(x)) %>% 
    mutate(hoover_index = round(hoover_index, digits = 2),
           hoover_adjust = round(hoover_adjust, digits = 2)) %>%
    mutate(hour = hour.labs[as.integer(PANEL)]) %>%
    select(hour, hoover_index, hoover_adjust, n) %>%
    mutate(text = paste0(hour, "\n", "Bias: ", hoover_index, " (n = ", n, ")\n", 
                         "Adjustment: ", hoover_adjust))
  
  label_text <- pull(label_text_df, text)
  names(label_text) <- 1:24
  
  plot_bias <- gg_plots[[i]] +
    facet_wrap(~ hour, ncol = 6, labeller = labeller(hour = label_text))
  
  gg_bias[[i]] <- plot_bias
  
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Hour/", gg_names[[i]], "_hoover.svg"),
    plot = plot_bias, width = 12, height = 5)
  ggsave(filename = paste0(
    "Figures/Spatial-Temporal/by Hour/", gg_names[[i]], "_hoover.png"),
    plot = plot_bias, width = 12, height = 5)
}

rm(label_text, label_text_df, i, plot_x, plot_all, plot_bias, plot_distribution)

spatialtemporalhour_plot <- ggarrange(
  plotlist = gg_bias, align = "hv", ncol = 1, nrow = 3)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Hour/SpatialTemporalHourBias.svg", 
       plot = spatialtemporalhour_plot, width = 9, height = 12)
ggsave(filename = 
         "Figures/Spatial-Temporal/by Hour/SpatialTemporalHourBias.png", 
       plot = spatialtemporalhour_plot, width = 9, height = 12)

## Save State ----

save(list = str_subset(ls(), "_plot"), file = "plots.RData")
