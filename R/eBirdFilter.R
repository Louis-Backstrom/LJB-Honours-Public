# Libraries ----

library(auk)
library(tidyverse)
library(sf)

# Spatial Boundary ----

bne <- read_sf("Data/NewSpatial/NewBrisbane.gpkg") %>% 
  st_transform(crs = 4326)

# Prepare Observation Data ----

observation_data <- auk_ebd(
  file = paste0("Data/eBird Basic Dataset January 2022 - Brisbane/", 
                "ebd_obs_2022-01-17_BNE.txt"),
  file_sampling = paste0("Data/eBird Basic Dataset January 2022 - Brisbane/",
                         "ebd_sampling_2022-01-17_BNE.txt")) %>% 
  auk_complete() %>%                                                            # only complete checklists
  auk_county("AU-QLD-BRI") %>%                                                  # only Brisbane
  auk_distance(distance = c(0, 8)) %>%                                          # checklists 0 (stationary) to 8km long
  auk_duration(duration = c(5, 300)) %>%                                        # checklists 5min to 5hr long
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%                     # standard protocols
  auk_year(year = 2012:2021) %>% 
  auk_filter(file = "Data/eBird Filtered Dataset/ebd_bne.txt",
             file_sampling = "Data/eBird Filtered Dataset/sam_bne.txt",
             filter_sampling = TRUE,
             overwrite = TRUE) %>% 
  read_ebd(rollup = TRUE, unique = TRUE)

sampling_data <- read_sampling("Data/eBird Filtered Dataset/sam_bne.txt",
                               unique = TRUE)

# Filter to Mainland Brisbane ----

mainland_checklists <- sampling_data %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_join(., bne) %>% 
  filter(brisbane == TRUE) %>% 
  st_drop_geometry() %>% 
  pull(checklist_id)

observation_data <- observation_data %>% 
  filter(checklist_id %in% mainland_checklists)

sampling_data <- sampling_data %>% 
  filter(checklist_id %in% mainland_checklists)

# Output ----

write_rds(sampling_data, "Data/eBird Filtered Dataset/Sampling.rds")    
write_rds(observation_data, "Data/eBird Filtered Dataset/Observation.rds")

rm(bne, observation_data, sampling_data, mainland_checklists)
