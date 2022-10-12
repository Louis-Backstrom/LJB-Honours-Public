# Libraries and Functions ----

library(tidyverse)
library(ggpubr)

source("R/CustomFunctions.R")

# Set-up ----

all_habitats <- c("red", "green", "blue")                                       # set up 3 habitats
all_time <- 1:3650                                                              # 10 years x 365 days
all_prob <- seq(from = 100/73, by = 40/5329, length.out = 3650)                 # temporal bias - linear: 1000 in year 1, 10000 in year 10

scenarios <- tibble(                                                            # set up 6 scenarios
  start_occupancy = list(c(0.7, 0.3, 0.3), c(0.7, 0.3, 0.3),                    # starting occupancy per-habitat per-scenario  
                         c(0.7, 0.3, 0.3), c(0.7, 0.3, 0.3), 
                         c(0.7, 0.3, 0.3), c(0.7, 0.3, 0.3)),
  delta_occupancy = list(c(0, 0, 0), c(0, -0.001, -0.001),                      # daily change in occupancy per-habitat per-scenario
                         c(0, 0, 0), c(0, 0, 0),
                         c(0, -0.001, -0.001),  c(0, -0.001, -0.001)),
  start_habbiases = list(c(1/3, 1/3, 1/3), c(1/3, 1/3, 1/3),                    # starting bias per-habitat
                         c(1/9, 3/9, 5/9), c(1/9, 3/9, 5/9),
                         c(1/9, 3/9, 5/9), c(1/9, 3/9, 5/9)),
  delta_habbiases = list(c(0, 0, 0), c(0, 0, 0),                                # total change in bias per-habitat
                         c(0, 0, 0), c(4/9, 0, -4/9),
                         c(0, 0, 0), c(4/9, 0, -4/9))
)

# Run Scenarios ----

scenario_plots <- vector(mode = "list", length = 6)

for (i in 1:nrow(scenarios)) {                                                  # run each scenario
  start_occupancy <- unlist(scenarios[i, 1])                                    # set up the starting occupancy for the scenario
  delta_occupancy <- unlist(scenarios[i, 2])                                    # set up the daily change in occupancy for the scenario
  start_habbiases <- unlist(scenarios[i, 3])                                    # set up the starting bias for the scenario
  delta_habbiases <- unlist(scenarios[i, 4])                                    # set up the total change in bias for the scenario
  
  r_bias <- seq(from = start_habbiases[1],                                      # red bias over time
                by = delta_habbiases[1]/3650, 
                length.out = 3650)
  g_bias <- seq(from = start_habbiases[2],                                      # green bias over time
                by = delta_habbiases[2]/3650, 
                length.out = 3650)
  b_bias <- seq(from = start_habbiases[3],                                      # blue bias over time
                by = delta_habbiases[3]/3650, 
                length.out = 3650)
  bias <- tibble(r_bias, g_bias, b_bias) %>%                                    # all biases over time
    rename(red = r_bias, green = g_bias, blue = b_bias)
  
  r_occupancy <- occ_change(start = start_occupancy[1],                         # red occupancy over time
                            change = delta_occupancy[1], 
                            length.out = 3650)
  g_occupancy <- occ_change(start = start_occupancy[2],                         # green occupancy over time
                            change = delta_occupancy[2], 
                            length.out = 3650)
  b_occupancy <- occ_change(start = start_occupancy[3],                         # blue occupancy over time
                            change = delta_occupancy[3], 
                            length.out = 3650)
  occupancy <- tibble(r_occupancy, g_occupancy, b_occupancy) %>%                # all occupancies over time
    rename(red = r_occupancy, green = g_occupancy, blue = b_occupancy)
  
  total_occupancy <- rowMeans(occupancy)                                        # total (average) occupancy over time
  
  simulation_occupancy <- occupancy %>%                                         # collate occupancies into format for simulations
    rowid_to_column(var = "time") %>% 
    pivot_longer(cols = c("red", "green", "blue"), 
                 names_to = "habitat", values_to = "occupancy")
  
  for (j in 1:100) {                                                            # run 100 simulations per scenario
    simdata <- tibble(                                                          # set up a tibble of simulated data
      time = sample(all_time, size = 55000, replace = TRUE, prob = all_prob),   # randomly sample the time period as per temporal bias
      r_bias = r_bias[time],                                                    # determine the red bias at each time
      g_bias = g_bias[time],                                                    # determine the green bias at each time
      b_bias = b_bias[time],                                                    # determine the blue bias at each time
      r_occupancy = r_occupancy[time],                                          # determine the red occupancy at each time
      g_occupancy = g_occupancy[time],                                          # determine the green occupancy at each time
      b_occupancy = b_occupancy[time]) %>%                                      # determine the blue occupancy at each time
      rowwise() %>% 
      mutate(habitat = sample(                                                  # randomly sample the habitat for each event as per spatial bias
        all_habitats, size = 1, prob = c(r_bias, g_bias, b_bias))) %>% 
      left_join(., simulation_occupancy, by = c("time", "habitat")) %>%         # append the occupancy values to each event
      mutate(present = sample(                                                  # randomly sample detection based on occupancy
        c(1, 0), size = 1, prob = c(occupancy, 1 - occupancy))) %>% 
      ungroup()  %>% 
      select(time, habitat, present)
    
    simdata_train <- simdata[seq(1, 55000, by = 2), ]                           # select odd rows for training the models
    
    nullinformed_glm <- glm(data = simdata_train,                               # run a basic, uninformed GLM
                            formula = present ~ time,
                            family = "binomial")
    semiinformed_glm <- glm(data = simdata_train,                               # run a simple, semi-informed GLM
                            formula = present ~ time + habitat,
                            family = "binomial")
    fullinformed_glm <- glm(data = simdata_train,                               # run a complex, fully-informed GLM
                            formula = present ~ time * habitat,
                            family = "binomial")
    
    simdata_tests <- simdata[seq(2, 55000, by = 2), ]                           # select even rows for testing the models
    
    simdata_tests$nullinformed_prediction <- nullinformed_glm %>%               # run uninformed model predictions on test data
      predict(., simdata_tests) %>% 
      plogis()
    simdata_tests$semiinformed_prediction <- semiinformed_glm %>%               # run semi-informed model predictions on test data
      predict(., simdata_tests) %>% 
      plogis()
    simdata_tests$fullinformed_prediction <- fullinformed_glm %>%               # run fully-informed model predictions on test data
      predict(., simdata_tests) %>% 
      plogis()
    
    simdata_tests <- simdata_tests %>% 
      mutate(nullinformed_int_error = present - nullinformed_prediction,        # assess model errors relative to observed values (0 or 1)
             semiinformed_int_error = present - semiinformed_prediction,
             fullinformed_int_error = present - fullinformed_prediction) %>% 
      left_join(., simulation_occupancy, by = c("time", "habitat")) %>%         # append actual values
      mutate(nullinformed_ext_error = occupancy - nullinformed_prediction,      # assess model errors relative to actual values (occupancy probability)
             semiinformed_ext_error = occupancy - semiinformed_prediction,
             fullinformed_ext_error = occupancy - fullinformed_prediction)
    
    run_performance = tibble(
      run = j,
      nullinformed_int_rmse = rmse(simdata_tests$nullinformed_int_error),
      semiinformed_int_rmse = rmse(simdata_tests$semiinformed_int_error),
      fullinformed_int_rmse = rmse(simdata_tests$fullinformed_int_error),
      nullinformed_ext_rmse = rmse(simdata_tests$nullinformed_ext_error),
      semiinformed_ext_rmse = rmse(simdata_tests$semiinformed_ext_error),
      fullinformed_ext_rmse = rmse(simdata_tests$fullinformed_ext_error),
      nullinformed_int_rss = rss(simdata_tests$nullinformed_int_error),
      semiinformed_int_rss = rss(simdata_tests$semiinformed_int_error),
      fullinformed_int_rss = rss(simdata_tests$fullinformed_int_error),
      nullinformed_ext_rss = rss(simdata_tests$nullinformed_ext_error),
      semiinformed_ext_rss = rss(simdata_tests$semiinformed_ext_error),
      fullinformed_ext_rss = rss(simdata_tests$fullinformed_ext_error)
    )
    
    simdata_plots <- expand_grid(time = all_time, habitat = all_habitats)       # set up a tibble to generate plot data
    
    simdata_plots$nullinformed_prediction <- nullinformed_glm %>%               # run uninformed model predictions on plot data
      predict(., simdata_plots) %>% 
      plogis()
    simdata_plots$semiinformed_prediction <- semiinformed_glm %>%               # run semi-informed model predictions on plot data
      predict(., simdata_plots) %>% 
      plogis()
    simdata_plots$fullinformed_prediction <- fullinformed_glm %>%               # run fully-informed model predictions on plot data
      predict(., simdata_plots) %>% 
      plogis()
    
    simdata_plots <- simdata_plots %>%                                          # append actual values
      left_join(., simulation_occupancy, by = c("time", "habitat"))
    
    simdata_plots <- simdata_plots %>%                                          # group by time and find average values
      group_by(time) %>% 
      summarise(nullinformed_prediction = mean(nullinformed_prediction),
                semiinformed_prediction = mean(semiinformed_prediction),
                fullinformed_prediction = mean(fullinformed_prediction),
                occupancy = mean(occupancy)) %>% 
      mutate(habitat = "average") %>% 
      rbind(simdata_plots, .) %>%                                               # join back so that there are 4 values per time (r/g/b & average)
      arrange(time, habitat) %>% 
      mutate(run = j)
    
    if (j == 1) {                                                               # save data from each run into a tibble
      scenario_performance = run_performance
      scenario_plotdata = simdata_plots
    } else {
      scenario_performance = rbind(scenario_performance, run_performance)
      scenario_plotdata = rbind(scenario_plotdata, simdata_plots)
    }
    
    if (j %% 10 == 0) {
      print(paste0("scenario ", i, ", run ", j))  
    }
    
  }
  
  plot_data <- scenario_plotdata %>%                                            # set up the dataframe used in the plots
    group_by(habitat, time) %>% 
    summarise(nullinformed_qmin = min(nullinformed_prediction),                 # get min, max, and 10, 25, 75, 90% quantiles of model runs
              nullinformed_q010 = quantile(nullinformed_prediction, 0.1),
              nullinformed_q025 = quantile(nullinformed_prediction, 0.25),
              nullinformed_q075 = quantile(nullinformed_prediction, 0.75),
              nullinformed_q090 = quantile(nullinformed_prediction, 0.90),
              nullinformed_qmax = max(nullinformed_prediction),
              semiinformed_qmin = min(semiinformed_prediction),
              semiinformed_q010 = quantile(semiinformed_prediction, 0.1),
              semiinformed_q025 = quantile(semiinformed_prediction, 0.25),
              semiinformed_q075 = quantile(semiinformed_prediction, 0.75),
              semiinformed_q090 = quantile(semiinformed_prediction, 0.90),
              semiinformed_qmax = max(semiinformed_prediction),
              fullinformed_qmin = min(fullinformed_prediction),
              fullinformed_q010 = quantile(fullinformed_prediction, 0.1),
              fullinformed_q025 = quantile(fullinformed_prediction, 0.25),
              fullinformed_q075 = quantile(fullinformed_prediction, 0.75),
              fullinformed_q090 = quantile(fullinformed_prediction, 0.90),
              fullinformed_qmax = max(fullinformed_prediction))
  
  simulation_occupancy <- simulation_occupancy %>% 
    group_by(time) %>% 
    summarise(occupancy = mean(occupancy)) %>% 
    mutate(habitat = "average") %>% 
    rbind(simulation_occupancy, .) %>% 
    arrange(time, habitat)
  
  plot_base <- ggplot(data = plot_data, aes(x = time, fill = habitat)) +        # set up the basic (shared) elements of the plots
    theme_classic() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("black", "blue", "green", "red")) +
    scale_x_continuous(breaks = seq(0, 3650, by = 365), labels = 0:10) +
    lims(y = c(0, 1)) +
    labs(x = "", y = "")
  
  plot_nullinformed <- plot_base +                                              # plot un-informed model
    geom_ribbon(alpha = 0.1,
                aes(ymin = nullinformed_qmin, ymax = nullinformed_qmax)) +
    geom_ribbon(alpha = 0.2, 
                aes(ymin = nullinformed_q010, ymax = nullinformed_q090)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = nullinformed_q025, ymax = nullinformed_q075))
  
  plot_semiinformed <- plot_base +                                              # plot semi-informed model
    geom_ribbon(alpha = 0.1,
                aes(ymin = semiinformed_qmin, ymax = semiinformed_qmax)) +
    geom_ribbon(alpha = 0.2, 
                aes(ymin = semiinformed_q010, ymax = semiinformed_q090)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = semiinformed_q025, ymax = semiinformed_q075))
  
  plot_fullinformed <- plot_base +                                              # plot fully-informed model
    geom_ribbon(alpha = 0.1,
                aes(ymin = fullinformed_qmin, ymax = fullinformed_qmax)) +
    geom_ribbon(alpha = 0.2, 
                aes(ymin = fullinformed_q010, ymax = fullinformed_q090)) +
    geom_ribbon(alpha = 0.3,
                aes(ymin = fullinformed_q025, ymax = fullinformed_q075))
  
  plot_truth <- ggplot(data = simulation_occupancy) +
    geom_line(aes(x = time, y = occupancy, colour = habitat), 
              size = 1, alpha = 0.5) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_colour_manual(values = c("black", "blue", "green", "red")) +
    scale_x_continuous(breaks = seq(0, 3650, by = 365), labels = 0:10) +
    lims(y = c(0, 1)) +
    labs(x = "", y = "")
  
  plot_combined <- ggarrange(                                                   # arrange plots
    plot_nullinformed, plot_semiinformed, plot_fullinformed, plot_truth,
    ncol = 4)                                              

  scenario_plots[[i]] <- plot_combined                                          # save plots to list
  
  scenario_performance <- scenario_performance %>%                              
    mutate(scenario = i)                                                        # add scenario number
  
  if (i == 1) {                                                                 # save model performance across all scenarios to large tibble
    scenario_summary <- scenario_performance  
  } else {
    scenario_summary = rbind(scenario_summary, scenario_performance)
  }
  
  print(paste0("scenario ", i, " - done!"))
}

# Outputs ----

full_plot <- ggarrange(plotlist = scenario_plots, nrow = 6)                     # create full plot
ggsave(filename = "Figures/SimulationStudy_ModelOutput.svg",                    # save full plot
       plot = full_plot, width = 8, height = 11)

summary <- scenario_summary %>%                                                 
  pivot_longer(cols = 2:13) %>%                                                 # pivot longer, to allow for easier summarising
  separate(name, into = c("model", "type", "error"), sep = "_")

plot_cor_rmse <- ggplot(                                                        # plot correlations between internal and external RMSE errors
  data = filter(pivot_wider(summary, names_from = "type"), error == "rmse"),
  aes(x = int, y = ext)) +
  geom_point(alpha = 0.5) +
  stat_cor(p.accuracy	= 0.001, r.accuracy = 0.001, 
           geom = "label", fill = "white") +
  facet_wrap(scenario ~ model, scales = "free", nrow = 6, ncol = 3) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Internal RMSE\n(Training-Test)",
       y = "External RMSE\n(Truth-Test)")

ggsave(filename = "Figures/SimulationStudy_CorrelationsRMSE.svg",
       plot = plot_cor_rmse, width = 9, height = 12)

plot_cor_rss <- ggplot(                                                         # plot correlations between internal and external RSS errors
  data = filter(pivot_wider(summary, names_from = "type"), error == "rss"),
  aes(x = int, y = ext)) +
  geom_point(alpha = 0.5) +
  stat_cor(p.accuracy	= 0.001, r.accuracy = 0.001, 
           geom = "label", fill = "white") +
  facet_wrap(scenario ~ model, scales = "free", nrow = 6, ncol = 3) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x = "Internal RSS\n(Training-Test)",
       y = "External RSS\n(Truth-Test)")

ggsave(filename = "Figures/SimulationStudy_CorrelationsRSS.svg",
       plot = plot_cor_rss, width = 9, height = 12)

simulation_performance <- summary %>%                                           # summarise model performance across runs
  group_by(scenario, model, type, error) %>% 
  summarise(mean = mean(value),
            sd = sqrt(var(value))) %>% 
  mutate(model_order = if_else(model == "nullinformed", 1,
                               if_else(model == "semiinformed", 2, 3)),
         type_order = if_else(type == "int", 1, 2)) %>% 
  arrange(error, type_order, scenario, model_order) %>% 
  select(error, type, scenario, model, mean, sd)

write_csv(simulation_performance, "Outputs/SimulationPerformance.csv")
