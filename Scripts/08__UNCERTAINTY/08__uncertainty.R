# ---
# title: "08__uncertainty.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-03_09"
# updated: "2025-12_30"
# ---


# WORK FLOW NOTES:

# WHAT
# - Visualises uncertainty within binary projections and predictions for each spp, time, scenario, threshold


#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#


### Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")


library(terra)
library(ggplot2)
library(viridis)
library(patchwork)
library(stringr)
library(dplyr)

#==============================================================================#
#                        ---- 1. Load in required data ----
#==============================================================================#


## Set up base folders ----
datafolder <- file.path("//xxx")

output_dir <- file.path(datafolder, "08__uncertainty")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Species occurrences
species_list <- c("calo", "vertlept")

# SSP126, SSP370, SSP585
climate_pathways <- c("ssp126", "ssp370", "ssp585")

# Create a nested list to store thresholds by species and threshold type
threshold_values <- list()
threshold_types <- c("opt2", "opt4", "opt12")

# Define paths
datafolder05_curr <- file.path(datafolder, "05__model_prediction", "output", 
                               "current_projections", "weighted_mean_ensembles", 
                               "extra_settings")
projections <- list.files(datafolder05_curr, pattern="current_prediction", 
                          full.names = TRUE)

# Load the two current projection rasters
calo_proj <- rast(projections[1])
vertlept_proj <- rast(projections[2])

# Inputs: Predictions for future climates
prediction_list <- list.files(file.path(datafolder, "05__model_prediction", "output", 
                                        "future_predictions", "extra_settings"), 
                              pattern="prediction", full.names = TRUE)

# Load threshold values
thresholds_dir <- file.path(datafolder, "05__model_prediction", "intermediate", "thresholds")
threshold_files <- list.files(thresholds_dir, full.names = TRUE)


#==============================================================================#
#                           2. Uncertainty function ----
#==============================================================================#

# all encompassing function that runs the uncertainty analysis, plots and saves
# Uncertainty visualisation function
uncertainty_viz <- function(prediction_raster, 
                            uncertainty_raster, 
                            threshold, 
                            threshold_type = "",
                            species_name = "",
                            time_period = "",
                            scenario = "",
                            uncertainty_percentile = 75,
                            plot_title = NULL,
                            border_sf = NULL) {
  
  # Create auto-generated title if not provided
  if (is.null(plot_title)) {
    species_display <- paste0(toupper(substr(species_name, 1, 1)), 
                              substr(species_name, 2, nchar(species_name)))
    
    if (time_period == "current") {
      time_scenario_text <- "Current Climate"
    } else {
      # Format scenario for display (e.g., "SSP1-2.6" instead of "ssp126")
      scenario_display <- toupper(substr(scenario, 1, 3)) # Extract "SSP"
      
      # Extract the numbers and format them
      numbers <- substr(scenario, 4, nchar(scenario))
      if (nchar(numbers) == 3) {
        # For scenarios like ssp126, ssp370
        scenario_display <- paste0(scenario_display, numbers[1], "-", substr(numbers, 2, 3))
      } else {
        scenario_display <- paste0(scenario_display, numbers)
      }
      
      time_scenario_text <- paste0(time_period, " (", scenario_display, ")")
    }
    
    # Create the full title
    plot_title <- paste0(species_display, " Distribution: ", time_scenario_text,
                         "\nBinary threshold: ", round(threshold, 3),
                         " (", threshold_type, ")")
  }
  
  # Convert rasters to data frames for ggplot
  pred_df <- as.data.frame(prediction_raster, xy = TRUE)
  names(pred_df)[3] <- "prediction"
  
  uncert_df <- as.data.frame(uncertainty_raster, xy = TRUE)
  names(uncert_df)[3] <- "uncertainty"
  
  # Combine data
  viz_df <- merge(pred_df, uncert_df)
  
  # Create binary classification
  viz_df$presence <- ifelse(viz_df$prediction > threshold, 1, 0)
  
  # Calculate uncertainty threshold based on percentile
  uncertainty_threshold <- stats::quantile(viz_df$uncertainty, uncertainty_percentile/100, na.rm = TRUE)
  
  # Check if the threshold is 0, which can happen with sparse uncertainty data
  if (uncertainty_threshold == 0) {
    # Find the first non-zero percentile
    for (test_percentile in seq(uncertainty_percentile + 5, 99, by = 5)) {
      test_threshold <- stats::quantile(viz_df$uncertainty, test_percentile/100, na.rm = TRUE)
      if (test_threshold > 0) {
        message(paste0("Original ", uncertainty_percentile, "th percentile gave threshold of 0. ",
                       "Using ", test_percentile, "th percentile instead."))
        uncertainty_percentile <- test_percentile
        uncertainty_threshold <- test_threshold
        break
      }
    }
    
    # If all percentiles give 0, use the mean or a small positive value
    if (uncertainty_threshold == 0) {
      mean_uncertainty <- mean(viz_df$uncertainty, na.rm = TRUE)
      if (mean_uncertainty > 0) {
        message(paste0("All percentiles gave threshold of 0. Using mean uncertainty: ", 
                       round(mean_uncertainty, 3)))
        uncertainty_threshold <- mean_uncertainty
      } else {
        # Last resort - use a small positive value
        uncertainty_threshold <- 0.05
        message(paste0("All percentiles and mean gave threshold of 0. Using default value: ", 
                       uncertainty_threshold))
      }
    }
  }
  
  # Print the calculated threshold
  message(paste0("Using ", uncertainty_percentile, "th percentile uncertainty threshold: ", 
                 round(uncertainty_threshold, 3)))
  
  # Calculate confidence percentage (for legend labeling)
  confidence_percent <- round(100 - (uncertainty_percentile))
  
  # Create confidence categories using both prediction and uncertainty
  viz_df$category <- dplyr::case_when(
    viz_df$presence == 1 & viz_df$uncertainty < uncertainty_threshold ~ 
      paste0("High confidence (", confidence_percent, "%+) presence"),
    viz_df$presence == 1 & viz_df$uncertainty >= uncertainty_threshold ~ 
      paste0("Low confidence (<", confidence_percent, "%) presence"),
    viz_df$presence == 0 & viz_df$uncertainty < uncertainty_threshold ~ 
      paste0("High confidence (", confidence_percent, "%+) absence"),
    viz_df$presence == 0 & viz_df$uncertainty >= uncertainty_threshold ~ 
      paste0("Low confidence (<", confidence_percent, "%) absence")
  )
  
  # Create category level definitions for consistent reference
  high_conf_presence <- paste0("High confidence (", confidence_percent, "%+) presence")
  low_conf_presence <- paste0("Low confidence (<", confidence_percent, "%) presence")
  high_conf_absence <- paste0("High confidence (", confidence_percent, "%+) absence")
  low_conf_absence <- paste0("Low confidence (<", confidence_percent, "%) absence")
  
  # Order factors for legend
  viz_df$category <- factor(viz_df$category, 
                            levels = c(high_conf_presence, 
                                       low_conf_presence,
                                       low_conf_absence, 
                                       high_conf_absence))
  
  # Define custom fill colours (colour-blind friendly)
  category_colours <- c()
  category_colours[high_conf_presence] <- "#1b9e77"
  category_colours[low_conf_presence] <- "#7570b3"
  category_colours[low_conf_absence] <- "#d95f02"
  category_colours[high_conf_absence] <- "#e7e7e7"
  
  # Create subtitle with uncertainty threshold info
  subtitle_text <- paste0("Uncertainty threshold: ", round(uncertainty_threshold, 3), 
                          " (", uncertainty_percentile, "th percentile)")
  
  # Create the main categorical plot
  cat_plot <- ggplot(viz_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = category)) +
    scale_fill_manual(values = category_colours, name = "Prediction Confidence") +
    coord_equal() +
    theme_minimal() +
    labs(
      title = plot_title,
      subtitle = subtitle_text
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Add border if provided
  if (!is.null(border_sf)) {
    cat_plot <- cat_plot + geom_sf(data = border_sf, fill = NA, colour = "black", size = 0.5)
  }
  
  # Create a continuous prediction with improved visibility
  cont_plot <- ggplot(viz_df, aes(x = x, y = y)) +
    # First add a base layer with just the prediction
    geom_raster(aes(fill = prediction)) +
    # Then add the prediction with uncertainty as transparency on top
    geom_raster(aes(fill = prediction, alpha = 1 - uncertainty)) +
    # Use plasma palette instead of viridis for more contrast
    scale_fill_viridis_c(name = "Habitat\nSuitability", option = "plasma") +
    # Increase the minimum alpha for better visibility
    scale_alpha_continuous(name = "Certainty", range = c(0.3, 0.9)) +
    coord_equal() +
    theme_minimal() +
    labs(title = "Continuous Prediction with Uncertainty") +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # Add border to continuous plot if provided
  if (!is.null(border_sf)) {
    cont_plot <- cont_plot + geom_sf(data = border_sf, fill = NA, colour = "black", size = 0.5)
  }
  
  # Calculate area statistics for the footer
  total_cells <- sum(!is.na(viz_df$presence))
  presence_cells <- sum(viz_df$presence == 1, na.rm = TRUE)
  high_conf_presence_cells <- sum(viz_df$category == high_conf_presence, na.rm = TRUE)
  
  # Create area stats text
  area_stats <- paste0(
    "Total area: ", format(total_cells, big.mark=","), " cells | ",
    "Suitable habitat: ", format(presence_cells, big.mark=","), " cells (", 
    round(presence_cells/total_cells*100, 1), "%) | ",
    "High confidence suitable: ", format(high_conf_presence_cells, big.mark=","), " cells (",
    round(high_conf_presence_cells/presence_cells*100, 1), "% of suitable)"
  )
  
  # Add caption to categorical plot
  cat_plot <- cat_plot + 
    labs(caption = area_stats) +
    theme(plot.caption = element_text(size = 9, hjust = 0))
  
  # Return as list
  plot_list <- list(
    categorical = cat_plot,
    continuous = cont_plot
  )
  
  return(plot_list)
} # End function uncertainty_viz

# Function to extract species and scenario info from filename
parse_filename <- function(filename) {
  # Extract parts from filename
  parts <- unlist(stringr::str_split(basename(filename), "_"))
  
  if (length(parts) >= 4) {
    species <- parts[1]
    
    # Check if it's a current or future projection
    if (parts[2] == "current") {
      time_period <- "current"
      scenario <- "current"
    } else {
      time_period <- parts[2]  # Year: 2040, 2060, 2080
      scenario <- parts[3]     # SSP: ssp126, ssp370, ssp585
    }
    
    return(list(
      species = species,
      time_period = time_period,
      scenario = scenario
    ))
  } else {
    return(NULL)
  }
}


#==============================================================================#
#                                   3. Run set up ----
#==============================================================================#

### Load all threshold values ----
for (species in species_list) {
  threshold_values[[species]] <- list()
  
  # Find threshold files for this species
  species_threshold_files <- grep(species, threshold_files, value = TRUE)
  
  for (threshold_type in threshold_types) {
    # Find threshold file for this type
    threshold_file <- grep(threshold_type, species_threshold_files, value = TRUE)
    
    if (length(threshold_file) > 0) {
      threshold_values[[species]][[threshold_type]] <- readRDS(threshold_file)
    } else {
      warning(paste("No threshold file found for", species, "with type", threshold_type))
    }
  }
}

### Create tracking data frame to record processing ----
results_tracker <- data.frame(
  species = character(),
  time_period = character(),
  scenario = character(),
  threshold_type = character(),
  suitable_area = numeric(),
  high_conf_suitable = numeric(),
  uncertainty_threshold = numeric(),
  stringsAsFactors = FALSE
)

# Main loop for processing all combinations ----
for (species in species_list) {
  message(paste("Processing species:", species))
  
  # Get the current projection for this species
  if (species == "calo") {
    current_proj <- calo_proj
  } else if (species == "vertlept") {
    current_proj <- vertlept_proj
  } else {
    next
  }
  
  # Get all future projections for this species
  species_predictions <- grep(species, prediction_list, value = TRUE)
  
  ### Create a list of all projections (current + future) ----
  all_projections <- c(
    # Add current projection
    list(list(
      file = paste0(species, "_current_current_prediction.tif"),
      raster = current_proj
    )),
    
    # Add future projections
    lapply(species_predictions, function(pred_file) {
      file_info <- parse_filename(pred_file)
      if (!is.null(file_info)) {
        return(list(
          file = basename(pred_file),
          raster = rast(pred_file),
          time_period = file_info$time_period,
          scenario = file_info$scenario
        ))
      }
      return(NULL)
    })
  )
  
  # Filter out any NULL entries
  all_projections <- all_projections[!sapply(all_projections, is.null)]
  
  # Process each projection with each threshold type ----
  for (proj in all_projections) {
    # Extract file info
    file_info <- parse_filename(proj$file)
    
    if (is.null(file_info)) {
      warning(paste("Could not parse filename:", proj$file))
      next
    }
    
    time_period <- file_info$time_period
    scenario <- file_info$scenario
    
    message(paste("  Time period:", time_period, "| Scenario:", scenario))
    
    # Process each threshold type
    for (threshold_type in threshold_types) {
      if (!is.null(threshold_values[[species]][[threshold_type]])) {
        threshold_value <- threshold_values[[species]][[threshold_type]]
        
        message(paste("    Using threshold type:", threshold_type, "=", round(threshold_value, 3)))
        
        # Extract prediction and uncertainty layers
        prediction_layer <- proj$raster[[1]]  # First layer is the prediction
        uncertainty_layer <- proj$raster[[2]]  # Second layer is uncertainty
        
        # Skip if either layer is missing
        if (is.null(prediction_layer) || is.null(uncertainty_layer)) {
          warning(paste("Missing layers for", proj$file))
          next
        }
        
        # Generate visualisations ----
        viz_plots <- uncertainty_viz(
          prediction_raster = prediction_layer,
          uncertainty_raster = uncertainty_layer,
          threshold = threshold_value,
          threshold_type = threshold_type,
          species_name = species,
          time_period = time_period,
          scenario = scenario,
          uncertainty_percentile = 75  # Using 75th percentile as recommended
        )
        
        # Create output file names
        output_base <- paste0(species, "_", time_period, "_", scenario, "_", threshold_type)
        cat_output <- file.path(output_dir, paste0(output_base, "_categorical.png"))
        cont_output <- file.path(output_dir, paste0(output_base, "_continuous.png"))
        
        # Save plots ----
        ggsave(cat_output, viz_plots$categorical, width = 8, height = 7, dpi = 300)
        ggsave(cont_output, viz_plots$continuous, width = 8, height = 7, dpi = 300)
        
        # Interpret for writing up ----
        # Extract area statistics from plot caption
        caption <- viz_plots$categorical$labels$caption
        
        # Parse area stats using regex
        suitable_area <- as.numeric(gsub(",", "", regmatches(
          caption, 
          regexpr("Suitable habitat: ([0-9,]+) cells", caption)
        )))
        
        high_conf_suitable <- as.numeric(gsub(",", "", regmatches(
          caption, 
          regexpr("High confidence suitable: ([0-9,]+) cells", caption)
        )))
        
        # Get uncertainty threshold from console output
        # We'll capture the actual values used in the visualisation
        uncertainty_values <- as.vector(as.matrix(uncertainty_layer))
        base_uncertainty_threshold <- stats::quantile(uncertainty_values, 0.75, na.rm = TRUE)
        
        # Determine final uncertainty threshold (following same logic as in uncertainty_viz function)
        final_uncertainty_percentile <- 75  # Start with recommended percentile
        final_uncertainty_threshold <- base_uncertainty_threshold
        
        # Check if the threshold is 0
        if (base_uncertainty_threshold == 0) {
          # Find first non-zero percentile
          for (test_percentile in seq(80, 99, by = 5)) {
            test_threshold <- stats::quantile(uncertainty_values, test_percentile/100, na.rm = TRUE)
            if (test_threshold > 0) {
              final_uncertainty_percentile <- test_percentile
              final_uncertainty_threshold <- test_threshold
              break
            }
          }
          
          # If all percentiles give 0, use mean or default
          if (final_uncertainty_threshold == 0) {
            mean_uncertainty <- mean(uncertainty_values, na.rm = TRUE)
            if (mean_uncertainty > 0) {
              final_uncertainty_threshold <- mean_uncertainty
              final_uncertainty_percentile <- NA  # Indicate we used mean instead
            } else {
              final_uncertainty_threshold <- 0.05  # Default value
              final_uncertainty_percentile <- NA  # Indicate we used default
            }
          }
        }
        # Record results ----
        # Record results with actual percentile used
        results_tracker <- rbind(results_tracker, data.frame(
          species = species,
          time_period = time_period,
          scenario = scenario,
          threshold_type = threshold_type,
          suitable_area = suitable_area,
          high_conf_suitable = high_conf_suitable,
          uncertainty_threshold = final_uncertainty_threshold,
          uncertainty_percentile = final_uncertainty_percentile,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Save all results tracker ----
write.csv(results_tracker, file.path(output_dir, "uncertainty_analysis_results.csv"), row.names = FALSE)

# Create summary visualisations ----
if (nrow(results_tracker) > 0) {
  # Load necessary libraries for summary plots
  library(ggplot2)
  library(dplyr)
  
  # Set factor levels for proper ordering
  results_tracker$time_period <- factor(
    results_tracker$time_period,
    levels = c("current", "2040", "2060", "2080")
  )
  
  results_tracker$scenario <- factor(
    results_tracker$scenario,
    levels = c("current", "ssp126", "ssp370", "ssp585")
  )
  
  # Calculate percentage of high confidence areas
  results_tracker$high_conf_percent <- results_tracker$high_conf_suitable / 
    results_tracker$suitable_area * 100
  
  # Confidence level change plot
  conf_change_plot <- ggplot(results_tracker, 
                             aes(x = time_period, y = high_conf_percent, 
                                 colour = scenario, group = interaction(scenario, threshold_type))) +
    geom_line() +
    geom_point(aes(shape = threshold_type), size = 3) +
    facet_wrap(~ species, scales = "free_y") +
    labs(title = "Change in High Confidence Predictions Over Time",
         x = "Time Period", 
         y = "High Confidence Area (%)",
         colour = "Climate Scenario",
         shape = "Threshold Type") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save summary plots
  ggsave(file.path(output_dir, "area_change_over_time.png"), 
         area_change_plot, width = 10, height = 7, dpi = 300)
  
  ggsave(file.path(output_dir, "confidence_change_over_time.png"), 
         conf_change_plot, width = 10, height = 7, dpi = 300)
}

message("Processing complete. Results saved to ", output_dir)

# Tidy up
rm(list = ls())
gc()

#==============================================================================#
#                             ---- End of workflow ----
#==============================================================================#