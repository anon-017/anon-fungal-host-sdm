# ---
# title: "07__visualisations.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-03-22"
# update: "2025-12-31"
# ---

#  Visualisations script that produces maps and charts from the outputs
# of 06_range_shift_metrics_consolidated.R script

#==============================================================================#
#                             0. Workspace setup ----
#==============================================================================#


## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

## Set up base folders ----
datafolder <- file.path("//xxx")


# Inputs: Range shift metric summary csv
data06_out <- file.path(datafolder, "06__rangeshift_metrics", "output")
range_shifts <- read.csv(file.path(data06_out, "metrics", "extra_settings", "summary", 
                                   "all_species_all_thresholds_metrics.csv"))

# Set up directories
data05_out_cur <- file.path(datafolder, "05__model_prediction", "output", "current_projections", "binary", "extra_settings")
data05_out_pred <- file.path(datafolder, "05__model_prediction", "output", "future_predictions", "binary", "extra_settings")
data07_out <- file.path(datafolder, "07__visualisations", "output")

# Create output figure directories
figure_path <- file.path(data07_out, "figures", "extra_settings")
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

# Create only the directories we need for our simplified visualisation
output_dirs <- c(
  "combined_metrics", "range_change", "bivariate"
)
for(dir in output_dirs) {
  dir.create(file.path(figure_path, dir), recursive = TRUE, showWarnings = FALSE)
}

# Set up logging
log_file <- file.path(data07_out, "processing_log.txt")
cat("Processing started at", Sys.time(), "\n", file=log_file)



# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggtext)
library(terra)
library(scales)
#library(purrr) # For functional programming

#==============================================================================#
#                       1. Configuration variables ----
#==============================================================================#

# Global variables
thresholds <- c("opt2", "opt4", "opt12") 
scenarios <- c("ssp126", "ssp370", "ssp585")
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

# Define species colours
calo_colour <- "#006400"      # Dark green for Calophyllum
vertlept_colour <- "#6A0DAD"  # Deep purple for Verticillium/Leptographium

# Define species labels with proper formatting
species_labels <- c(
  "calo" = expression(italic("Calophyllum paniculatum")),
  "vertlept" = expression(italic("Verticillium")~"/"~italic("Leptographium"))
)

# Define threshold pairs for bivariate analysis
threshold_pairs <- data.frame(
  calo = c("opt2", "opt12"),
  vertlept = c("opt4", "opt12")
)

#==============================================================================#
#                             2. Helper functions ----
#==============================================================================#

# Function to format species names in italics for ggplot
italic_species <- function(x) {
  # Only italicize genus+species names
  if(grepl(" ", x)) {
    # For species names with spaces, wrap in italics
    return(paste0("italic('", x, "')"))
  } else {
    # For codes or names without spaces, return as is
    return(x)
  }
}

# Function to save plots with consistent naming
save_plot <- function(plot, name, figure_path, width = 8, height = 6) {
  ggsave(
    filename = file.path(figure_path, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  # Log information
  cat(paste0("Saved plot: ", name, ".png\n"), file = log_file, append = TRUE)
}

# Function to get threshold description
get_threshold_description <- function(threshold_type) {
  switch(threshold_type,
         "opt2" = "Maximise Sum of Sensitivity and Specificity",
         "opt12" = "10th Percentile Training Presence",
         "opt4" = "Minimum Distance to ROC Curve",
         paste0("Unknown threshold (", threshold_type, ")"))
}

# Function to list and categorize raster files
list_raster_files <- function(dir_path) {
  # Log the function call
  cat("Listing raster files in:", dir_path, "\n", file = log_file, append = TRUE)
  
  # List all files in the directory with .tif extension
  all_files <- list.files(dir_path, pattern = "\\.tif$|\\.asc$|\\.bin$", 
                          full.names = TRUE, recursive = TRUE)
  
  # Check if any files were found
  if (length(all_files) == 0) {
    cat("No raster files found in:", dir_path, "\n", file = log_file, append = TRUE)
    return(data.frame())
  }
  
  # Log the number of files found
  cat("Found", length(all_files), "raster files\n", file = log_file, append = TRUE)
  
  # Extract metadata from filenames
  file_info <- data.frame(filepath = all_files, stringsAsFactors = FALSE)
  
  # Extract just the filename without the path
  file_info$filename <- basename(file_info$filepath)
  
  # Extract species information
  file_info$species <- ifelse(grepl("calo", file_info$filename, ignore.case = TRUE), 
                              "calo", "vertlept")
  
  # Extract threshold information
  file_info$threshold <- NA
  file_info$threshold[grepl("opt2", file_info$filename, ignore.case = TRUE)] <- "opt2"
  file_info$threshold[grepl("opt4", file_info$filename, ignore.case = TRUE)] <- "opt4"
  file_info$threshold[grepl("opt12", file_info$filename, ignore.case = TRUE)] <- "opt12"
  
  # Extract scenario information
  file_info$scenario <- NA
  file_info$scenario[grepl("ssp126", file_info$filename, ignore.case = TRUE)] <- "ssp126"
  file_info$scenario[grepl("ssp370", file_info$filename, ignore.case = TRUE)] <- "ssp370"
  file_info$scenario[grepl("ssp585", file_info$filename, ignore.case = TRUE)] <- "ssp585"
  file_info$scenario[!grepl("ssp", file_info$filename, ignore.case = TRUE)] <- "current"
  
  # Extract time period information
  file_info$time_period <- NA
  file_info$time_period[grepl("2011-2040", file_info$filename, ignore.case = TRUE)] <- "2011-2040"
  file_info$time_period[grepl("2041-2070", file_info$filename, ignore.case = TRUE)] <- "2041-2070"
  file_info$time_period[grepl("2071-2100", file_info$filename, ignore.case = TRUE)] <- "2071-2100"
  
  # Log a summary of what we found
  cat("Categorized rasters by species:\n", table(file_info$species), "\n", 
      file = log_file, append = TRUE)
  cat("Categorized rasters by threshold:\n", table(file_info$threshold), "\n", 
      file = log_file, append = TRUE)
  cat("Categorized rasters by scenario:\n", table(file_info$scenario), "\n", 
      file = log_file, append = TRUE)
  
  return(file_info)
}

# Function to find and load the best matching raster
find_best_raster <- function(file_list, species, threshold, scenario = NULL, time_period = NULL) {
  # Log what we're looking for
  cat("Looking for raster: species=", species, 
      ", threshold=", threshold, 
      ", scenario=", ifelse(is.null(scenario), "current", scenario), 
      ", time_period=", ifelse(is.null(time_period), "NA", time_period), 
      "\n", file = log_file, append = TRUE)
  
  # Filter file list to find matching files
  matches <- file_list %>%
    filter(
      species == !!species,
      threshold == !!threshold
    )
  
  # Add scenario filter if provided
  if (!is.null(scenario)) {
    matches <- matches %>% filter(scenario == !!scenario)
  }
  
  # Add time period filter if provided
  if (!is.null(time_period)) {
    matches <- matches %>% filter(time_period == !!time_period)
  }
  
  # Check if we found any matches
  if (nrow(matches) == 0) {
    cat("No matching raster found\n", file = log_file, append = TRUE)
    return(NULL)
  }
  
  # If multiple matches, use the first one but log a warning
  if (nrow(matches) > 1) {
    cat("Warning: Found", nrow(matches), "matching rasters, using the first one\n", 
        file = log_file, append = TRUE)
    cat("Matches:", paste(matches$filename, collapse=", "), "\n", 
        file = log_file, append = TRUE)
  }
  
  # Get the file path of the best match
  best_match <- matches$filepath[1]
  cat("Loading raster from:", best_match, "\n", file = log_file, append = TRUE)
  
  # Try to load the raster
  tryCatch({
    # Load the raster at its full extent
    raster <- terra::rast(best_match)
    
    # Log raster properties
    cat("Loaded raster with dimensions:", dim(raster)[1], "x", dim(raster)[2], "\n", 
        file = log_file, append = TRUE)
    
    # Return the raster and metadata
    return(list(
      raster = raster,
      filepath = best_match,
      filename = basename(best_match)
    ))
  }, 
  error = function(e) {
    cat("Error loading raster:", e$message, "\n", file = log_file, append = TRUE)
    return(NULL)
  })
}

# Standardized data frame conversion function
convert_raster_to_df <- function(raster, col_name = "value") {
  if (is.null(raster)) {
    return(data.frame())
  }
  
  df <- as.data.frame(raster, xy = TRUE)
  names(df)[3] <- col_name
  return(df)
}

# Function to ensure consistent geometry between rasters
ensure_compatible_geometry <- function(reference_raster, target_raster) {
  if (is.null(reference_raster) || is.null(target_raster)) {
    return(target_raster)
  }
  
  if (!terra::compareGeom(reference_raster, target_raster)) {
    return(terra::resample(target_raster, reference_raster))
  }
  return(target_raster)
}

# Function to calculate overlap change
calculate_overlap_change <- function(calo_current, vertlept_current, 
                                     calo_future, vertlept_future) {
  # Make sure all rasters have the same extent and resolution
  calo_future <- ensure_compatible_geometry(calo_current, calo_future)
  vertlept_current <- ensure_compatible_geometry(calo_current, vertlept_current)
  vertlept_future <- ensure_compatible_geometry(calo_current, vertlept_future)
  
  # Calculate overlap for current and future
  current_overlap <- calo_current * vertlept_current
  future_overlap <- calo_future * vertlept_future
  
  # Calculate change in overlap
  overlap_change <- future_overlap - current_overlap
  
  return(overlap_change)
}

# Standardized plotting theme
get_map_theme <- function() {
  theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.1, "cm"),
      legend.text = element_text(size = 9)
    )
}

# Get metrics theme for bar charts
get_metrics_theme <- function() {
  theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom",
      legend.text = element_text(face = "italic"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# Create a faceted grid of overlap changes
create_overlap_change_facet_grid <- function(calo_current, vertlept_current,
                                             future_file_list, scenarios, time_periods,
                                             threshold_calo, threshold_vertlept,
                                             output_path) {
  # Create base mask for plotting background
  mask_df <- convert_raster_to_df(calo_current)
  mask_df$is_valid <- !is.na(mask_df$value)
  
  # Background colour
  bg_colour <- "#D0D0D0"
  
  # Calculate current overlap
  current_overlap <- calo_current * vertlept_current
  
  # Create a data frame to hold all overlap changes
  all_overlap_changes <- data.frame()
  
  # Process each scenario and time period
  for (scenario in scenarios) {
    for (time_period in time_periods) {
      # Find and load the rasters
      calo_result <- find_best_raster(future_file_list, "calo", threshold_calo, scenario, time_period)
      vertlept_result <- find_best_raster(future_file_list, "vertlept", threshold_vertlept, scenario, time_period)
      
      # Skip this scenario/time period if files not found
      if (is.null(calo_result) || is.null(vertlept_result)) {
        cat(paste0("Skipping overlap change for scenario ", scenario, 
                   ", time period ", time_period, " due to missing files.\n"), 
            file = log_file, append = TRUE)
        next
      }
      
      calo_future <- calo_result$raster
      vertlept_future <- vertlept_result$raster
      
      # Ensure compatible geometry
      calo_future <- ensure_compatible_geometry(calo_current, calo_future)
      vertlept_future <- ensure_compatible_geometry(calo_current, vertlept_future)
      
      # Calculate overlap and change
      future_overlap <- calo_future * vertlept_future
      overlap_change <- future_overlap - current_overlap
      
      # Convert to dataframe
      change_df <- convert_raster_to_df(overlap_change, "change")
      
      # Add scenario and time period info
      change_df$scenario <- scenario
      change_df$time_period <- time_period
      
      # Only include non-zero changes
      change_df <- subset(change_df, change != 0)
      
      # Append to the main dataframe
      all_overlap_changes <- rbind(all_overlap_changes, change_df)
    }
  }
  
  # Create the facet plot of overlap changes
  change_facet_plot <- ggplot() +
    # Background
    geom_raster(data = mask_df, aes(x = x, y = y), fill = bg_colour) +
    # Overlap changes
    geom_raster(data = all_overlap_changes, 
                aes(x = x, y = y, fill = factor(change))) +
    # Facets
    facet_grid(scenario ~ time_period) +
    # colours
    scale_fill_manual(
      values = c("-1" = "#D55E00", "1" = "#0072B2"),
      name = "Overlap Change",
      labels = c("Lost Overlap", "New Overlap")
    ) +
    labs(title = "Overlap Changes") +
    get_map_theme() +
    coord_equal()
  
  # Create caption text with thresholds 
  caption_text <- paste(
    "Thresholds:",
    paste0("Calophyllum: ", get_threshold_description(threshold_calo)),
    paste0("Verticillium: ", get_threshold_description(threshold_vertlept)),
    sep = "  "
  )
  
  # Add caption
  change_facet_plot_with_caption <- change_facet_plot +
    labs(caption = caption_text) +
    theme(
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10))
    )
  
  # Save the facet plot
  filename <- paste0("overlap_changes_facet_", 
                     threshold_calo, "_", threshold_vertlept, ".png")
  ggsave(
    file.path(output_path, filename),
    change_facet_plot_with_caption,
    width = 12,
    height = 12,
    dpi = 300
  )
  
  cat(paste0("Created overlap change facet plot: ", filename, "\n"), 
      file = log_file, append = TRUE)
  
  return(change_facet_plot_with_caption)
}

# Main function to process all bivariate visualisations
run_bivariate_analysis <- function(current_dir, future_dir, output_dir, 
                                   threshold_pairs, scenarios, time_periods) {
  # Record start time for logging
  start_time <- Sys.time()
  cat(paste("Starting bivariate analysis at", start_time, "\n"), 
      file = log_file, append = TRUE)
  
  # List all raster files in both directories
  current_file_list <- list_raster_files(current_dir)
  future_file_list <- list_raster_files(future_dir)
  
  # Validate file lists
  if (nrow(current_file_list) == 0) {
    stop("No current raster files found in the specified directory.")
  }
  
  if (nrow(future_file_list) == 0) {
    stop("No future raster files found in the specified directory.")
  }
  
  # Log file counts
  cat(paste0("Found ", nrow(current_file_list), " current raster files.\n"), 
      file = log_file, append = TRUE)
  cat(paste0("Found ", nrow(future_file_list), " future raster files.\n"), 
      file = log_file, append = TRUE)
  
  # Create bivariate output directory if it doesn't exist
  bivariate_dir <- file.path(output_dir, "bivariate")
  dir.create(bivariate_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each threshold pair
  for (i in 1:nrow(threshold_pairs)) {
    calo_threshold <- threshold_pairs[i, "calo"]
    vertlept_threshold <- threshold_pairs[i, "vertlept"]
    
    cat(paste0("\nProcessing threshold combination: Calo=", calo_threshold, 
               ", Vertlept=", vertlept_threshold, "\n"), 
        file = log_file, append = TRUE)
    
    # Find and load current distribution rasters
    calo_current_result <- find_best_raster(current_file_list, "calo", calo_threshold, "current")
    vertlept_current_result <- find_best_raster(current_file_list, "vertlept", vertlept_threshold, "current")
    
    if (is.null(calo_current_result) || is.null(vertlept_current_result)) {
      cat("ERROR: Could not find current rasters for thresholds: ", 
          calo_threshold, ", ", vertlept_threshold, "\n", 
          file = log_file, append = TRUE)
      next
    }
    
    # Load the full raster data
    calo_current <- calo_current_result$raster
    vertlept_current <- vertlept_current_result$raster
    
    # Make sure we're working with the full extent
    cat("Current raster dimensions - Calo:", dim(calo_current)[1], "x", dim(calo_current)[2], 
        ", Vertlept:", dim(vertlept_current)[1], "x", dim(vertlept_current)[2], "\n", 
        file = log_file, append = TRUE)
    
    # 1. Create individual plots for each scenario and time period
    for (scenario in scenarios) {
      for (time_period in time_periods) {
        # Find and load future rasters
        calo_result <- find_best_raster(future_file_list, "calo", calo_threshold, scenario, time_period)
        vertlept_result <- find_best_raster(future_file_list, "vertlept", vertlept_threshold, scenario, time_period)
        
        # Skip this scenario/time period if files not found
        if (is.null(calo_result) || is.null(vertlept_result)) {
          cat(paste0("Skipping visualisation for scenario ", scenario, 
                     ", time period ", time_period, " due to missing files.\n"), 
              file = log_file, append = TRUE)
          next
        }
        
        # Load the full future rasters
        calo_future <- calo_result$raster
        vertlept_future <- vertlept_result$raster
        
        # Log the dimensions of the future rasters
        cat("Future raster dimensions -", scenario, time_period, "- Calo:", 
            dim(calo_future)[1], "x", dim(calo_future)[2], 
            ", Vertlept:", dim(vertlept_future)[1], "x", dim(vertlept_future)[2], "\n", 
            file = log_file, append = TRUE)
        
        # Make sure geometries match (don't crop or reduce resolution)
        if (!terra::compareGeom(calo_current, calo_future)) {
          cat("Resampling future calo raster to match current geometry\n", file = log_file, append = TRUE)
          calo_future <- terra::resample(calo_future, calo_current, method = "near")
        }
        
        if (!terra::compareGeom(calo_current, vertlept_future)) {
          cat("Resampling future vertlept raster to match current geometry\n", file = log_file, append = TRUE)
          vertlept_future <- terra::resample(vertlept_future, calo_current, method = "near")
        }
        
        # Create the three-panel plot
        cat("Creating bivariate map for", scenario, time_period, "\n", file = log_file, append = TRUE)
        create_bivariate_overlap_map(
          calo_current, vertlept_current,
          calo_future, vertlept_future,
          scenario, time_period,
          calo_threshold, vertlept_threshold,
          bivariate_dir
        )
      }
    }
    
    # 2. Create facet grids for this threshold combination
    cat("Creating future distribution facet grid\n", file = log_file, append = TRUE)
    create_future_distribution_facet_grid(
      calo_current, vertlept_current,
      future_file_list, scenarios, time_periods,
      calo_threshold, vertlept_threshold,
      bivariate_dir
    )
    
    cat("Creating overlap change facet grid\n", file = log_file, append = TRUE)
    create_overlap_change_facet_grid(
      calo_current, vertlept_current,
      future_file_list, scenarios, time_periods,
      calo_threshold, vertlept_threshold,
      bivariate_dir
    )
  }
  
  # Record end time
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat(paste0("\nBivariate analysis completed in ", round(elapsed, 2), " minutes.\n"), 
      file = log_file, append = TRUE)
  
  return(bivariate_dir)
}


#==============================================================================#
#                            3. Data preparation ----
#==============================================================================#

# Clean and prepare range shifts data
prepare_range_shifts_data <- function(range_shifts) {
  # Clean and prepare data
  range_shifts <- range_shifts %>%
    # Remove rows with missing values in key columns
    filter(!is.na(percent_stability) & !is.na(percent_expansion) & 
             !is.na(percent_contraction) & !is.na(habitat_distance_km)) %>%
    # Remove problematic columns
    dplyr::select(-matches("centroid_shift_velocity_km_yr")) %>%
    dplyr::select(-matches("centroid_shift_velocity_km_dec")) %>%
    # Convert factors for better plotting
    mutate(
      species = factor(species),
      species_fullname = factor(species_fullname),
      scenario = factor(scenario),
      time_period = factor(time_period),
      threshold = factor(threshold)
    )
  
  return(range_shifts)
}

# Create threshold lookup table
create_threshold_lookup <- function(range_shifts) {
  data.frame(
    threshold = unique(range_shifts$threshold)
  ) %>%
    mutate(threshold_desc = sapply(as.character(threshold), get_threshold_description))
}

# Apply data preparation
range_shifts <- prepare_range_shifts_data(range_shifts)
threshold_lookup <- create_threshold_lookup(range_shifts)

#==============================================================================#
#                         4. Range change visualisation ----
#==============================================================================#

# Create range change comparison plot
create_range_change_plot <- function(range_shifts, threshold_lookup) {
  range_shifts %>%
    # Join with threshold lookup to get descriptions
    left_join(threshold_lookup, by = "threshold") %>%
    dplyr::select(species_fullname, scenario, time_period, threshold, threshold_desc, 
                  percent_stability, percent_expansion, percent_contraction) %>%
    pivot_longer(cols = c(percent_stability, percent_expansion, percent_contraction),
                 names_to = "component", values_to = "percentage") %>%
    dplyr::mutate(component = factor(component, 
                                     levels = c("percent_contraction", "percent_stability", "percent_expansion"),
                                     labels = c("Contraction", "Stability", "Expansion")),
                  # Create facet label with species and threshold
                  facet_label = paste0("<i>", species_fullname, "</i><br><br>", "Threshold: ", threshold_desc)) %>%
    ggplot(aes(x = scenario, y = percentage, fill = component)) +
    geom_bar(stat = "identity") +
    facet_grid(time_period ~ facet_label) +
    # colour scheme: red for contraction, grey for stability, blue for expansion
    scale_fill_manual(values = c("Contraction" = "#D55E00", 
                                 "Stability" = "lightyellow",
                                 "Expansion" = "#0072B2"),
                      name = "Range Component") +
    labs(title = "Range Change Components by Species, Threshold and Scenario",
         x = "Climate Scenario", 
         y = "Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_markdown(size = 8, lineheight = 1.2),
          panel.spacing.y = unit(1, "lines"),
          legend.position = "bottom")
}

#==============================================================================#
#                       5. Metrics visualisation ----
#==============================================================================#

# Flexible function to create metrics visualisations
# Corrected flexible function to create metrics visualisations
create_metrics_visualisation <- function(data, threshold_val = NULL, 
                                         scenario_val = NULL, time_period_val = NULL,
                                         facet_by = "threshold") {
  # Filter data based on provided filters
  filtered_data <- data
  
  if (!is.null(threshold_val)) {
    filtered_data <- filtered_data %>% filter(threshold == threshold_val)
  }
  
  if (!is.null(scenario_val)) {
    filtered_data <- filtered_data %>% filter(scenario == scenario_val)
  }
  
  if (!is.null(time_period_val)) {
    filtered_data <- filtered_data %>% filter(time_period == time_period_val)
  }
  
  # Join with threshold lookup for better labels
  filtered_data <- filtered_data %>%
    left_join(threshold_lookup, by = "threshold")
  
  # Determine the x-axis based on faceting choice
  x_var <- switch(facet_by,
                  "threshold" = "scenario",
                  "scenario" = "threshold",
                  "time_period" = "scenario",
                  "scenario") # Default
  
  # Select the columns we need for the metrics
  metrics_data <- filtered_data %>%
    dplyr::select(species, !!sym(x_var), time_period, threshold_desc, 
                  habitat_distance_km, habitat_exposure_km2, 
                  spatial_disruption, aoo_percent_change)
  
  # Determine facet variables and formula after we've already selected the base data
  if (facet_by == "threshold") {
    facet_formula <- "metric ~ time_period"
    title_suffix <- paste("Threshold:", unique(filtered_data$threshold_desc)[1])
  } else if (facet_by == "scenario") {
    facet_formula <- "metric ~ threshold_desc"
    title_suffix <- paste("Scenario:", scenario_val, "| Time Period:", time_period_val)
  } else if (facet_by == "time_period") {
    facet_formula <- "metric ~ threshold_desc"
    title_suffix <- paste("Time Period:", time_period_val)
  }
  
  # Transform the data
  metrics_data <- metrics_data %>%
    # Scale metrics to make them comparable on visualisation
    mutate(
      habitat_distance_km = habitat_distance_km,
      habitat_exposure_km2 = habitat_exposure_km2,
      spatial_disruption = spatial_disruption * 100, # Convert to percentage
      aoo_percent_change = aoo_percent_change
    ) %>%
    # Convert to long format
    pivot_longer(
      cols = c(habitat_distance_km, habitat_exposure_km2, 
               spatial_disruption, aoo_percent_change),
      names_to = "metric",
      values_to = "value"
    ) %>%
    # Create better labels for metrics
    mutate(
      metric = factor(
        metric,
        levels = c("habitat_distance_km", "habitat_exposure_km2", 
                   "spatial_disruption", "aoo_percent_change"),
        labels = c("Habitat Distance (km)", "Habitat Exposure (km²)", 
                   "Spatial Disruption (%)", "AOO Change (%)")
      )
    )
  
  # Build the facet formula from the facet variables
  facet_formula <- as.formula(facet_formula)
  
  # Create the plot
  metrics_plot <- ggplot(metrics_data, aes(x = !!sym(x_var), y = value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(facet_formula, scales = "free_y") +
    scale_fill_manual(values = c("calo" = calo_colour, "vertlept" = vertlept_colour), 
                      name = "Species",
                      labels = species_labels) +
    labs(title = paste("Key Range Shift Metrics -", title_suffix),
         x = ifelse(x_var == "threshold", "Threshold Type", 
                    ifelse(x_var == "scenario", "Climate Scenario", "Time Period")), 
         y = "Value") +
    get_metrics_theme()
  
  return(metrics_plot)
}

# Create metrics plots for different views
create_all_metrics_plots <- function(range_shifts, threshold_lookup) {
  # For each threshold
  walk(unique(range_shifts$threshold), function(thresh) {
    metrics_plot <- create_metrics_visualisation(
      range_shifts, 
      threshold_val = thresh,
      facet_by = "threshold"
    )
    
    save_plot(metrics_plot, 
              paste0("metrics_by_threshold_", thresh), 
              file.path(figure_path, "combined_metrics"), 
              width = 10, height = 12)
  })
  
  # For each scenario and time period combination
  for (scenario_val in unique(range_shifts$scenario)) {
    for (time_period_val in unique(range_shifts$time_period)) {
      metrics_plot <- create_metrics_visualisation(
        range_shifts,
        scenario_val = scenario_val,
        time_period_val = time_period_val,
        facet_by = "scenario"
      )
      
      save_plot(metrics_plot,
                paste0("metrics_", scenario_val, "_", time_period_val),
                file.path(figure_path, "combined_metrics"),
                width = 10, height = 8)
    }
  }
}



#==============================================================================#
#                          6. Summary visualisation ----
#==============================================================================#


# Create a summary visualisation that shows key metrics across all scenarios
create_summary_plot <- function(data) {
  # Prepare data for summary plot
  summary_data <- data %>%
    group_by(species, scenario, time_period) %>%
    summarize(
      mean_stability = mean(percent_stability, na.rm = TRUE),
      mean_disruption = mean(spatial_disruption * 100, na.rm = TRUE),
      mean_aoo_change = mean(aoo_percent_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = starts_with("mean_"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("mean_stability", "mean_disruption", "mean_aoo_change"),
        labels = c("Range Stability (%)", 
                   "Spatial Disruption (%)", "AOO Change (%)")
      )
    )
  
  # Create summary plot
  summary_plot <- ggplot(summary_data, aes(x = scenario, y = value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(metric ~ time_period, scales = "free_y") +
    scale_fill_manual(
      values = c("calo" = calo_colour, "vertlept" = vertlept_colour),
      name = "Species",
      labels = species_labels
    ) +
    labs(
      title = "Summary of Range Shift Metrics Across All Thresholds",
      subtitle = "Averaged across all threshold types",
      x = "Climate Scenario",
      y = "Value"
    ) +
    get_metrics_theme()
  
  # Save the summary plot
  save_plot(
    summary_plot,
    "summary_all_metrics",
    figure_path,
    width = 12,
    height = 10
  )
  
  return(summary_plot)
}


#==============================================================================#
#                 7. Bivariate distribution Maps ----
#==============================================================================#

# Create a plot for current distribution
create_current_distribution_plot <- function(calo_current, vertlept_current, title = "Current Distribution", 
                                             show_legend = TRUE) {
  # Calculate overlap
  current_overlap <- calo_current * vertlept_current
  
  # Convert to data frames with consistent naming
  # Use entire extent instead of creating small squares
  mask_df <- convert_raster_to_df(calo_current)
  mask_df$is_valid <- !is.na(mask_df$value)
  
  calo_df <- convert_raster_to_df(calo_current)
  calo_df <- subset(calo_df, value > 0)
  
  vertlept_df <- convert_raster_to_df(vertlept_current)
  vertlept_df <- subset(vertlept_df, value > 0)
  
  overlap_df <- convert_raster_to_df(current_overlap)
  overlap_df <- subset(overlap_df, value > 0)
  
  # Background colour
  bg_colour <- "#D0D0D0"
  
  # Create the plot
  plot <- ggplot() +
    # Background/mask with darker grey
    geom_raster(data = mask_df, aes(x = x, y = y), fill = bg_colour) +
    # Calo distribution (green)
    geom_raster(data = calo_df, aes(x = x, y = y), fill = "#006400", alpha = 0.6) +
    # Vertlept distribution (purple)  
    geom_raster(data = vertlept_df, aes(x = x, y = y), fill = "#6A0DAD", alpha = 0.6) + 
    # Current overlap (black)
    geom_raster(data = overlap_df, aes(x = x, y = y), fill = "#000000")
  
  # Add legend only if requested
  if (show_legend) {
    plot <- plot +
      # Create invisible points for legend with proper order
      geom_point(data = data.frame(
        species = factor(c("Calophyllum", "Verticillium", "Overlap"),
                         levels = c("Calophyllum", "Verticillium", "Overlap")),
        x = -Inf, y = -Inf),
        aes(x = x, y = y, colour = species), alpha = 0) +
      scale_colour_manual(
        name = "Distribution",
        values = c("Calophyllum" = "#006400", 
                   "Verticillium" = "#6A0DAD", 
                   "Overlap" = "#000000"),
        labels = c(
          expression(italic("Calophyllum paniculatum")),
          expression(italic("Verticillium/Leptographium")),
          "Overlap"
        ),
        guide = guide_legend(override.aes = list(
          alpha = c(0.6, 0.6, 1),
          size = 5,
          shape = 15
        ))
      )
  }
  
  # Add title and theme
  plot <- plot +
    labs(title = title) +
    get_map_theme() +
    coord_equal()
  
  return(plot)
}

# Create a plot for future distribution
create_future_distribution_plot <- function(calo_future, vertlept_future, scenario, time_period,
                                            show_legend = TRUE) {
  # Calculate overlap
  future_overlap <- calo_future * vertlept_future
  
  # Convert to data frames with consistent naming
  mask_df <- convert_raster_to_df(calo_future)
  mask_df$is_valid <- !is.na(mask_df$value)
  
  calo_df <- convert_raster_to_df(calo_future)
  calo_df <- subset(calo_df, value > 0)
  
  vertlept_df <- convert_raster_to_df(vertlept_future)
  vertlept_df <- subset(vertlept_df, value > 0)
  
  overlap_df <- convert_raster_to_df(future_overlap)
  overlap_df <- subset(overlap_df, value > 0)
  
  # Background colour
  bg_colour <- "#D0D0D0"
  
  # Create the plot
  plot <- ggplot() +
    geom_raster(data = mask_df, aes(x = x, y = y), fill = bg_colour) +
    geom_raster(data = calo_df, aes(x = x, y = y), fill = "#006400", alpha = 0.6) +
    geom_raster(data = vertlept_df, aes(x = x, y = y), fill = "#6A0DAD", alpha = 0.6) + 
    geom_raster(data = overlap_df, aes(x = x, y = y), fill = "#000000")
  
  # Add legend only if requested
  if (show_legend) {
    plot <- plot +
      geom_point(data = data.frame(
        species = factor(c("Calophyllum", "Verticillium", "Overlap"),
                         levels = c("Calophyllum", "Verticillium", "Overlap")),
        x = -Inf, y = -Inf),
        aes(x = x, y = y, colour = species), alpha = 0) +
      scale_colour_manual(
        name = "Distribution",
        values = c("Calophyllum" = "#006400", 
                   "Verticillium" = "#6A0DAD", 
                   "Overlap" = "#000000"),
        labels = c(
          expression(italic("Calophyllum paniculatum")),
          expression(italic("Verticillium/Leptographium")),
          "Overlap"
        ),
        guide = guide_legend(override.aes = list(
          alpha = c(0.6, 0.6, 1),
          size = 5,
          shape = 15
        ))
      )
  }
  
  # Add title and theme
  plot <- plot +
    labs(title = paste("Future Distribution -", scenario, time_period)) +
    get_map_theme() +
    coord_equal()
  
  return(plot)
}

# Create a plot for overlap change
create_overlap_change_plot <- function(overlap_change, scenario, time_period) {
  # Convert to data frame
  mask_df <- convert_raster_to_df(overlap_change)
  mask_df$is_valid <- !is.na(mask_df$value)
  
  change_df <- convert_raster_to_df(overlap_change, "change")
  change_df <- subset(change_df, change != 0)
  
  # Background colour
  bg_colour <- "#D0D0D0"
  
  # Create the plot
  ggplot() +
    geom_raster(data = mask_df, aes(x = x, y = y), fill = bg_colour) +
    geom_raster(data = change_df, 
                aes(x = x, y = y, fill = factor(change))) +
    scale_fill_manual(
      values = c("-1" = "#D55E00", "1" = "#0072B2"),
      name = "Overlap Change",
      labels = c("Lost Overlap", "New Overlap")
    ) +
    labs(title = paste("Overlap Change -", scenario, time_period)) +
    get_map_theme() +
    coord_equal()
}

# Create three-panel visualisation for a single scenario and time period
create_bivariate_overlap_map <- function(calo_current, vertlept_current, 
                                         calo_future, vertlept_future,
                                         scenario, time_period,
                                         threshold_calo, threshold_vertlept,
                                         output_path) {
  # Ensure all rasters have the same extent and resolution
  calo_future <- ensure_compatible_geometry(calo_current, calo_future)
  vertlept_current <- ensure_compatible_geometry(calo_current, vertlept_current)
  vertlept_future <- ensure_compatible_geometry(calo_current, vertlept_future)
  
  # Calculate overlap change
  overlap_change <- calculate_overlap_change(calo_current, vertlept_current, 
                                             calo_future, vertlept_future)
  
  # Create the three individual plots - only show legend on the first plot
  current_map <- create_current_distribution_plot(calo_current, vertlept_current, show_legend = TRUE)
  future_map <- create_future_distribution_plot(calo_future, vertlept_future, scenario, time_period, show_legend = FALSE)
  change_map <- create_overlap_change_plot(overlap_change, scenario, time_period)
  
  # Create caption text with thresholds 
  caption_text <- paste(
    "Thresholds:",
    paste0("Calophyllum: ", get_threshold_description(threshold_calo)),
    paste0("Verticillium: ", get_threshold_description(threshold_vertlept)),
    sep = "  "
  )
  
  # Combine all three plots horizontally
  combined_map <- current_map + future_map + change_map +
    plot_layout(ncol = 3, widths = c(1, 1, 1)) +
    plot_annotation(
      title = paste("Distribution Comparison -", scenario, time_period),
      caption = caption_text,
      theme = theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10))
      )
    )
  
  # Save the combined map
  filename <- paste0(scenario, "_", time_period, "_", 
                     threshold_calo, "_", threshold_vertlept, "_overlap.png")
  ggsave(
    file.path(output_path, filename),
    combined_map,
    width = 15,  # Wider to accommodate horizontal layout
    height = 6,   # Shorter since plots are side by side
    dpi = 300
  )
  
  cat(paste0("Created bivariate map: ", filename, "\n"), file = log_file, append = TRUE)
  
  return(combined_map)
}

# Create a faceted grid of future distributions
create_future_distribution_facet_grid <- function(calo_current, vertlept_current,
                                                  future_file_list, scenarios, time_periods,
                                                  threshold_calo, threshold_vertlept,
                                                  output_path) {
  # Create base mask for plotting background - use full extent
  mask_df <- convert_raster_to_df(calo_current)
  mask_df$is_valid <- !is.na(mask_df$value)
  
  # Background colour
  bg_colour <- "#D0D0D0"
  
  # Create a data frame to hold all future distributions
  all_future_distributions <- data.frame()
  
  # Process each scenario and time period
  for (scenario in scenarios) {
    for (time_period in time_periods) {
      # Find and load the rasters
      calo_result <- find_best_raster(future_file_list, "calo", threshold_calo, scenario, time_period)
      vertlept_result <- find_best_raster(future_file_list, "vertlept", threshold_vertlept, scenario, time_period)
      
      # Skip this scenario/time period if files not found
      if (is.null(calo_result) || is.null(vertlept_result)) {
        cat(paste0("Skipping future distribution for scenario ", scenario, 
                   ", time period ", time_period, " due to missing files.\n"), 
            file = log_file, append = TRUE)
        next
      }
      
      calo_future <- calo_result$raster
      vertlept_future <- vertlept_result$raster
      
      # Ensure compatible geometry
      calo_future <- ensure_compatible_geometry(calo_current, calo_future)
      vertlept_future <- ensure_compatible_geometry(calo_current, vertlept_future)
      
      # Calculate overlap
      future_overlap <- calo_future * vertlept_future
      
      # Convert to dataframes using full extent
      calo_df <- convert_raster_to_df(calo_future)
      vertlept_df <- convert_raster_to_df(vertlept_future)
      overlap_df <- convert_raster_to_df(future_overlap)
      
      # Add scenario and time period info
      calo_df$scenario <- scenario
      calo_df$time_period <- time_period
      calo_df$type <- "Calophyllum"
      
      vertlept_df$scenario <- scenario
      vertlept_df$time_period <- time_period
      vertlept_df$type <- "Verticillium"
      
      overlap_df$scenario <- scenario
      overlap_df$time_period <- time_period
      overlap_df$type <- "Overlap"
      
      # Combine all data
      scenario_data <- rbind(
        calo_df,
        vertlept_df,
        overlap_df
      )
      
      # Append to the main dataframe
      all_future_distributions <- rbind(all_future_distributions, scenario_data)
    }
  }
  
  # Create separate dataframes for each type for plotting
  calo_future_all <- subset(all_future_distributions, 
                            type == "Calophyllum" & value > 0)
  vertlept_future_all <- subset(all_future_distributions, 
                                type == "Verticillium" & value > 0)
  overlap_future_all <- subset(all_future_distributions, 
                               type == "Overlap" & value > 0)
  
  # Create the facet plot of future distributions
  future_facet_plot <- ggplot() +
    # Background
    geom_raster(data = mask_df, aes(x = x, y = y), fill = bg_colour) +
    # Calophyllum distribution (green)
    geom_raster(data = calo_future_all, 
                aes(x = x, y = y), fill = "#006400", alpha = 0.6) +
    # Verticillium distribution (purple)
    geom_raster(data = vertlept_future_all, 
                aes(x = x, y = y), fill = "#6A0DAD", alpha = 0.6) +
    # Overlap (black)
    geom_raster(data = overlap_future_all, 
                aes(x = x, y = y), fill = "#000000") +
    # Facets
    facet_grid(scenario ~ time_period) +
    # Legend (invisible points for legend)
    geom_point(data = data.frame(
      species = factor(c("Calophyllum", "Verticillium", "Overlap"),
                       levels = c("Calophyllum", "Verticillium", "Overlap")),
      x = -Inf, y = -Inf,
      scenario = scenarios[1],
      time_period = time_periods[1]),
      aes(x = x, y = y, colour = species), alpha = 0) +
    scale_colour_manual(
      name = "Future Distributions",
      values = c("Calophyllum" = "#006400", 
                 "Verticillium" = "#6A0DAD", 
                 "Overlap" = "#000000"),
      labels = c(
        expression(italic("Calophyllum paniculatum")),
        expression(italic("Verticillium/Leptographium")),
        "Overlap"
      ),
      guide = guide_legend(override.aes = list(
        alpha = c(0.6, 0.6, 1),
        size = 5,
        shape = 15
      ))
    ) +
    labs(title = "Future Species Distributions") +
    get_map_theme() +
    coord_equal()
  
  # Create caption text with thresholds 
  caption_text <- paste(
    "Thresholds:",
    paste0("Calophyllum: ", get_threshold_description(threshold_calo)),
    paste0("Verticillium: ", get_threshold_description(threshold_vertlept)),
    sep = "  "
  )
  
  # Add caption
  future_facet_plot_with_caption <- future_facet_plot +
    labs(caption = caption_text) +
    theme(
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10))
    )
  
  # Save the facet plot
  filename <- paste0("future_distributions_facet_", 
                     threshold_calo, "_", threshold_vertlept, ".png")
  ggsave(
    file.path(output_path, filename),
    future_facet_plot_with_caption,
    width = 12,
    height = 12,
    dpi = 300
  )
  
  cat(paste0("Created future distribution facet plot: ", filename, "\n"), 
      file = log_file, append = TRUE)
  
  return(future_facet_plot_with_caption)
}
  
  
#==============================================================================#
#                                  Run ----
#==============================================================================#
  
# Record overall start time
overall_start_time <- Sys.time()
cat(paste("Overall visualisation process started at", overall_start_time, "\n"), 
    file = log_file, append = TRUE)
  
# Step 1: Generate range change plots
cat("\nGenerating range change visualisations...\n", file = log_file, append = TRUE)
range_change_plot <- create_range_change_plot(range_shifts, threshold_lookup)
save_plot(range_change_plot, "range_change_comparison", 
    file.path(figure_path, "range_change"), width = 12, height = 8)
  
# Step 2: Generate metrics visualisations
cat("\nGenerating metrics visualisations...\n", file = log_file, append = TRUE)
create_all_metrics_plots(range_shifts, threshold_lookup)
  
# Step 3: Generate summary plot
cat("\nGenerating summary visualisations...\n", file = log_file, append = TRUE)
summary_plot <- create_summary_plot(range_shifts)
  
# Step 4: Run bivariate analysis
cat("\nRunning bivariate distribution analysis...\n", file = log_file, append = TRUE)
bivariate_output <- run_bivariate_analysis(
    data05_out_cur,
    data05_out_pred,
    figure_path,
    threshold_pairs,
    scenarios,
    time_periods
  )
  
# Record overall completion
overall_end_time <- Sys.time()
overall_elapsed <- difftime(overall_end_time, overall_start_time, units = "mins")
  
# Finish log
cat("\n=============================================================================\n", 
    file = log_file, append = TRUE)
cat(paste("All visualisations completed at", Sys.time(), "\n"), 
    file = log_file, append = TRUE)
cat(paste0("Total processing time: ", round(overall_elapsed, 2), " minutes\n"), 
    file = log_file, append = TRUE)
cat(paste("Output visualisations saved to:", figure_path, "\n"), 
    file = log_file, append = TRUE)
cat("=============================================================================\n", 
    file = log_file, append = TRUE)
  
# Print completion message to console
cat("\n✓ visualisation script completed successfully!\n")
cat(paste0("Total processing time: ", round(overall_elapsed, 2), " minutes\n"))
cat("Output files saved to:", figure_path, "\n")
  
# Print a summary of outputs created
cat("\nvisualisation outputs summary:\n")
cat(paste0("- Range change plots: ", file.path(figure_path, "range_change"), "\n"))
cat(paste0("- Metrics visualisations: ", file.path(figure_path, "combined_metrics"), "\n"))
cat(paste0("- Bivariate distribution maps: ", file.path(figure_path, "bivariate"), "\n"))
cat(paste0("- Summary visualisations: ", figure_path, "\n"))# 07 visualisations - Optimized workflow
# Date: 2025-03-22
  
  
# Tidy up
rm(list = ls())
gc()


#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#
