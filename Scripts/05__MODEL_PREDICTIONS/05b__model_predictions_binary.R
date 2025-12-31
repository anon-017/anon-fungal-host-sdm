# ---
# title: "05b__model_predictions_binary"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-03-01"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Uses threshold values from 05__model_prediction > intermediate > thresholds to make binary versions of model predictions
#
# HOW
# - convert_to_binary function.
# 
# WHY
# - To make predictions of presence-absence rather than having a continuous habitat suitability grid output
# - Output binary predictions can be used to calculate range shift metrics
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#


## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

## Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

# Remove function that are not required for this script
rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, 
   align.forest.raster, align.raster, validate.processed.climate.files, calc_buffer_distance, calculate_metrics,
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, 
   random_bg_whole_area, rename_stack_layers, select07_cv, 
   analyze_elev_distribution, generate_buf_bg_pts, generate_dispersal_bg_pts, generate_elev_bg_pts, plot_bg_comparison)


## Set up base folders ----

# Required derived data for binary predictions:
## model files (top 1 per spp)
## future environmental data cropped to MDG

# Base folders
datafolder <- file.path("//xxx")

# Inputs: Occurrence data
# Data folder relating to script 00 (species occurrences)
datafolder00 <- file.path(datafolder, "00__species_occurrence_cleaning")
pres_list <- list.files(file.path(datafolder00, "output", "presence", "dupfree"), full.names = T)
load(pres_list[1]) # loads as "calonodupsp" (duplicate free presences used to train models)
load(pres_list[2]) # loads as "vert_lept_nodupsp" (duplicate free presences used to train models)
 
# Vectorise for plotting later
pres_calo <- terra::vect(calonodupsp, geom=c("x", "y"), crs="ESRI:102022")
pres_vertlept <- terra::vect(vert_lept_nodupsp, geom=c("x", "y"), crs="ESRI:102022")
rm(pres_list, calonodupsp, vert_lept_nodupsp) # keep tidy

# Inputs: Projections for current climate (05__ensemble_projection_feasibility.R)
datafolder05_curr <- file.path(datafolder, "05__model_prediction", "output", "current_projections","weighted_mean_ensembles", "extra_settings")
projections <- list.files(datafolder05_curr, pattern="current_prediction", full.names = T)

# Load the two current projection rasters
calo_proj <- rast(projections[1])
vertlept_proj <- rast(projections[2])

# Inputs: Predictions for future climates (05__model_prediction) for each species and year - loads nine model predictions per species (3 years x 3 scenarios each)
prediction_list <- list.files(file.path(datafolder, "05__model_prediction", "output", 
                                    "future_predictions", "extra_settings"), pattern="prediction", full.names = T)

# Inputs: Threshold values to make binary predictions
# 05__ ensemble prediction feasibility intermediate files
thresholds_dir <- file.path(datafolder, "05__model_prediction", "intermediate", "thresholds")
threshold_files <- list.files(thresholds_dir, full.names = T)

# Outputs: Binary prediction outputs
# binary output for current projections
data__05b_out_cur <- file.path(datafolder, "05__model_prediction", "output", "current_projections", "binary", "extra_settings")
dir.create(file.path(data__05b_out_cur), recursive = TRUE, showWarnings = FALSE)

# Outputs: Binary raster output for future predictions/forecasting
data__05b_out <- file.path(datafolder, "05__model_prediction", "output", "future_predictions", "binary", "extra_settings")
dir.create(file.path(data__05b_out), recursive = TRUE, showWarnings = FALSE)

# Outputs: Binary prediction figures
data__05b_out_fig <- file.path(datafolder, "05__model_prediction", "output", "figures", "binary", "extra_settings")
dir.create(file.path(data__05b_out_fig), recursive = TRUE, showWarnings = FALSE)


#==============================================================================#
#                        1. Load in required data ----
#==============================================================================#


# Species occurrences
species_list <- c("calo", "vertlept")

# Set up time spans
# 2040, 2070, 2100 - 3 time stamps to predict
# Define future time periods
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

# Set up the climate pathway versions
# SSP126, SSP370, SSP585
climate_pathways <- c("ssp126", "ssp370", "ssp585")

# Load all thresholds for each species
# Create a nested list to store thresholds by species and threshold type
threshold_values <- list()
threshold_types <- c("opt2", "opt4", "opt12")

# Loop through species to organize threshold values
for (species in species_list) {
  species_thresh_files <- grep(species, threshold_files, value = TRUE)
  threshold_values[[species]] <- list()
  
  # Match threshold files with threshold types
  for (thresh_type in threshold_types) {
    thresh_file <- grep(thresh_type, species_thresh_files, value = TRUE)
    if (length(thresh_file) > 0) {
      threshold_values[[species]][[thresh_type]] <- readRDS(thresh_file)
    }
  }
}


#==============================================================================#
#                        2. Helper functions ----
#==============================================================================#


# Function to extract scenario information from raster path
extract_scenario_info <- function(raster_path) {
  file_name <- basename(raster_path)
  
  # Extract species information
  species <- NA
  for (sp in species_list) {
    if (grepl(sp, file_name)) {
      species <- sp
      break
    }
  }
  # Extract climate scenario information
  is_current <- grepl("current_current", file_name)
  scenario <- NA
  time_period <- NA
  
  if (is_current) {
    scenario <- "current"
    time_period <- "current"
  } else {
    # Extract SSP scenario
    for (ssp in climate_pathways) {
      if (grepl(ssp, file_name)) {
        scenario <- ssp
        break
      }
    }
    # Extract time period
    for (period in time_periods) {
      if (grepl(period, file_name)) {
        time_period <- period
        break
      }
    }
  }
  return(list(
    path = raster_path,
    species = species,
    scenario = scenario,
    time_period = time_period,
    is_current = is_current
  ))
}

# Modified function to convert continuous prediction to binary based on threshold
convert_to_binary <- function(prediction_raster, threshold_value) {
  # Check if raster has values
  if (is.null(prediction_raster) || terra::global(prediction_raster, "isNA", na.rm=FALSE)$isNA == terra::ncell(prediction_raster)) {
    stop("Input raster has no valid values")
  }
  
  # Create a copy of the input raster
  binary_raster <- prediction_raster
  
  # Check for NA or missing values and ensure proper conversion
  # Use terra's reclassification which handles NAs better
  m <- c(-Inf, threshold_value, 0,
         threshold_value, Inf, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  
  # Apply reclassification
  binary_raster <- terra::classify(binary_raster, m, include.lowest=TRUE)
  
  return(binary_raster)
}


# Function to safely plot SpatRaster layers with custom palettes
# Addresses the .default.pal error
safe_plot <- function(rast_layer, title, pal_fun = hcl.colors, 
                      n_colors = 100, zlim = NULL, filename = NULL) {
  # Create color palette - using base R palettes to avoid dependency issues
  if (is.function(pal_fun)) {
    col_palette <- pal_fun(n_colors)
  } else if (pal_fun == "terrain") {
    col_palette <- terrain.colors(n_colors)
  } else if (pal_fun == "heat") {
    col_palette <- heat.colors(n_colors)
  } else if (pal_fun == "topo") {
    col_palette <- topo.colors(n_colors)
  } else if (pal_fun == "cm") {
    col_palette <- cm.colors(n_colors)
  } else {
    # Default to terrain colors if palette function not recognized
    col_palette <- terrain.colors(n_colors)
  }
  # Open device if filename is provided
  if (!is.null(filename)) {
    png(filename, width = 900, height = 650, res = 100)
  }
  
  # Plot with explicit parameters to avoid .default.pal error
  plot(rast_layer, 
       main = title,
       col = col_palette,
       zlim = zlim,
       legend = TRUE)
  
  # Close device if needed
  if (!is.null(filename)) {
    dev.off()
  }
}

# Function to create a standardized binary prediction plot with threshold info
create_binary_plot <- function(raster_data, species_name, scenario_info = NULL, 
                               presence_points, threshold_info, 
                               output_path, width = 1200, height = 1000) {
  
  # Extract threshold value, type, and description from the threshold_info parameter
  threshold_value <- threshold_info$value
  threshold_type <- threshold_info$type
  threshold_desc <- threshold_info$desc
  
  # Start PNG device with explicit parameters
  png(output_path, width = width, height = height, res = 120)
  
  # Create a fixed layout with two rows
  # Row 1: Title area (smaller)
  # Row 2: Plot area (larger)
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(1, 4))
  
  # --- PANEL 1: TITLE AREA ---
  par(mar = c(0, 0, 1, 0))
  plot.new()
  
  # Add species name in italics
  text(0.5, 0.7, species_name, font = 3, cex = 1.5)
  
  # Add scenario information
  if (is.null(scenario_info) || scenario_info$is_current) {
    scenario_text <- "Binary current projection (2010)"
  } else {
    scenario_text <- paste0(toupper(scenario_info$scenario), " ", scenario_info$time_period)
  }
  text(0.5, 0.4, scenario_text, cex = 1.3)
  
  # Add threshold information
  threshold_text <- paste0("Threshold: ", threshold_desc, " (", round(threshold_value, 4), ")")
  text(0.5, 0.1, threshold_text, cex = 1.1)
  
  # --- PANEL 2: MAP AREA ---
  # Set margins with extra space on right for legend
  par(mar = c(1, 1, 1, 8), xpd = TRUE)
  
  # Define colors based on species name
  if (grepl("vertlept|Verticillium/Leptographium", species_name, ignore.case = TRUE)) {
    # For vertlept species (purple theme)
    my_colors <- c("#E0E0E0", "#8A4D9E")  # Medium-light grey for unsuitable, medium purple for suitable
  } else {
    # For calo species (green theme)
    my_colors <- c("#E0E0E0", "#1E8E4C")  # Medium-light grey for unsuitable, medium green for suitable
  }
  # Plot the raster
  plot(raster_data,
       main = "",
       col = my_colors,
       legend = FALSE,
       axes = FALSE,
       box = FALSE)
  
  # Extract raster extent for legend positioning
  r_ext <- as.vector(ext(raster_data))
  
  # Position legend explicitly using coordinates
  # Place it outside the right edge of the plot
  legend(x = r_ext[2] + (r_ext[2] - r_ext[1]) * 0.01,  # Just outside right edge
         y = r_ext[4] - (r_ext[4] - r_ext[3]) * 0.1,   # Near the top
         legend = c("Unsuitable habitat", "Suitable habitat", "Occurrence records"),
         fill = c(my_colors, NA),
         pch = c(NA, NA, 21),
         pt.bg = c(NA, NA, "white"),
         col = c(NA, NA, "black"),
         pt.cex = c(NA, NA, 1.2),
         border = c("black", "black", NA),
         bty = "n",
         cex = 1.1)
  
  # Add presence points
  points(presence_points, pch = 21, bg = "white", col = "black", cex = 0.8)
  
  # Close device
  dev.off()
  
  return(TRUE)
}


process_prediction_raster <- function(prediction_file, layer_index = 1) {
  # Check if file exists
  if (!file.exists(prediction_file)) {
    stop(paste("Prediction file not found:", prediction_file))
  }
  
  # Load raster with try-catch to handle errors
  prediction_raster <- tryCatch({
    rast(prediction_file)
  }, error = function(e) {
    stop(paste("Error loading raster:", prediction_file, "-", e$message))
  })
  
  # Check if raster has data
  if (is.null(prediction_raster) || terra::global(prediction_raster, "isNA", na.rm=FALSE)$isNA == terra::ncell(prediction_raster)) {
    stop(paste("Loaded raster has no valid values:", prediction_file))
  }
  
  # Check if it's multi-layer and select first layer if needed
  if (nlyr(prediction_raster) > 1) {
    # For binary prediction, we use only the first layer (weighted prediction)
    return(prediction_raster[[layer_index]])
  } else {
    return(prediction_raster)
  }
}


# Function to get appropriate species display name
get_species_display_name <- function(species_code) {
  if (species_code == "calo") {
    return("Calophyllum paniculatum")
  } else if (species_code == "vertlept") {
    return("Verticillium/Leptographium")
  } else {
    return(species_code)
  }
}

# Function to get threshold description
get_threshold_description <- function(threshold_type) {
  switch(threshold_type,
         "opt2" = "Maximise Sum of Sensitivity and Specificity",
         "opt12" = "10th Percentile Training Presence",
         "opt4" = "Minimum Distance to ROC Curve",
         paste0("Unknown threshold (", threshold_type, ")"))
}



#==============================================================================#
#                        3. Run binary predictions ----
#==============================================================================#


## Process each species ----
for (species in species_list) {
  cat("\n\n=== Processing", species, "===\n")
  
  #### Get the projection rasters (current climate projections) ----
  if (species == "calo") {
    projection <- calo_proj
    projection_path <- projections[1]
    pres_points <- pres_calo
  } else {
    projection <- vertlept_proj
    projection_path <- projections[2]
    pres_points <- pres_vertlept
  }
  
  # Print some debug info about the projection raster
  cat("Projection raster summary:\n")
  cat("  - Path:", projection_path, "\n")
  cat("  - Number of layers:", nlyr(projection), "\n")
  
  #### Process each threshold type ----
  for (thresh_type in threshold_types) {
    # Skip if this threshold type doesn't exist for this species
    if (is.null(threshold_values[[species]][[thresh_type]])) {
      cat("Skipping threshold type", thresh_type, "for", species, "- not found\n")
      next
    }
    
    # Get threshold value for current species and threshold type
    threshold_value <- threshold_values[[species]][[thresh_type]]
    cat("Using threshold type:", thresh_type, "value:", threshold_value, "for", species, "\n")
    
    # Process current projection
    tryCatch({
      # Extract first layer from multi-layer raster
      proj_layer <- projection[[1]]
      cat("Selected layer 1 from raster\n")

      #### Convert current projection raster to binary ----
      # Create binary raster using terra's classify function
      cat("Converting to binary raster...\n")
      m <- matrix(c(-Inf, threshold_value, 0, 
                    threshold_value, Inf, 1), 
                  ncol=3, byrow=TRUE)
      
      proj_bin <- terra::classify(proj_layer, m, include.lowest=TRUE)
      cat("Successfully converted to binary raster\n")
      
      #### Write binary raster projection ----
      binary_current_file <- file.path(data__05b_out_cur, paste0(species, "_", thresh_type, "_bin_projection.tif"))
      cat("Writing binary raster to:", binary_current_file, "\n")
      writeRaster(proj_bin, filename = binary_current_file, overwrite = TRUE)
      cat("Successfully wrote binary raster\n")
      gc()
      
      #### Plot binary projection map ----
      # Create a quick plot for visual inspection of current binary projection
      figure_path <- file.path(data__05b_out_fig,
                               paste0(species, "_", thresh_type, "current_projection.png"))
      
      # Define species name for plotting
      species_name_str <- get_species_display_name(species)
      
      # Create threshold info object with description
      threshold_info <- list(
        value = threshold_value,
        type = thresh_type,
        desc = get_threshold_description(thresh_type)
      )
      
      cat("Creating binary plot for current projection...\n")
      # Create the plot using plot function "create_binary_plot"
      create_binary_plot(
        raster_data = proj_bin,
        species_name = species_name_str,
        scenario_info = list(is_current = TRUE),
        presence_points = pres_points,
        threshold_info = threshold_info,
        output_path = figure_path
      )
      
      cat("Figure saved to:", figure_path, "\n")
      
    }, error = function(e) {
      cat("ERROR processing current projection for", species, "with threshold", thresh_type, ":", e$message, "\n")
      # Continue with next threshold type
    })
    
    ### Binary pres/abs outputs for each prediction raster ----
    # Filter prediction list for species
    species_predictions <- grep(species, prediction_list, value = TRUE)
    
    cat("Found", length(species_predictions), "future prediction files for", species, "\n")
    if (length(species_predictions) == 0) {
      cat("Debug: All prediction files (first 5):\n")
      for (i in 1:min(5, length(prediction_list))) {
        cat("  -", basename(prediction_list[i]), "\n")
      }
    }
    
    for (prediction_file in species_predictions) {
      tryCatch({
        # Extract scenario information
        scenario_info <- extract_scenario_info(prediction_file)
        
        # Skip if this prediction is not for the current species or if scenario info is incomplete
        if (is.na(scenario_info$species) || scenario_info$species != species) {
          cat("Skipping", basename(prediction_file), "- not for", species, "\n")
          next
        }
        if (is.na(scenario_info$scenario) || is.na(scenario_info$time_period)) {
          cat("Skipping", basename(prediction_file), "- incomplete scenario information\n")
          next
        }
        cat("\n Creating Binary Prediction for:",
            "\n  Species:", species,
            "\n  Scenario:", scenario_info$scenario,
            "\n  Time period:", scenario_info$time_period,
            "\n  Threshold type:", thresh_type,
            "\n  Threshold value:", threshold_value, "\n")
        
        #### Load prediction raster ----
        cat("Loading prediction raster...\n")
        # Load and process the prediction raster
        prediction_raster <- terra::rast(prediction_file)
        prediction_layer <- prediction_raster[[1]]
        cat("Successfully loaded prediction raster (layer 1)\n")
        
        # Create output filename
        output_dir <- ifelse(scenario_info$is_current,
                             file.path(data__05b_out_cur),
                             file.path(data__05b_out))
        
        output_filename <- paste0(species, "_", scenario_info$scenario, "_",
                                  scenario_info$time_period, "_", thresh_type,"_bin.tif")
        output_predict_path <- file.path(output_dir, output_filename)
        
        #### Make binary prediction using thresholds per spp ----
        cat("Creating binary prediction...\n")
        # Use classification matrix
        m <- matrix(c(-Inf, threshold_value, 0, 
                      threshold_value, Inf, 1), 
                    ncol=3, byrow=TRUE)
        
        bin_predict <- terra::classify(prediction_layer, m, include.lowest=TRUE)
        cat("Successfully created binary prediction\n")
        
        #### Write binary prediction to file ----
        cat("Writing binary prediction to:", output_predict_path, "\n")
        writeRaster(bin_predict, filename = output_predict_path, overwrite = TRUE)
        cat("Successfully wrote binary prediction file\n")
        gc()
        
        #### Plot binary prediction ----
        figure_path <- file.path(data__05b_out_fig,
                                 paste0(species, "_", scenario_info$scenario, "_",
                                        scenario_info$time_period, thresh_type,"_bin.png"))
        
        # Define species name for plotting
        species_name_str <- get_species_display_name(species)
        
        # Create threshold info object
        threshold_info <- list(
          value = threshold_value,
          type = thresh_type,
          desc = get_threshold_description(thresh_type)
        )
        
        cat("Creating binary plot...\n")
        # Create the plot using our standard function
        create_binary_plot(
          raster_data = bin_predict,
          species_name = species_name_str,
          scenario_info = scenario_info,
          presence_points = pres_points,
          threshold_info = threshold_info,
          output_path = figure_path
        )
        
        cat("Figure saved to:", figure_path, "\n")
        
      }, error = function(e) {
        cat("ERROR processing future prediction", basename(prediction_file), "for", species, 
            "with threshold", thresh_type, ":", e$message, "\n")
        # Continue with next prediction file
      })
    }
  }
}
cat("\n=== Binary Prediction script completed successfully ===\n")


#==============================================================================#
#                        ----  End of workflow ----
#==============================================================================#

# # Clean up
gc()
rm(list = ls())