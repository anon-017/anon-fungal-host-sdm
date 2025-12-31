# ---
# title: "functions.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-03_09"
# updated: "2025-12_30"
# ---


#==============================================================================#
#                       LIST OF FUNCTIONS ----
#==============================================================================#

# General helper functions ----

# (for use in all scripts within project starting in 00__, 01__, 02__, 03__, 04__, 05__,06__)
# Function that installs and loads required packages
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos = 'http://cran.us.r-project.org')
  require(x, character.only = TRUE)
} # add names of the packages required as character objects
package_vec <- c(
  "rgbif", "foreach", "doParallel", "terra", "sf", "sfheaders", 
  "ggplot2", "readr", "lattice", "dplyr", 
  "gridExtra","corrplot", 
  "lattice", "cowplot" , "gdalUtilities","extrafont", 
 "stringr", "future.apply", "data.table", "sdm"
)
sapply(package_vec, install.load.package)


# (00__script specific) ----

## Coordinate cleaning ----
# Function to swap coordinates if they're in the wrong columns
swap_coords <- function(lat, lon) {
  if (lat >= lon_range[1] && lat <= lon_range[2] && 
      lon >= lat_range[1] && lon <= lat_range[2]) {
    return(c(lon, lat))
  } else {
    return(c(lat, lon))
  }
}


# spatial thinning (for use during 00__)
## Spatial thinning ----
thin <- function(sf, thin_dist = 2000, runs = 1, ncores = 1){
  
  require(sf, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(furrr, quietly = TRUE)
  
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  sf_buffer <- st_buffer(sf, thin_dist)
  buff_int <- st_intersects(sf, sf_buffer)
  buff_int <- setNames(buff_int, 1:length(buff_int))
  
  n_int <- map_dbl(buff_int, length)
  
  plan(multisession, workers = ncores)
  
  seeds <- sample.int(n = runs)
  results_runs <- future_map(seeds, function(i){
    
    set.seed(i)
    while (max(n_int) > 1) {
      max_neighbors <- names(which(n_int == max(n_int)))
      
      # remove point with max neighbors
      sampled_id <- sample.vec(max_neighbors, 1)
      
      pluck(buff_int, sampled_id) <- NULL
      buff_int <- map(buff_int, function(x) setdiff(x, as.numeric(sampled_id)))
      n_int <- map_dbl(buff_int, length)
    }
    
    unlist(buff_int) %>% unique()
    
  })
  
  lengths <- map_dbl(results_runs, length)
  
  selected_run <- results_runs[[sample.vec(which(lengths == max(lengths)), 1)]]
  
  out <- sf[selected_run,]
  
  out <- sf_to_df(out)[,3:4] %>%
    rename("lon" = "x", "lat" = "y")
  
  out
}

## Background points creation - CALO ----
### Ecoregion limitation ----
# Ecoregion-limited background points
generate_ecoregion_bg_pts <- function(presence_points, study_area, ecoregions, species_name,
                                      n_points = 1000) {
  # Extract eco IDs
  eco_ids <- unique(terra::extract(ecoregions, presence_points)[["ECO_ID"]])
  
  # Subset without aggregating
  background_area <- terra::intersect(
    ecoregions[ecoregions$ECO_ID %in% eco_ids,],
    study_area)
  
  # Sample points within ecoregion mask
  background_points <- terra::spatSample(background_area, size = n_points, method = "random")
  
  # Get coordinates of bg pts, then add to df
  coords <- terra::crds(background_points)
  bg_ecor_df <- data.frame(x = coords[,1], y = coords[,2])
  
  # Write points and buffer area
  base_name <- paste0(species_name, "_n", n_points)
  # Write vector shp
  writeVector(background_area, file.path(bg_out_dir, paste0(base_name, "_ecoregion_area.shp")), overwrite=TRUE)
  
  # Write df (RData)
  save(bg_ecor_df, file = file.path(bg_out_dir, paste0(base_name, "_ecor_bg_pts.RData")))
  
  # Write CSV file
  write.csv(bg_ecor_df[, c("x", "y")],  # Select only x,y coordinates for CSV
            file = file.path(bg_out_dir, paste0(base_name, "_ecor_bg_pts.csv")),
            row.names = FALSE)
  # Tidy up
  rm(background_area, coords, bg_ecor_df, base_name, eco_ids)
  gc()
  return(background_points)
}

# # # THIS FUNCTION WOULD NOT PLOT AFR - TOO LARGE/DO IN GIS
#### Plot to check bg point methods ----
# plot_bg_comparison <- function(species_points, study_area, 
#                                bg_buffer, bg_elev, bg_disper, bg_eco,
#                                species_name) {
#   
#   # start png device for writing to file
#   png(file.path(plot_out_dir, paste0(species_name, "_bg_comparison.png")),
#       width = 8, height = 8, units = "in", res = 300)
#   
#   # Set up a 2x2 plotting layout with minimal margins
#   par(mfrow = c(2, 2),
#       # Reduce all margins (bottom, left, top, right) for each plot
#       mar = c(0, 0, 1.5, 0),
#       # Reduce outer margins
#       oma = c(0, 0, 2, 0))
#   
#   # Enhanced base plotting function with axis removal
#   plot_base <- function(main) {
#     plot(study_area, 
#          main = main,
#          col = "#F5F5F5", 
#          border = "grey",
#          axes = FALSE,  # No axes
#          xlab = "",    # No x label
#          ylab = "")    # No y label
#     plot(species_points, 
#          add = TRUE, 
#          col = "black", 
#          pch = 16, 
#          cex = 0.5)
#   }
#   
#   # Plot each background point method
#   plot_base("Buffer constraint (ratio to study area)")
#   plot(bg_buffer, add = TRUE, col = alpha("navy", 0.3), border = NA)
#   
#   plot_base("Elevation constraint")
#   plot(bg_elev, add = TRUE, col = alpha("darkgreen", 0.3), border = NA)
#   
#   plot_base("Dispersal limitation")
#   plot(bg_disper, add = TRUE, col = alpha("purple", 0.3), border = NA)
#   
#   plot_base("Eco-region constraint")
#   plot(bg_eco, add = TRUE, col = alpha("orange", 0.3), border = NA)
#   
#   # Add overall title if species name is provided
#   if (!missing(species_name)) {
#     mtext(species_name, outer = TRUE, line = 0.5, cex = 1.2)
#   }
#   dev.off() # close png device
# }


####

## Background points creation - VERT-LEPT ----
### 2° far method limitation ----

# STEPS
# 1. Convert 2 degrees into a geographical distance (km/m) based on latitudinal position
# 2. Function that finds clusters of species presences and assigns differing buffer values
# based on their latitudinal extent
# 3. Erases those buffer areas from study extent
# 4 Generates background points from within the "cookie-cut" area left 
# (areas within study extent, but outside of 2deg buffers)

# Function to calculate / convert 
# similar to : https://www.opendem.info/arc2meters.html 

# The required inputs are:
# latitude value in decimal degrees
calc_buffer_distance <- function(lat) {
  # Earth's radius in meters
  R <- 6371000
  
  # Calculate the length of 2 degrees of latitude at this position
  # Using haversine formula simplified for north-south distance
  buffer_dist <- R * (2 * pi/180) # 2 degrees converts to radians with 2 * pi/180
  return(buffer_dist) # returns a buffer distance of approx. 222.6 km
  }


process_2degfar_bg_pts <- function(points, study_area, n_clusters = 3, n_random = 340) {
  # Store the target CRS (from study_area)
  target_crs <- crs(study_area)
  
  # If points are in a different CRS, project to study_area
  if (crs(points) != target_crs) {
    points <- project(points, target_crs)
  }
  # Extract coordinates and latitudes - must be geographic coordinates for latitude calcs
  points_geo <- project(points, "EPSG:4326")
  coords <- crds(points_geo)
  lats <- coords[, 2]
  
  # Perform k-means clustering on latitudes
  kmeans_result <- kmeans(lats, centers = n_clusters)
  
  # Add cluster information to points (in original projection)
  points$cluster <- kmeans_result$cluster
  
  # Initialize all_buffers with the first cluster
  first_cluster <- points[points$cluster == 1,]
  first_cluster_geo <- points_geo[points_geo$cluster == 1,]
  mean_lat <- mean(crds(first_cluster_geo)[, 2])
  buffer_dist <- calc_buffer_distance(mean_lat)
  all_buffers <- buffer(first_cluster, width = buffer_dist)
  
  # Process remaining clusters
  if(n_clusters > 1) {
    for(i in 2:n_clusters) {
      cluster_points <- points[points$cluster == i,]
      cluster_points_geo <- points_geo[points_geo$cluster == i,]
      mean_lat <- mean(crds(cluster_points_geo)[, 2])
      buffer_dist <- calc_buffer_distance(mean_lat)
      current_buffer <- buffer(cluster_points, width = buffer_dist)
      
      # Merge with existing buffers
      all_buffers <- terra::union(all_buffers, current_buffer)
    }
  }

  # Mask buffers to study area (crops to land - no buffer overlapping coastal regions)
  all_buffers <- terra::intersect(all_buffers, study_area)
  
  # Get coordinates of bg pts, then add to df
  coords <- terra::crds(random_points)
  bg_2far_df <- data.frame(x = coords[,1], y = coords[,2])
  
  # Write points (RData & csv) and buffer area (shp)
  base_name <- paste0("vert_lept_2far", "_n", n_random)
  # Write vector shp
  writeVector(remaining_area, file.path(bg_out_dir, paste0(base_name, ".shp")), overwrite=TRUE)
  
  # Write df (RData)
  save(bg_2far_df, file = file.path(bg_out_dir, paste0(base_name, "_bg_pts.RData")))
  
  # Write CSV file
  write.csv(bg_2far_df[, c("x", "y")],  # Select only x,y coordinates for CSV
            file = file.path(bg_out_dir, paste0(base_name, "_bg_pts.csv")),
            row.names = FALSE)
  
  # Return results
  return(list(
    original_points = points,
    clustered_points = points,
    buffers = all_buffers,
    remaining_area = remaining_area,
    random_points = random_points
  ))
}

  
  
####

# (01__script specific) ----
## Validate TIFF file for climate data ----
validate.tiff <- function(file_path) {
  tryCatch({
    gdalinfo(file_path)
    return(TRUE)
  }, error = function(e) {
    write_log(paste("TIFF validation failed for:", file_path))
    write_log(paste("Error message:", e$message))
    return(FALSE)
  })
}

## Repair TIFF file for corrupted climate data ----
repair.tiff <- function(input_file, output_file) {
  tryCatch({
    gdal_translate(input_file, output_file, of = "GTiff", co = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    return(TRUE)
  }, error = function(e) {
    write_log(paste("TIFF repair failed for:", input_file))
    write_log(paste("Error message:", e$message))
    return(FALSE)
  })
}

####

# (01__) 
## Align rasters ----
# Function to process a single global climate raster file and align it with chosen tif file (currentforest - forestatrisk)
align.raster <- function(input_file, output_file) {
  tryCatch({
    # Read the input raster
    input_raster <- rast(input_file)
    
    # Create a temporary file for the projected raster to not to use too much memory
    temp_projected_file <- file.path(temp_dir, paste0("temp_projected_", basename(input_file)))
    
    # Check if the projected file already exists so not processed twice
    if (file.exists(temp_projected_file)) {
      write_log(paste("Using existing projected file:", temp_projected_file))
      projected_raster <- rast(temp_projected_file)
    } else {
      # Project to match current_forest CRS 
      projected_raster <- terra::project(input_raster, current_forest, filename = temp_projected_file)
    }
    
    # Crop to current_forest extent (this is going from global to MDG extent)
    cropped_raster <- crop(projected_raster, current_forest, mask=TRUE)
    
    # Ensure exact same extent as current_forest
    combined_extent <- ext(
      min(ext(cropped_raster)[1], ext(current_forest)[1]),
      max(ext(cropped_raster)[2], ext(current_forest)[2]),
      min(ext(cropped_raster)[3], ext(current_forest)[3]),
      max(ext(cropped_raster)[4], ext(current_forest)[4])
    )
    
    # Extend the cropped raster to the combined extent
    extended_raster <- extend(cropped_raster, combined_extent)
    
    # Resample to match current_forest resolution - SDMs require exact same cell position and resolution
    final_raster <- resample(extended_raster, current_forest, method = 'bilinear')
    
    # Write the processed raster to file
    writeRaster(final_raster, filename = output_file, overwrite = TRUE)
    
    # Remove temporary projected file to clean up
    if (file.exists(temp_projected_file)) {
      file.remove(temp_projected_file)
    }
    
    return(TRUE)
  }, error = function(e) {
    write_log(paste("Error processing file:", input_file))
    write_log(paste("Error message:", e$message))
    return(FALSE)
  })
}


## Align forest rasters ----
# Function to align forest rasters with current_forest raster 4326 (originally African extent and projection)
align.forest.raster <- function(input_file, output_file, reference_raster, temp_dir) {
  tryCatch({
    # Input validation
    if (!file.exists(input_file)) stop("Input file does not exist")
    if (!dir.exists(dirname(output_file))) stop("Output directory does not exist")
    
    # Read the input raster
    input_raster <- rast(input_file)
    
    # Create a temporary file for the projected raster
    temp_projected_file <- file.path(temp_dir, paste0("temp_projected_", basename(input_file)))
    
    # Check if the projected file already exists
    if (file.exists(temp_projected_file)) {
      message(paste("Using existing projected file:", temp_projected_file))
      projected_raster <- rast(temp_projected_file)
    } else {
      # Project to match reference_raster CRS 
      projected_raster <- terra::project(input_raster, reference_raster, filename = temp_projected_file)
    }
    
    # Crop to reference_raster extent and resample in one step
    final_raster <- terra::crop(projected_raster, reference_raster, mask=TRUE, snap="near")
    final_raster <- terra::resample(final_raster, reference_raster, method = 'bilinear')
    
    # Write the processed raster to file
    terra::writeRaster(final_raster, filename = output_file, overwrite = TRUE)
    
    # Remove temporary projected file
    if (file.exists(temp_projected_file)) {
      file.remove(temp_projected_file)
    }
    
    return(TRUE)
  }, error = function(e) {
    message(paste("Error processing file:", input_file))
    message(paste("Error message:", e$message))
    return(FALSE)
  })
}



## Process bioclim ----
# Function to process climate data (for use within 01__)
process.climate.data <- function(scenario, year_range) {
  # base path to where the climate data is stored, first by year range, then in each sub folder by scenario
  base_path <- here("CHELSA", "RAW", "futures", year_range, paste0("SSP", scenario))
  
  # Load climate data
  climate_files <- list.files(base_path, pattern = "\\.tif$", full.names = TRUE)
  climate_files <- climate_files[!grepl("fcc_", climate_files)]  # Exclude forest cover files 
  
  # Process climate data file by file
  for (file in climate_files) {
    # Extract bioclimatic number from OG file name and create simplified updated filename for output later
    bio_num <- gsub(".*bio(\\d+).*", "\\1", basename(file))
    output_filename <- paste0("bio", bio_num, "_", year_range, "_ssp", scenario, ".tif")
    output_path <- file.path(output_dir, output_filename)
    
    # Check if the file has already been processed
    if (file.exists(output_path)) {
      cat("Skipping (already processed):", output_filename, "\n")
      next
    }
    
    cat("Processing:", output_filename, "\n") # adding for visual updates while processing
    
    # Validate TIFF file in case corrupted while downloading original climate data
    if (!validate.tiff(file)) {
      # Attempt to repair the file
      repaired_file <- file.path(temp_dir, paste0("repaired_", basename(file)))
      if (repair.tiff(file, repaired_file)) {
        file <- repaired_file
      } else {
        write_log(paste("Skipping file due to validation and repair failure:", file))
        next
      }
    }
    
    # Process the file if successfully validated and/or repaired
    success <- align.raster(file, output_path)
    
    if (!success) {
      write_log(paste("Failed to process file:", file))
    } else {
      write_log(paste("Successfully processed file:", file)) # update log file for each climate file 
    }
    
    # Clear memory to save space for next file
    gc()
  }
  
  return(output_dir)
}

## Validate processed climate data ----
validate.processed.climate.files <- function(output_dir, current_forest) {
  # Get list of all processed files
  processed_files <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # Initialize counters
  total_files <- length(processed_files)
  successful_files <- 0
  aligned_files <- 0
  
  # Create a log file for validation results
  validation_log <- file.path(output_dir, "validation_log.txt")
  write_validation_log <- function(message) {
    cat(paste0(Sys.time(), " - ", message, "\n"), file = validation_log, append = TRUE)
  }
  
  write_validation_log(paste("Validating", total_files, "processed files"))
  
  for (file in processed_files) {
    tryCatch({
      # Try to read the file
      r <- rast(file)
      successful_files <- successful_files + 1
      
      # Check alignment with current forest layer
      if (identical(ext(r), ext(current_forest)) &&
          identical(res(r), res(current_forest)) &&
          identical(crs(r), crs(current_forest))) {
        aligned_files <- aligned_files + 1
        write_validation_log(paste("File aligned correctly:", basename(file)))
      } else {
        write_validation_log(paste("File not aligned:", basename(file)))
        if (!identical(ext(r), ext(current_forest))) {
          write_validation_log("  Extent mismatch")
        }
        if (!identical(res(r), res(current_forest))) {
          write_validation_log("  Resolution mismatch")
        }
        if (!identical(crs(r), crs(current_forest))) {
          write_validation_log("  CRS mismatch")
        }
      }
    }, error = function(e) {
      write_validation_log(paste("Error reading file:", basename(file)))
      write_validation_log(paste("Error message:", e$message))
    })
  }
  
  # Calculate percentages
  percent_successful <- (successful_files / total_files) * 100
  percent_aligned <- (aligned_files / total_files) * 100
  
  # Write summary to log
  write_validation_log("\nValidation Summary:")
  write_validation_log(paste("Total files:", total_files))
  write_validation_log(paste("Successfully processed files:", successful_files, sprintf("(%.2f%%)", percent_successful)))
  write_validation_log(paste("Correctly aligned files:", aligned_files, sprintf("(%.2f%%)", percent_aligned)))
  
  # Print summary to console
  cat("\nValidation Summary:\n")
  cat(paste("Total files:", total_files, "\n"))
  cat(paste("Successfully processed files:", successful_files, sprintf("(%.2f%%)\n", percent_successful)))
  cat(paste("Correctly aligned files:", aligned_files, sprintf("(%.2f%%)\n", percent_aligned)))
  
  # Return results as a list
  return(list(
    total_files = total_files,
    successful_files = successful_files,
    aligned_files = aligned_files,
    percent_successful = percent_successful,
    percent_aligned = percent_aligned
  ))
}

# (03__script specific) ----


## Ensemble SDM (SDM package) ----
create_sdm_ensemble <- function(species_name, presence_points, env_data, 
                                bg_points_list, n_replicates = 10) {
  
  ensemble_results <- list()
  performance_metrics <- list()
  
  for(bg_set in 1:length(bg_points_list)) {
    bg_points <- bg_points_list[[bg_set]]
    replicate_models <- list()
    replicate_metrics <- list()
    
    for(rep in 1:n_replicates) {
      # Data partitioning
      train_indices <- createDataPartition(presence_points$presence, p = 0.7, list = FALSE)
      train_data <- presence_points[train_indices,]
      test_data <- presence_points[-train_indices,]
      
      # Convert to terra format
      train_spatial <- vect(train_data, geom = c("x", "y"))
      test_spatial <- vect(test_data, geom = c("x", "y"))
      bg_spatial <- vect(bg_points, geom = c("x", "y"))
      
      # Create sdmData object with terra objects
      sdm_data <- sdmData(presence ~ ., train_spatial, 
                          predictors = env_data,
                          bg = bg_spatial)
      
      # Fit models
      sdm_models <- sdm(presence ~ .,
                        data = sdm_data,
                        methods = c('rf', 'brt', 'maxent'),
                        replication = 'cv',
                        cv.folds = 5)
      
      # Calculate metrics for each algorithm
      model_metrics <- list()
      for(method in c('rf', 'brt', 'maxent')) {
        pred <- predict(sdm_models, test_spatial, method)
        metrics <- calculate_metrics(test_data$presence, pred)
        model_metrics[[method]] <- metrics
      }
      
      replicate_models[[rep]] <- sdm_models
      replicate_metrics[[rep]] <- model_metrics
    }
    
    # Create ensemble
    ensemble_model <- ensemble(replicate_models,
                               setting = list(method = 'weighted',
                                              stat = 'AUC',
                                              weight = TRUE))
    
    # Calculate ensemble metrics
    ensemble_pred <- predict(ensemble_model, test_spatial)
    ensemble_metrics <- calculate_metrics(test_data$presence, ensemble_pred)
    
    # Store results
    ensemble_results[[bg_set]] <- ensemble_model
    performance_metrics[[bg_set]] <- list(
      individual_models = replicate_metrics,
      ensemble = ensemble_metrics
    )
  }
  
  return(list(
    models = ensemble_results,
    metrics = performance_metrics
  ))
}



## SDM performance metrics (SDM package) ----

# Function to calculate SDM performance metrics
calculate_metrics <- function(observed, predicted) {
  # Calculate confusion matrix
  cm <- confusionMatrix(factor(predicted > 0.5), factor(observed))
  
  # Calculate additional metrics
  thresh <- threshold(observed, predicted)
  boyce <- ecospat.boyce(observed, predicted)
  
  metrics <- list(
    AUC = auc(observed, predicted),
    TSS = max(sensitivity(observed, predicted) + specificity(observed, predicted) - 1),
    Kappa = cm$overall["Kappa"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Threshold = thresh$threshold,
    Boyce = boyce$cor
  )
  return(metrics)
}







# ## RF pdps ----
# # Function to create and save a partial dependence plot (random forest specific)
# create.rf.pdp <- function(var_name) {
#   pdp_data <- partial(calo_rf_cv, pred.var = var_name, prob = TRUE, which.class = "Class1")
# 
#   # Ensure column names are valid for ggplot
#   names(pdp_data) <- make.names(names(pdp_data))
# 
#   # Get the full name of the variable
#   full_name <- env_vars_full[var_name]
#   if (is.na(full_name)) full_name <- var_name  # Use original name if not in env_vars_full
# 
#   plot <- ggplot(pdp_data, aes(x = .data[[make.names(var_name)]], y = yhat)) +
#     geom_line() +
#     geom_rug(sides = "b", alpha = 0.2) +
#     theme_minimal() +
#     labs(
#       #title = paste("Partial Dependence Plot for", full_name),
#       x = full_name,
#       y = "Partial dependence"
#     )+
#     theme(text = element_text(family = "Times New Roman", size = 14)
#           )
# 
#   return(plot)
# }
# 
# 
# ## BRT PDPs ----
# 
# # Function to plot BRT partial dependence and save to PDF with debugging
# 
# create.brt.pdp <- function(var_name, model) {
#   # Generate partial dependence data for the boosted regression tree model
#   pdp_data <- partial(
#     calo_brt_cv,
#     pred.var = var_name,
#     prob = TRUE,
#     which.class = "Class1",  # For classification problems
#     plot = FALSE,            # Ensures only data is returned, not a plot
#     rug = TRUE               # If you want to include rug plots
#   )
# 
#   # Ensure column names are valid for ggplot
#   names(pdp_data) <- make.names(names(pdp_data))
# 
#   # Get the full name of the variable
#   full_name <- env_vars_full[var_name]
#   if (is.na(full_name)) full_name <- var_name  # Use original name if not in env_vars_full
# 
#   # Create the plot with the specified theme
#   plot <- ggplot(pdp_data, aes(x = .data[[make.names(var_name)]], y = yhat)) +
#     geom_line() +
#     geom_rug(sides = "b", alpha = 0.2) +
#     theme_minimal() +
#     labs(
#       x = full_name,
#       y = "Partial dependence"
#     ) +
#     theme(
#       text = element_text(family = "Times New Roman", size = 14)
#     )
# 
#   return(plot)
# }
# 
# 
# 
# ## Maxent PDP ----
# 
# # Function to create a set of Maxent PDP plots
# 
# create.mxt.pdp <- function(var_name, model, train_data) {
#   # Generate partial dependence data for the MaxEnt model
#   mxt_pdp_data <- partial(
#     model,
#     pred.var = var_name,
#     train = train_data,          # Explicitly pass the training data
#     prob = TRUE,
#     which.class = "presence",  # For MaxEnt, it is named "presence" for species distribution
#     plot = FALSE,              # Ensures only data is returned, not a plot
#     rug = TRUE                 # If you want to include rug plots
#   )
# 
#   # Ensure column names are valid for ggplot
#   names(mxt_pdp_data) <- make.names(names(mxt_pdp_data))
# 
#   # Get the full name of the variable
#   full_name <- env_vars_full[var_name]
#   if (is.na(full_name)) full_name <- var_name  # Use original name if not in env_vars_full
# 
#   # Create the plot with the specified theme
#   plot <- ggplot(mxt_pdp_data, aes(x = .data[[make.names(var_name)]], y = yhat)) +
#     geom_line() +
#     geom_rug(sides = "b", alpha = 0.2) +
#     theme_minimal() +
#     labs(
#       x = full_name,
#       y = "Partial dependence"
#     ) +
#     theme(
#       text = element_text(family = "Times New Roman", size = 14)
#     )
# 
#   return(plot)
# }
# 
# ## Maxent manual PDP ----
# 
# # Function: Manual PDF creation as MaxEnt algorithm does not handle predictions the same way as RF or BRT
# 
# create.mxt.pdp.manual <- function(var_name, model, train_data, grid_resolution = 100) {
#   # Check if the model and training data are valid
#   if (missing(model) || missing(train_data)) {
#     stop("Model and training data must be provided.")
#   }
# 
#   # Create a sequence of values for the variable of interest
#   grid_vals <- seq(min(train_data[[var_name]], na.rm = TRUE),
#                    max(train_data[[var_name]], na.rm = TRUE),
#                    length.out = grid_resolution)
# 
#   # Create a df where the variable of interest varies, but others are held constant
#   mxt_pdp_data <- train_data[1, , drop = FALSE]  # Use the first row as a template, keep it as a df
#   mxt_pdp_data <- mxt_pdp_data[rep(1, grid_resolution), ]  # Repeat the template row
#   mxt_pdp_data[[var_name]] <- grid_vals  # Replace the variable of interest with the grid values
# 
#   # Predict using the MaxEnt model
#   predictions <- predict(model, mxt_pdp_data, type = "cloglog")
# 
#   # Check if the predictions are a matrix or vector
#   if (is.matrix(predictions)) {
#     yhat <- predictions[, 1]  # Extract the first column if it's a matrix
#   } else {
#     yhat <- predictions  # If it's a vector, use it directly
#   }
# 
#   # Add the predictions to the df
#   mxt_pdp_data$yhat <- yhat
# 
#   # Return the PDP data
#   return(mxt_pdp_data)
# }
# 
# 
# ## Plot ROC curve ----
# 
# # Utility function to plot ROC curve
# plot.roc.curves <- function(curves, models) {
#   valid_curves <- !sapply(curves, is.null)
#   curves <- curves[valid_curves]
#   models <- models[valid_curves]
# 
#   # Create a data frame with all ROC curves
#   roc_data <- do.call(rbind, lapply(seq_along(curves), function(i) {
#     data.frame(
#       FPR = 1 - curves[[i]]$specificities,
#       TPR = curves[[i]]$sensitivities,
#       Model = models[i]
#     )
#   }))
# 
#   # Convert Model to a factor with levels in the desired order
#   roc_data$Model <- factor(roc_data$Model, levels = models)
# 
#   ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#     geom_line() +
#     labs(x = "False Positive Rate",
#          y = "True Positive Rate") +
#     theme_minimal(base_family = "Times New Roman", base_size = 14) +
#     scale_color_brewer(palette = "Dark2") +
#     theme(legend.position = "bottom")
# }
# 
# ## Plot Precision-Recall curve ----
# 
# # Utility function to plot Precision-Recall curve
# plot.pr.curves <- function(metrics_list, models) {
#   # Function to calculate precision and recall at various thresholds
#   calc_pr_curve <- function(metrics) {
#     pred <- prediction(metrics$probabilities, metrics$true_labels)
#     perf <- performance(pred, "prec", "rec")
#     df <- data.frame(
#       recall = perf@x.values[[1]],
#       precision = perf@y.values[[1]]
#     )
#     df <- df[!is.na(df$precision) & !is.na(df$recall), ]  # Remove NA values
#     return(df)
#   }
# 
#   # Create a data frame with all PR curves
#   pr_data <- do.call(rbind, lapply(seq_along(metrics_list), function(i) {
#     pr <- calc_pr_curve(metrics_list[[i]])
#     pr$Model <- models[i]
#     return(pr)
#   }))
# 
#   # Convert Model to a factor with levels in the desired order of models in legend
#   pr_data$Model <- factor(pr_data$Model, levels = models)
# 
#   # Create the plot
#   ggplot(pr_data, aes(x = recall, y = precision, color = Model)) +
#     geom_line() +
#     labs(x = "Recall",
#          y = "Precision") +
#     theme_minimal(base_family = "Times New Roman", base_size = 14) +
#     scale_color_brewer(palette = "Dark2") +
#     theme(legend.position = "bottom")  # Remove legend
# }
# 
# 
# 
# 
# ## Summarize and plot metrics ----
# summarize.and.plot.metrics <- function(metrics_df, model_abbreviations) {
#   # Remove any row names that might be causing issues
#   rownames(metrics_df) <- NULL
# 
#   # Create a named vector for mapping full names to abbreviations
#   full_names <- unique(metrics_df$Model)
#   name_to_abbrev <- setNames(model_abbreviations, full_names)
# 
#   # Replace full model names with abbreviations
#   metrics_df$Model <- name_to_abbrev[metrics_df$Model]
# 
#   # Ensure the Model column is a factor with the specified order
#   metrics_df$Model <- factor(metrics_df$Model, levels = model_abbreviations)
# 
#   # Melt the dataframe to long format
#   metrics_long <- melt(metrics_df, id.vars = "Model", variable.name = "Metric", value.name = "Value")
# 
#   # Remove any rows with NA values
#   metrics_long <- na.omit(metrics_long)
# 
#   # Function to create a plot for a single metric
#   create_metric_plot <- function(metric_data, metric_name) {
#     ggplot(metric_data, aes(x = Model, y = Value, fill = Model)) +
#       geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
#       labs(x = NULL, y = NULL) +
#       ggtitle(metric_name) +
#       theme_minimal() +
#       theme(
#         text = element_text(family = "Times New Roman", size = 12),
#         axis.text = element_text(size = 10),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
#         legend.position = "none",
#         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
#       ) +
#       scale_fill_brewer(palette = "Set2") +
#       coord_cartesian(ylim = c(0, 1))
#   }
# 
#   # Create a list of plots, one for each metric
#   metric_names <- setdiff(colnames(metrics_df), "Model")
#   plot_list <- lapply(metric_names, function(metric) {
#     metric_data <- metrics_long %>% filter(Metric == metric)
#     create_metric_plot(metric_data, metric)
#   })
# 
#   # Arrange the plots in a grid
#   n_metrics <- length(metric_names)
#   n_cols <- min(3, n_metrics)
#   n_rows <- ceiling(n_metrics / n_cols)
# 
#   combined_plot <- plot_grid(plotlist = plot_list, ncol = n_cols, nrow = n_rows)
# 
#   return(combined_plot)
# }
# 
# 
# ## Create metrics df ----
# create.metrics.df <- function(predictions, probabilities, true_labels, model_name) {
#   # Ensure factors are properly set
#   predictions <- factor(predictions, levels = c("Class0", "Class1"))
#   true_labels <- factor(true_labels, levels = c("Class0", "Class1"))
# 
#   # Calculate metrics
#   cm <- confusionMatrix(predictions, true_labels)
#   roc_obj <- roc(true_labels, probabilities)
#   kappa_tss <- calc_kappa_tss(true_labels, probabilities)
# 
#   metrics <- data.frame(
#     Model = model_name,
#     AUC = as.numeric(auc(roc_obj)),
#     Accuracy = cm$overall["Accuracy"],
#     Kappa = kappa_tss$kappa,
#     Sensitivity = cm$byClass["Sensitivity"],
#     Specificity = cm$byClass["Specificity"],
#     TSS = kappa_tss$tss,
#     Precision = cm$byClass["Pos Pred Value"],
#     Recall = cm$byClass["Sensitivity"],
#     F1_Score = cm$byClass["F1"]
#   )
# 
#   return(list(metrics = metrics, roc = roc_obj))
# }
# 
# 
# ## Calculate Kappa and TSS
# calc.kappa.tss <- function(true_labels, probabilities) {
#   predictions <- factor(as.numeric(probabilities > 0.5), levels = c(0, 1), labels = c("Class0", "Class1"))
#   true_labels <- factor(true_labels, levels = c("Class0", "Class1"))
# 
#   cm <- confusionMatrix(predictions, true_labels)
#   kappa <- cm$overall["Kappa"]
# 
#   sensitivity <- cm$byClass["Sensitivity"]
#   specificity <- cm$byClass["Specificity"]
#   tss <- sensitivity + specificity - 1
# 
#   return(list(kappa = kappa, tss = tss))
# }

# (04__script specific) ----




# (05__ script specific) ----



# (06__ script specific) ----



# End of functions ----