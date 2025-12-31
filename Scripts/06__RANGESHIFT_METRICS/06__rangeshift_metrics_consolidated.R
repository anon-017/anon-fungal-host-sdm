# ---
# title: "06__range_shift_metrics_consolidated.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-03-04"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Calculates range shift metrics between all current and future climate pairs of predictions
# - Based on binarised predictions calculated in previous 05__ codes
#
# HOW
# - Habitat distance
#  Ref (Choe et al., 2017)
#  Median distance between the centroid of the current range and all points in the future range. 
#  It measures how far, on average, the species will need to disperse to reach suitable future habitat@
#  median_distance <- calculate_habitat_distance(current_binary, future_binary)
# - Habitat exposure
#  Ref (Choe et al., 2017)
#  calculated as the difference between the current area range and the area of intersection between current 
#  and future range, divided by current range. 
#  It determines how different the future climate conditions are compared to those in the current projected 
# habitat
# - Spatial disruption
#  Ref (Yesuf et al., 2021)
#  Measures the relative change in range size compared to the total area ever occupied:
#  disruption <- calculate_spatial_disruption(current_binary, future_binary)
# - IUCN PAOO
#  IUCN Predicted Area of Occupancy (PAOO)
#  Follows IUCN guidelines for calculating Area of Occupancy using standardized 2×2 km grid cells:
#  paoo <- calculate_iucn_paoo(current_binary, future_binary, grid_size = 2)
# 
# WHY
# - Quantify risk to Calophyllum into the future: habitat loss/climate suitability loss?
# - Area overlap with Verticillium/Leptographium 
# - Find out the most interesting/worrying patterns from all results
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                         0. Workspace set up ----
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


#==============================================================================#
#                           1. Directory setup ----                                  
#==============================================================================#

## Set up base folders ----

# Base folder
datafolder <- file.path("//xxx")

# Input directories for binary projections
binary_cur_dir <- file.path(datafolder, "05__model_prediction", "output", "current_projections", "binary", "extra_settings")
binary_fut_dir <- file.path(datafolder, "05__model_prediction", "output", "future_predictions", "binary", "extra_settings")

# Template and elevation rasters
template_raster <- terra::rast(file.path(datafolder, "06__rangeshift_metrics", "input", "AOI", "MDG_aea.tif"))
elev <- terra::rast(file.path(datafolder, "06__rangeshift_metrics", "input", "elevation", "MDG_elv_msk.tif"))
elev_MDG <- terra::project(elev, template_raster) |> terra::resample(template_raster, method="bilinear")

# Output directory
data06_out <- file.path(datafolder, "06__rangeshift_metrics")
#dir.create(file.path(data06_out, "output/metrics/extra_settings/summary/"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(data06_out, "output/OikosMS"), recursive = TRUE, showWarnings = FALSE)


#==============================================================================#
#                             2.  Parameter Setup ----                                    
#==============================================================================#

# Species information
species_list <- c("calo", "vertlept")
species_fullnames <- c("C. paniculatum", "L. calophylli")
names(species_fullnames) <- species_list

# Time periods, climate pathways, thresholds
time_periods <- c("2011-2040", "2041-2070", "2071-2100")
climate_pathways <- c("SSP126", "SSP370", "SSP585")
#threshold_types <- c("opt2", "opt4","opt12")
threshold_types <- c("opt2", "opt4")


#==============================================================================#
#                         3. Helper Functions ----                                     
#==============================================================================#

# Extracts metadata from raster filenames
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
  
  # Extract threshold information - NEW
  threshold <- NA
  for (thresh in threshold_types) {
    if (grepl(thresh, file_name)) {
      threshold <- thresh
      break
    }
  }
  
  # Extract climate scenario information
  is_current <- grepl("projection.tif", file_name)
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
    threshold = threshold,
    scenario = scenario,
    time_period = time_period,
    is_current = is_current
  ))
}


#==============================================================================#
#                             4. Overlap metrics ----
#==============================================================================#

calculate_overlap_metrics <- function(calo_current, vertlept_current, 
                                      calo_future, vertlept_future,
                                      res_factor = 1) {
  # Function to calculate cell area in km²
  # (simplification - using mean resolution)
  cell_area_km2 <- mean(res(calo_current)) * mean(res(calo_current)) * res_factor
  
  # Current overlap
  current_overlap <- (calo_current > 0) & (vertlept_current > 0)
  current_overlap_cells <- sum(values(current_overlap), na.rm = TRUE)
  current_overlap_area <- current_overlap_cells * cell_area_km2
  
  # Future overlap
  future_overlap <- (calo_future > 0) & (vertlept_future > 0)
  future_overlap_cells <- sum(values(future_overlap), na.rm = TRUE)
  future_overlap_area <- future_overlap_cells * cell_area_km2
  
  # Species ranges
  calo_current_cells <- sum(values(calo_current > 0), na.rm = TRUE)
  vertlept_current_cells <- sum(values(vertlept_current > 0), na.rm = TRUE)
  calo_future_cells <- sum(values(calo_future > 0), na.rm = TRUE)
  vertlept_future_cells <- sum(values(vertlept_future > 0), na.rm = TRUE)
  
  # Areas in km²
  calo_current_area <- calo_current_cells * cell_area_km2
  vertlept_current_area <- vertlept_current_cells * cell_area_km2
  calo_future_area <- calo_future_cells * cell_area_km2
  vertlept_future_area <- vertlept_future_cells * cell_area_km2
  
  # Jaccard similarity index
  current_union <- (calo_current > 0) | (vertlept_current > 0)
  current_union_cells <- sum(values(current_union), na.rm = TRUE)
  current_jaccard <- current_overlap_cells / current_union_cells
  
  future_union <- (calo_future > 0) | (vertlept_future > 0)
  future_union_cells <- sum(values(future_union), na.rm = TRUE)
  future_jaccard <- future_overlap_cells / future_union_cells
  
  # Overlap percentage (of total range)
  current_overlap_percent <- 100 * current_overlap_cells / current_union_cells
  future_overlap_percent <- 100 * future_overlap_cells / future_union_cells
  
  # Overlap change
  overlap_change_percent <- ((future_overlap_percent - current_overlap_percent) / 
                               current_overlap_percent) * 100
  
  # Create metrics dataframe
  overlap <- data.frame(
    Metric = c(
      "Calophyllus Range (km²)",
      "Verticillium Range (km²)",
      "Overlap Area (km²)",
      "Overlap (% of total range)",
      "Jaccard Similarity Index"
    ),
    Current = c(
      round(calo_current_area, 1),
      round(vertlept_current_area, 1),
      round(current_overlap_area, 1),
      round(current_overlap_percent, 1),
      round(current_jaccard, 3)
    ),
    Future = c(
      round(calo_future_area, 1),
      round(vertlept_future_area, 1),
      round(future_overlap_area, 1),
      round(future_overlap_percent, 1),
      round(future_jaccard, 3)
    ),
    Change = c(
      paste0(round(((calo_future_area - calo_current_area) / calo_current_area) * 100, 1), "%"),
      paste0(round(((vertlept_future_area - vertlept_current_area) / vertlept_current_area) * 100, 1), "%"),
      paste0(round(((future_overlap_area - current_overlap_area) / current_overlap_area) * 100, 1), "%"),
      paste0(round(overlap_change_percent, 1), "%"),
      paste0(round(((future_jaccard - current_jaccard) / current_jaccard) * 100, 1), "%")
    )
  )
  return(overlap)
}


#==============================================================================#
#                       5. Metric calculation functions ----
#==============================================================================#


# Calculate area-based metrics
calculate_area_metrics <- function(current_binary, future_binary) {
  # Convert to terra objects if needed
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  
  # Ensure rasters have same geometry
  if (!terra::compareGeom(current_binary_terra, future_binary_terra, stopOnError = FALSE)) {
    future_binary_terra <- terra::resample(future_binary_terra, current_binary_terra)
  }
  
  # Calculate area of a single cell in square kilometers
  cell_area_m2 <- terra::res(current_binary_terra)[1] * terra::res(current_binary_terra)[2]
  cell_area_km2 <- cell_area_m2 / 1000000  # Convert to km²
  
  # Count cells for each category
  current_cells <- terra::global(current_binary_terra > 0, "sum", na.rm=TRUE)$sum
  future_cells <- terra::global(future_binary_terra > 0, "sum", na.rm=TRUE)$sum
  
  # Calculate overlap, expansion, and contraction
  overlap_mask <- (current_binary_terra > 0) & (future_binary_terra > 0)
  expansion_mask <- (future_binary_terra > 0) & (current_binary_terra == 0)
  contraction_mask <- (current_binary_terra > 0) & (future_binary_terra == 0)
  union_mask <- (current_binary_terra > 0) | (future_binary_terra > 0)
  
  overlap_cells <- terra::global(overlap_mask, "sum", na.rm=TRUE)$sum
  expansion_cells <- terra::global(expansion_mask, "sum", na.rm=TRUE)$sum
  contraction_cells <- terra::global(contraction_mask, "sum", na.rm=TRUE)$sum
  union_cells <- terra::global(union_mask, "sum", na.rm=TRUE)$sum
  
  # Calculate areas in square kilometers
  current_area <- current_cells * cell_area_km2
  future_area <- future_cells * cell_area_km2
  overlap_area <- overlap_cells * cell_area_km2
  expansion_area <- expansion_cells * cell_area_km2
  contraction_area <- contraction_cells * cell_area_km2
  union_area <- union_cells * cell_area_km2
  
  # Calculate derived metrics
  percent_change <- ifelse(current_area > 0, (future_area - current_area) / current_area * 100, 0)
  percent_stability <- ifelse(current_area > 0, overlap_area / current_area * 100, 0)
  percent_expansion <- ifelse(future_area > 0, expansion_area / future_area * 100, 0)
  percent_contraction <- ifelse(current_area > 0, contraction_area / current_area * 100, 0)
  spatial_disruption <- ifelse(union_area > 0, (current_area - future_area) / union_area, 0)
  
  # Similarity indices
  sorensen <- ifelse((current_area + future_area) > 0, 
                     2 * overlap_area / (current_area + future_area), 0)
  jaccard <- ifelse(union_area > 0, overlap_area / union_area, 0)
  
  return(list(
    "current_range_area" = current_area,
    "future_range_area" = future_area,
    "stable_presence_area" = overlap_area,
    "expansion_area" = expansion_area,
    "contraction_area" = contraction_area,
    "net_area_change" = future_area - current_area,
    "percent_change" = percent_change, 
    "percent_stability" = percent_stability,
    "percent_expansion" = percent_expansion,
    "percent_contraction" = percent_contraction,
    "spatial_disruption" = spatial_disruption,
    "sorensen_similarity" = sorensen,
    "jaccard_similarity" = jaccard,
    "range_stability" = percent_stability/100
  ))
}

# Create range change raster for visualization
calculate_range_change_raster <- function(current_binary, future_binary) {
  # Ensure rasters have same extent and resolution
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  
  if (!terra::compareGeom(current_binary_terra, future_binary_terra, stopOnError = FALSE)) {
    future_binary_terra <- terra::resample(future_binary_terra, current_binary_terra)
  }
  
  # Create a range change raster where:
  # 0 = Absence in both time periods
  # 1 = Expansion (absent in current, present in future)
  # 2 = Present in both time periods
  # 3 = Contraction (present in current, absent in future)
  range_change <- current_binary_terra + 2 * future_binary_terra
  range_change[range_change == 3] <- 2  # 1+2 = both present
  range_change[range_change == 1] <- 3  # Only present in current = contraction
  
  return(list("range_change_raster" = range_change))
}

# Calculate centroid shift and velocity
calculate_centroid_shift <- function(current_binary, future_binary, years_difference = NULL) 
{
    # Convert to terra objects if needed
    current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
    future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
    
    # Get coordinates for presence areas
    current_cells <- terra::as.data.frame(current_binary_terra, xy=TRUE)
    current_cells <- current_cells[current_cells[,3] > 0, c(1,2)]
    
    future_cells <- terra::as.data.frame(future_binary_terra, xy=TRUE)
    future_cells <- future_cells[future_cells[,3] > 0, c(1,2)]
    
    # Check if we have presence points
    if (nrow(current_cells) == 0 || nrow(future_cells) == 0) {
      return(list(
        "current_centroid" = c(NA, NA),
        "future_centroid" = c(NA, NA),
        "distance" = NA,
        "direction_degrees" = NA,
        "direction_cardinal" = NA,
        "velocity_km_yr" = NA,
        "velocity_km_decade" = NA
      ))
    }
    
    # Calculate centroids
    current_centroid <- c(mean(current_cells[,1]), mean(current_cells[,2]))
    future_centroid <- c(mean(future_cells[,1]), mean(future_cells[,2]))
    
    # Calculate distance between centroids (in km)
    distance_m <- sqrt(sum((future_centroid - current_centroid)^2))
    distance_km <- distance_m / 1000
    
    # Calculate direction
    angle <- atan2(future_centroid[2] - current_centroid[2], future_centroid[1] - current_centroid[1])
    direction <- (angle * 180 / pi) %% 360
    
    # Determine cardinal direction
    cardinal <- cut(direction, 
                    breaks = c(0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360),
                    labels = c("E", "NE", "N", "NW", "W", "SW", "S", "SE", "E"), 
                    include.lowest = TRUE)
    
    # Calculate velocity if years_difference is provided
    velocity_km_yr <- NA
    velocity_km_decade <- NA
    if(!is.null(years_difference) && !is.na(distance_km) && years_difference > 0) {
      velocity_km_yr <- distance_km / years_difference
      velocity_km_decade <- velocity_km_yr * 10
    }
    
    # Return list in the exact format needed
    return(list(
      "current_centroid" = current_centroid,
      "future_centroid" = future_centroid,
      "distance" = distance_km,
      "direction_degrees" = direction,
      "direction_cardinal" = as.character(cardinal),
      "velocity_km_yr" = velocity_km_yr,
      "velocity_km_decade" = velocity_km_decade
    ))
  }

# Calculate altitudinal shift
calculate_altitude_shift <- function(current_binary, future_binary, dem_raster) {
  # Convert to terra objects if needed
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  dem_raster_terra <- if(class(dem_raster)[1] == "RasterLayer") terra::rast(dem_raster) else dem_raster
  
  # Ensure matching projections
  if (terra::crs(current_binary_terra) != terra::crs(dem_raster_terra)) {
    dem_raster_terra <- terra::project(dem_raster_terra, current_binary_terra)
  }
  
  # Ensure matching geometries
  current_ext <- terra::ext(current_binary_terra)
  dem_ext <- terra::ext(dem_raster_terra)
  current_res <- terra::res(current_binary_terra)
  dem_res <- terra::res(dem_raster_terra)
  
  if (!all(current_ext == dem_ext) || !all(current_res == dem_res)) {
    dem_raster_terra <- terra::resample(dem_raster_terra, current_binary_terra)
  }
  
  # Extract elevation values at presence locations
  current_presence <- terra::as.data.frame(current_binary_terra, xy=TRUE)
  current_presence <- current_presence[current_presence[,3] > 0, c(1,2)]
  colnames(current_presence) <- c("x", "y")
  
  future_presence <- terra::as.data.frame(future_binary_terra, xy=TRUE)
  future_presence <- future_presence[future_presence[,3] > 0, c(1,2)]
  colnames(future_presence) <- c("x", "y")
  
  # Check for empty presence areas
  if (nrow(current_presence) == 0 || nrow(future_presence) == 0) {
    return(list(
      "current_mean_elevation" = NA,
      "future_mean_elevation" = NA,
      "current_min_elevation" = NA,
      "future_min_elevation" = NA,
      "current_max_elevation" = NA,
      "future_max_elevation" = NA,
      "mean_elevation_shift" = NA,
      "min_elevation_shift" = NA,
      "max_elevation_shift" = NA
    ))
  }
  
  # Extract elevation values
  current_elev <- terra::extract(dem_raster_terra, current_presence)
  future_elev <- terra::extract(dem_raster_terra, future_presence)
  
  # Check if extractions succeeded
  if (nrow(current_elev) == 0 || all(is.na(current_elev[,2])) ||
      nrow(future_elev) == 0 || all(is.na(future_elev[,2]))) {
    return(list(
      "current_mean_elevation" = NA,
      "future_mean_elevation" = NA,
      "current_min_elevation" = NA,
      "future_min_elevation" = NA,
      "current_max_elevation" = NA,
      "future_max_elevation" = NA,
      "mean_elevation_shift" = NA,
      "min_elevation_shift" = NA,
      "max_elevation_shift" = NA
    ))
  }
  
  # Calculate elevation statistics
  current_mean_elev <- mean(current_elev[,2], na.rm=TRUE)
  current_min_elev <- min(current_elev[,2], na.rm=TRUE)
  current_max_elev <- max(current_elev[,2], na.rm=TRUE)
  
  future_mean_elev <- mean(future_elev[,2], na.rm=TRUE)
  future_min_elev <- min(future_elev[,2], na.rm=TRUE)
  future_max_elev <- max(future_elev[,2], na.rm=TRUE)
  
  # Calculate shifts
  mean_shift <- future_mean_elev - current_mean_elev
  min_shift <- future_min_elev - current_min_elev
  max_shift <- future_max_elev - current_max_elev
  
  return(list(
    "current_mean_elevation" = current_mean_elev,
    "future_mean_elevation" = future_mean_elev,
    "current_min_elevation" = current_min_elev,
    "future_min_elevation" = future_min_elev,
    "current_max_elevation" = current_max_elev,
    "future_max_elevation" = future_max_elev,
    "mean_elevation_shift" = mean_shift,
    "min_elevation_shift" = min_shift,
    "max_elevation_shift" = max_shift
  ))
}

# Calculate Habitat Distance 
calculate_habitat_distance <- function(current_binary, future_binary) {
  # Convert to terra objects if needed
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  
  # Get coordinates of presence cells
  current_presence <- terra::as.data.frame(current_binary_terra, xy=TRUE)
  current_presence <- current_presence[current_presence[,3] > 0, c(1,2)]
  colnames(current_presence) <- c("x", "y")
  
  future_presence <- terra::as.data.frame(future_binary_terra, xy=TRUE)
  future_presence <- future_presence[future_presence[,3] > 0, c(1,2)]
  colnames(future_presence) <- c("x", "y")
  
  # Check if we have presence points in both distributions
  if (nrow(current_presence) == 0 || nrow(future_presence) == 0) {
    return(NA)
  }
  
  # Calculate centroid of current range
  current_centroid <- c(mean(current_presence$x), mean(current_presence$y))
  
  # Calculate distances from each future presence point to the current centroid
  distances <- numeric(nrow(future_presence))
  for (i in 1:nrow(future_presence)) {
    dx <- future_presence$x[i] - current_centroid[1]
    dy <- future_presence$y[i] - current_centroid[2]
    distances[i] <- sqrt(dx^2 + dy^2)
  }
  
  # Convert to kilometers (assuming projection in meters)
  distances_km <- distances / 1000
  
  # Calculate median distance
  median_distance <- median(distances_km)
  
  return(median_distance)
}

# Calculate habitat exposure (based on Choe et al. (2017))
calculate_habitat_exposure <- function(current_binary, future_binary) {
  # Validate input rasters
  if (!inherits(current_binary, "SpatRaster") || 
      !inherits(future_binary, "SpatRaster")) {
    stop("Inputs must be SpatRaster objects")
  }
  
  # Ensure rasters are identical in extent and resolution
  if (!compareGeom(current_binary, future_binary)) {
    stop("Current and future rasters must have identical extent and resolution")
  }
  
  # Calculate intersection of current and future ranges
  intersection <- current_binary * future_binary
  
  # Calculate total current range area
  # Using global() instead of cellStats() in terra
  current_range_area <- global(current_binary, fun = "sum", na.rm = TRUE)[[1]]
  
  # Calculate intersection area
  intersection_area <- global(intersection, fun = "sum", na.rm = TRUE)[[1]]
  
  # Calculate habitat exposure
  # Habitat Exposure = (Current Range Area - Intersection Area) / Current Range Area
  habitat_exposure <- (current_range_area - intersection_area) / current_range_area
  
  # Return the habitat exposure value
  return(habitat_exposure)
}


# Calculate IUCN Predicted Area of Occupancy (PAOO)
calculate_iucn_paoo <- function(current_binary, future_binary, grid_size = 2) {
  # Convert to terra objects if needed
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  
  # Convert grid_size from km to the units of the raster (usually meters)
  grid_size_m <- grid_size * 1000
  
  # Get the extent of the study area
  current_ext <- terra::ext(current_binary_terra)
  future_ext <- terra::ext(future_binary_terra)
  
  # Create a combined extent for both distributions
  combined_ext <- terra::ext(min(current_ext[1], future_ext[1]),
                             max(current_ext[2], future_ext[2]),
                             min(current_ext[3], future_ext[3]),
                             max(current_ext[4], future_ext[4]))
  
  # Create a standard IUCN grid with 2x2 km cells (or user-specified size)
  template <- terra::rast(extent = combined_ext, 
                          resolution = c(grid_size_m, grid_size_m), 
                          crs = terra::crs(current_binary_terra))
  
  # Resample current and future distribution to the standard grid
  current_resampled <- terra::resample(current_binary_terra, template, method = "max")
  future_resampled <- terra::resample(future_binary_terra, template, method = "max")
  
  # Convert to binary (any cell with value > 0 becomes 1)
  current_resampled <- current_resampled > 0
  future_resampled <- future_resampled > 0
  
  # Count occupied grid cells
  current_occupied_cells <- sum(terra::values(current_resampled), na.rm = TRUE)
  future_occupied_cells <- sum(terra::values(future_resampled), na.rm = TRUE)
  
  # Calculate AOO (each cell is grid_size × grid_size km²)
  current_aoo <- current_occupied_cells * (grid_size * grid_size)
  future_aoo <- future_occupied_cells * (grid_size * grid_size)
  
  # Calculate AOO change
  aoo_change <- future_aoo - current_aoo
  aoo_percent_change <- ifelse(current_aoo > 0, 
                               (aoo_change / current_aoo) * 100, 
                               NA)
  
  return(list(
    "current_aoo" = current_aoo,
    "future_aoo" = future_aoo,
    "aoo_change" = aoo_change,
    "aoo_percent_change" = aoo_percent_change
  ))
}

# Calculate range boundaries and their movements
calculate_range_boundaries <- function(current_binary, future_binary) {
  # Convert to terra objects if needed
  current_binary_terra <- if(class(current_binary)[1] == "RasterLayer") terra::rast(current_binary) else current_binary
  future_binary_terra <- if(class(future_binary)[1] == "RasterLayer") terra::rast(future_binary) else future_binary
  
  # Get coordinates of presence cells
  current_presence <- terra::as.data.frame(current_binary_terra, xy=TRUE)
  current_presence <- current_presence[current_presence[,3] > 0, c(1,2)]
  colnames(current_presence) <- c("x", "y")
  
  future_presence <- terra::as.data.frame(future_binary_terra, xy=TRUE)
  future_presence <- future_presence[future_presence[,3] > 0, c(1,2)]
  colnames(future_presence) <- c("x", "y")
  
  # Calculate current range boundaries
  if (nrow(current_presence) > 0) {
    current_n <- max(current_presence$y)
    current_s <- min(current_presence$y)
    current_e <- max(current_presence$x)
    current_w <- min(current_presence$x)
  } else {
    return(list(
      "n_boundary_shift" = 0,
      "s_boundary_shift" = 0,
      "e_boundary_shift" = 0,
      "w_boundary_shift" = 0,
      "leading_edge_dist" = NA,
      "leading_edge_direction" = NA,
      "trailing_edge_dist" = NA,
      "trailing_edge_direction" = NA
    ))
  }
  
  # Calculate future range boundaries
  if (nrow(future_presence) > 0) {
    future_n <- max(future_presence$y)
    future_s <- min(future_presence$y)
    future_e <- max(future_presence$x)
    future_w <- min(future_presence$x)
  } else {
    return(list(
      "n_boundary_shift" = 0,
      "s_boundary_shift" = 0,
      "e_boundary_shift" = 0,
      "w_boundary_shift" = 0,
      "leading_edge_dist" = NA,
      "leading_edge_direction" = NA,
      "trailing_edge_dist" = NA,
      "trailing_edge_direction" = NA
    ))
  }
  
  # Calculate boundary shifts in kilometers
  unit_conversion <- 1000  # Convert meters to kilometers
  
  # Boundary movements (positive values indicate expansion, negative values indicate contraction)
  n_shift <- (future_n - current_n) / unit_conversion  # North boundary shift
  s_shift <- (future_s - current_s) / unit_conversion  # South boundary shift
  e_shift <- (future_e - current_e) / unit_conversion  # East boundary shift
  w_shift <- (future_w - current_w) / unit_conversion  # West boundary shift
  
  # Identify leading edge (expansion front) - the boundary with the greatest positive shift
  shifts <- c(n_shift, s_shift, e_shift, w_shift)
  directions <- c("North", "South", "East", "West")
  
  # Find maximum expansion (leading edge)
  max_expansion <- max(shifts)
  if (max_expansion > 0) {
    leading_edge_idx <- which.max(shifts)
    leading_edge_dist <- shifts[leading_edge_idx]
    leading_edge_direction <- directions[leading_edge_idx]
  } else {
    leading_edge_dist <- NA
    leading_edge_direction <- NA
  }
  
  # Find maximum contraction (trailing edge) - the boundary with the greatest negative shift
  min_shift <- min(shifts)
  if (min_shift < 0) {
    trailing_edge_idx <- which.min(shifts)
    trailing_edge_dist <- abs(shifts[trailing_edge_idx])  # Convert to positive value for distance
    trailing_edge_direction <- directions[trailing_edge_idx]
  } else {
    trailing_edge_dist <- NA
    trailing_edge_direction <- NA
  }
  
  return(list(
    "n_boundary_shift" = n_shift,
    "s_boundary_shift" = s_shift,
    "e_boundary_shift" = e_shift,
    "w_boundary_shift" = w_shift,
    "leading_edge_dist" = leading_edge_dist,
    "leading_edge_direction" = leading_edge_direction,
    "trailing_edge_dist" = trailing_edge_dist,
    "trailing_edge_direction" = trailing_edge_direction
  ))
}

# Calculate rate of change (km/decade)
calculate_rate_of_boundary_change <- function(boundary_shifts, years_difference) {
  # Convert years to decades
  decades <- years_difference / 10
  
  # Calculate rates for each boundary
  n_rate <- ifelse(is.na(boundary_shifts$n_boundary_shift) || decades == 0, 
                   0, boundary_shifts$n_boundary_shift / decades)
  
  s_rate <- ifelse(is.na(boundary_shifts$s_boundary_shift) || decades == 0, 
                   0, boundary_shifts$s_boundary_shift / decades)
  
  e_rate <- ifelse(is.na(boundary_shifts$e_boundary_shift) || decades == 0, 
                   0, boundary_shifts$e_boundary_shift / decades)
  
  w_rate <- ifelse(is.na(boundary_shifts$w_boundary_shift) || decades == 0, 
                   0, boundary_shifts$w_boundary_shift / decades)
  
  # Leading and trailing edge rates
  leading_edge_rate <- NA
  trailing_edge_rate <- NA
  
  if (!is.na(boundary_shifts$leading_edge_dist) && decades > 0) {
    leading_edge_rate <- boundary_shifts$leading_edge_dist / decades
  }
  
  if (!is.na(boundary_shifts$trailing_edge_dist) && decades > 0) {
    trailing_edge_rate <- boundary_shifts$trailing_edge_dist / decades
  }
  
  return(list(
    "n_rate" = n_rate,
    "s_rate" = s_rate,
    "e_rate" = e_rate,
    "w_rate" = w_rate,
    "leading_edge_rate" = leading_edge_rate,
    "trailing_edge_direction" = boundary_shifts$trailing_edge_direction,
    "leading_edge_direction" = boundary_shifts$leading_edge_direction,
    "trailing_edge_rate" = trailing_edge_rate
  ))
}


#==============================================================================#
#                         6. Combined analysis function ----
#==============================================================================#


analyze_range_shifts <- function(species, current_raster, future_raster, 
                                 threshold, scenario_info, dem_raster = NULL) {
  
  # Get full species name 
  species_fullname <- if(species %in% names(species_fullnames)) species_fullnames[species] else species
  
  cat("\nAnalyzing range shifts for", species_fullname, "\n")
  cat("Threshold:", threshold, "\n")
  cat("Scenario:", scenario_info$scenario, "\n")
  cat("Time period:", scenario_info$time_period, "\n")
  
  # Calculate basic range change raster
  range_shifts <- calculate_range_change_raster(current_raster, future_raster)
  
  # Calculate centroid shift
  centroid_shift <- calculate_centroid_shift(current_raster, future_raster)
  
  # Calculate altitude shift if DEM is provided
  alt_shift <- if (!is.null(dem_raster)) {
    calculate_altitude_shift(current_raster, future_raster, dem_raster)
  } else NULL
  
  habitat_exposure <- calculate_habitat_exposure(current_raster, future_raster)
  
  # Calculate range boundary movements
  boundary_shifts <- calculate_range_boundaries(current_raster, future_raster)
  
  # Calculate rates of change
  base_year <- 2010  # Assuming current predictions are for ~2010
  
  # Extract end year from time period (format: "2011-2040", "2041-2070", etc.)
  end_year <- NA
  if (scenario_info$time_period != "current") {
    time_parts <- strsplit(scenario_info$time_period, "-")[[1]]
    if (length(time_parts) == 2) {
      end_year <- as.numeric(time_parts[2])
    }
  }
  
  # Default to 30 years if can't parse
  years_difference <- ifelse(!is.na(end_year), end_year - base_year, 30)
  
  # Calculate rates of change
  change_rates <- calculate_rate_of_boundary_change(boundary_shifts, years_difference)
  
  # Calculate area metrics
  area_metrics <- calculate_area_metrics(current_raster, future_raster)
  
  # Calculate habitat distance
  habitat_distance <- calculate_habitat_distance(current_raster, future_raster)
  
  # Calculate IUCN PAOO
  paoo <- calculate_iucn_paoo(current_raster, future_raster)
  
  # Combine all metrics into a data frame
  metrics <- data.frame(
    # Metadata
    species = species,
    species_fullname = species_fullname,
    threshold = threshold,
    scenario = scenario_info$scenario,
    time_period = scenario_info$time_period,
    
    # Area metrics
    current_range_area = area_metrics$current_range_area,
    future_range_area = area_metrics$future_range_area,
    stable_presence_area = area_metrics$stable_presence_area,
    expansion_area = area_metrics$expansion_area,
    contraction_area = area_metrics$contraction_area,
    net_area_change = area_metrics$net_area_change,
    percent_change = area_metrics$percent_change,
    percent_stability = area_metrics$percent_stability,
    percent_expansion = area_metrics$percent_expansion,
    percent_contraction = area_metrics$percent_contraction,
    
    # Centroid metrics
    centroid_shift_distance = centroid_shift$distance,
    centroid_shift_direction_degrees = centroid_shift$direction_degrees,
    centroid_shift_direction = centroid_shift$direction_cardinal,
    centroid_shift_velocity_km_yr = centroid_shift$velocity_km_yr,
    centroid_shift_velocity_km_dec = centroid_shift$velocity_km_decade,

    # Habitat distance metric
    habitat_distance_km = habitat_distance,
    
    # Habitat exposure metric
    habitat_exposure_km2 = habitat_exposure,
    
    # Spatial disruption
    spatial_disruption = area_metrics$spatial_disruption,
    
    # IUCN PAOO metrics
    current_aoo_km2 = paoo$current_aoo,
    future_aoo_km2 = paoo$future_aoo,
    aoo_change = paoo$aoo_change,
    aoo_percent_change = paoo$aoo_percent_change,
    
    # Boundary shifts
    north_boundary_shift = boundary_shifts$n_boundary_shift,
    south_boundary_shift = boundary_shifts$s_boundary_shift,
    east_boundary_shift = boundary_shifts$e_boundary_shift,
    west_boundary_shift = boundary_shifts$w_boundary_shift,
    leading_edge_distance = boundary_shifts$leading_edge_dist,
    leading_edge_direction = boundary_shifts$leading_edge_direction,
    trailing_edge_distance = boundary_shifts$trailing_edge_dist,
    trailing_edge_direction = boundary_shifts$trailing_edge_direction,
    
    # Rate metrics
    years_projected = years_difference,
    decades_projected = years_difference / 10,
    north_shift_rate = change_rates$n_rate,
    south_shift_rate = change_rates$s_rate,
    east_shift_rate = change_rates$e_rate,
    west_shift_rate = change_rates$w_rate,
    leading_edge_rate = change_rates$leading_edge_rate,
    trailing_edge_rate = change_rates$trailing_edge_rate,
    
    # Similarity metrics
    sorensen_similarity = area_metrics$sorensen_similarity,
    jaccard_similarity = area_metrics$jaccard_similarity,
    range_stability = area_metrics$range_stability,
    stringsAsFactors = FALSE
  )
  
  # Add altitude metrics if available
  if (!is.null(alt_shift)) {
    altitude_data <- data.frame(
      current_mean_elevation = alt_shift$current_mean_elevation,
      future_mean_elevation = alt_shift$future_mean_elevation,
      current_min_elevation = alt_shift$current_min_elevation,
      future_min_elevation = alt_shift$future_min_elevation,
      current_max_elevation = alt_shift$current_max_elevation,
      future_max_elevation = alt_shift$future_max_elevation,
      mean_elevation_shift = alt_shift$mean_elevation_shift,
      min_elevation_shift = alt_shift$min_elevation_shift,
      max_elevation_shift = alt_shift$max_elevation_shift,
      stringsAsFactors = FALSE
    )
    metrics <- cbind(metrics, altitude_data)
    
  }
  
  # Generate output file names
  metrics_dir <- file.path(data06_out, "output","metrics", "extra_settings", "summary")
  
  # Create file name for metrics that includes threshold type
  metrics_file <- file.path(metrics_dir, 
                            paste0(species, "_", threshold, "_", 
                                   scenario_info$scenario, "_", 
                                   scenario_info$time_period, "_metrics.csv"))
  
  # Create file name for metrics that includes threshold type
  overlap_file <- file.path(metrics_dir, 
                            paste0(species, "_", threshold, "_", 
                                   scenario_info$scenario, "_", 
                                   scenario_info$time_period, "_overlap.csv"))
  
  # Save metrics to CSV
  write.csv(metrics, metrics_file, row.names = FALSE)
  cat("Metrics saved to:", metrics_file, "\n")
  
  return(metrics)
}


#==============================================================================#
#                        7. Main function to process all  ----    
#==============================================================================#


process_all_combinations <- function() {
  # Initialize empty data frame to store all metrics
  all_metrics <- data.frame()
  
  # Get lists of current and future binary rasters
  current_files <- list.files(binary_cur_dir, pattern=".tif", full.names = TRUE)
  future_files <- list.files(binary_fut_dir, pattern=".tif", full.names = TRUE)
  
  # Organise current projections by species and threshold
  current_rasters <- list()
  for (cur_file in current_files) {
    info <- extract_scenario_info(cur_file)
    species <- info$species
    threshold <- info$threshold
    
    if (!is.na(species) && !is.na(threshold) && info$is_current) {
      key <- paste(species, threshold, sep="_")
      current_rasters[[key]] <- terra::rast(cur_file)
      cat("Loaded current projection for", species, "with threshold", threshold, "\n")
    }
  }
  
  # Process each future prediction
  for (future_file in future_files) {
    # Extract scenario information
    scenario_info <- extract_scenario_info(future_file)
    
    # Skip if scenario info is incomplete
    if (is.na(scenario_info$species) || 
        is.na(scenario_info$scenario) || 
        is.na(scenario_info$time_period) ||
        is.na(scenario_info$threshold) ||
        scenario_info$is_current) {
      cat("Skipping", basename(future_file), "- incomplete info or is current projection\n")
      next
    }
    
    # Construct key to find matching current raster
    species <- scenario_info$species
    threshold <- scenario_info$threshold
    key <- paste(species, threshold, sep="_")
    
    # Skip if no current raster for this species-threshold combination
    if (!key %in% names(current_rasters)) {
      cat("Skipping", basename(future_file), "- no matching current projection available\n")
      next
    }
    
    # Get current raster for this species-threshold combination
    current_raster <- current_rasters[[key]]
    
    # Load future prediction raster
    future_raster <- terra::rast(future_file)
    
    # Calculate metrics
    metrics <- analyze_range_shifts(
      species = species,
      current_raster = current_raster,
      future_raster = future_raster,
      threshold = threshold,
      scenario_info = scenario_info,
      dem_raster = elev_MDG
    )
    
    # Append to all metrics
    all_metrics <- rbind(all_metrics, metrics)
  }
  
  # Save combined metrics
  if (nrow(all_metrics) > 0) {
    all_metrics_file <- file.path(data06_out, "output","metrics", "extra_settings", "summary", "all_species_all_thresholds_metrics.csv")
    write.csv(all_metrics, all_metrics_file, row.names = FALSE)
    cat("\nAll metrics saved to:", all_metrics_file, "\n")
    
    # Create summary metrics by species, threshold and scenario
    summary_metrics <- aggregate(
      cbind(percent_change, percent_stability, percent_expansion, percent_contraction, 
            centroid_shift_distance, habitat_distance_km, habitat_exposure_km2, spatial_disruption) ~ 
        species + threshold + scenario, 
      data = all_metrics, 
      FUN = mean, 
      na.rm = TRUE
    )
    
    # Save summary metrics
    summary_file <- file.path(data06_out, "output", "metrics", "extra_settings","summary", "summary_by_scenario_threshold.csv")
    write.csv(summary_metrics, summary_file, row.names = FALSE)
    cat("Summary metrics saved to:", summary_file, "\n")
  } else {
    cat("No metrics were calculated. Check input files.\n")
  }
  
  return(all_metrics)
  
  calculate_species_overlaps <- function() {
    cat("\nCalculating species overlap metrics...\n")
    
    # Get lists of current and future binary rasters
    current_files <- list.files(binary_cur_dir, pattern=".tif", full.names = TRUE)
    future_files <- list.files(binary_fut_dir, pattern=".tif", full.names = TRUE)
    
    # Organiserasters by species, threshold, scenario, and time period
    all_rasters <- list()
    
    # Process current files
    for (file in current_files) {
      info <- extract_scenario_info(file)
      if (!is.na(info$species) && !is.na(info$threshold)) {
        key <- paste(info$species, info$threshold, "current", "current", sep="_")
        all_rasters[[key]] <- terra::rast(file)
      }
    }
    
    # Process future files
    for (file in future_files) {
      info <- extract_scenario_info(file)
      if (!is.na(info$species) && !is.na(info$threshold) && 
          !is.na(info$scenario) && !is.na(info$time_period)) {
        key <- paste(info$species, info$threshold, info$scenario, info$time_period, sep="_")
        all_rasters[[key]] <- terra::rast(file)
      }
    }
    
    # Define the combinations to calculate
    # You can define different threshold combinations here
    threshold_combinations <- list(
      list(calo="opt2", vertlept="opt4"),
      list(calo="opt12", vertlept="opt12")
      # Add more combinations as needed
    )
    
    # Process each scenario and time period
    scenarios <- unique(sapply(strsplit(names(all_rasters), "_"), function(x) x[3]))
    scenarios <- scenarios[scenarios != "current"]
    
    time_periods <- unique(sapply(strsplit(names(all_rasters), "_"), function(x) x[4]))
    time_periods <- time_periods[time_periods != "current"]
    
    all_overlaps <- data.frame()
    
    for (scenario in scenarios) {
      for (time_period in time_periods) {
        for (thresh_combo in threshold_combinations) {
          calo_thresh <- thresh_combo$calo
          vertlept_thresh <- thresh_combo$vertlept
          
          # Get current rasters
          calo_current_key <- paste("calo", calo_thresh, "current", "current", sep="_")
          vertlept_current_key <- paste("vertlept", vertlept_thresh, "current", "current", sep="_")
          
          # Get future rasters
          calo_future_key <- paste("calo", calo_thresh, scenario, time_period, sep="_")
          vertlept_future_key <- paste("vertlept", vertlept_thresh, scenario, time_period, sep="_")
          
          # Check if all required rasters exist
          if (all(c(calo_current_key, vertlept_current_key, calo_future_key, vertlept_future_key) %in% names(all_rasters))) {
            cat("Calculating overlap for:", calo_thresh, vertlept_thresh, scenario, time_period, "\n")
            
            # Get the rasters
            calo_current <- all_rasters[[calo_current_key]]
            vertlept_current <- all_rasters[[vertlept_current_key]]
            calo_future <- all_rasters[[calo_future_key]]
            vertlept_future <- all_rasters[[vertlept_future_key]]
            
            # Calculate overlap
            overlap <- calculate_overlap_metrics(
              calo_current, vertlept_current,
              calo_future, vertlept_future
            )
            
            # Add metadata
            overlap$scenario <- scenario
            overlap$time_period <- time_period
            overlap$calo_threshold <- calo_thresh
            overlap$vertlept_threshold <- vertlept_thresh
            
            # Save this specific overlap combination
            overlap_file <- file.path(
              data06_out, "output", "metrics", "extra_settings", "summary",
              paste0("overlap_", calo_thresh, "_", vertlept_thresh, "_", 
                     scenario, "_", time_period, ".csv")
            )
            write.csv(overlap, overlap_file, row.names = FALSE)
            
            # Append to all overlaps
            all_overlaps <- rbind(all_overlaps, overlap)
          } else {
            cat("Skipping overlap for:", calo_thresh, vertlept_thresh, scenario, time_period, 
                "- missing one or more required rasters\n")
          }
        }
      }
    }
    
    # Save all overlaps to a single file
    if (nrow(all_overlaps) > 0) {
      all_overlaps_file <- file.path(
        data06_out, "output", "metrics", "extra_settings", "summary",
        "all_species_overlap_metrics.csv"
      )
      write.csv(all_overlaps, all_overlaps_file, row.names = FALSE)
      cat("All overlap metrics saved to:", all_overlaps_file, "\n")
    }
  }
}

# Run the analysis
cat("\n=== Starting Range Shift Metrics Analysis ===\n")
all_metrics <- process_all_combinations()
cat("\n=== Range Shift Metrics Analysis Completed ===\n")

# # Clean up
gc()
rm(list = ls())


#==============================================================================#
#                           ----  End of workflow ----
#==============================================================================#