# title: "c01__extract_species_enviro_data.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-01-12"
# update: "2025-12-31"
# aim: processing and extracting environmental variables for each species occurrence, 
#      year-range, and climate pathway
# N.B. the 4. Extraction section is using outputs from:
  # 00__species_occurrence_cleaning\output\combined_occs 
# which are pre-processed presence-background data (currently with more than one background point method)
# This code has been adapted for use within the xxx Linux cluster (University of xxx)

#==============================================================================#
#                              0. Workspace set up ----
#==============================================================================#

# individual packages that require installing directly in R interactive session while on cluster
library(terra)
library(doParallel)
library(foreach)
library(stringr)
library(gdalUtilities)
library(future.apply)
library(data.table)


# File paths on cluster and mnt from transfer folder
mnt_folder <- "/mnt/xxx/01__extract_species_enviro_data" # T drive mount on cluster where data is stored
import_dir <- "/import/xxx/01__extract_species_enviro_data" # cluster folder where AOI and template is stored (updated for ecoc9)
input_dir <- file.path(import_dir, "input")
terra_tmp <- "/import/01__extract_species_enviro_data/tmp"

template_path <- file.path(import_dir, "input", "AOI", "AFR_MDG_aea.tif") # this has been updated land mass = 1, rest is NA
enviro_base_path <- file.path(mnt_folder, "intermediate", "processed_environmental") # current or future > year range > bioclim > scenario(sspxyz)
#forest_path <- file.path(mnt_folder, "intermediate", "processed_forest") # 4 processed forest files (2020, 2040, 2070, 2100)
aspect_file <- file.path(mnt_folder, "intermediate", "processed_aspect", "aspect_processed.tif") # a single processed aspect file

# Create output directories
import <- "/import/xxx/01__extract_species_enviro_data/output"
output_dir <- file.path(import, "16_02_2025_extract")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "stacked"), showWarnings = FALSE) # fully stacked rasters of all env data
dir.create(file.path(output_dir, "extracted"), showWarnings = FALSE) # csv table of fully stacked/extracted for each spp

# limit memory
memory.limit(size = 128 * 1024)  # Set memory limit to 128GB

# # Set up logging to help debug
# log_file <- file.path(output_dir, "processing_log.txt")
# cat("Processing log\n\n", file = log_file)


# Define constants
scenarios <- c("current","ssp126", "ssp370", "ssp585")
year_ranges <- list(
  current = "1981-2010",
  future = c("2011-2040", "2041-2070","2071-2100")
)

species <- c("Calo_occs_combined_ecor_619", "Calo_occs_combined_whole_619", 
             "Vert_Lept_occs_2far_374", "Vert_Lept_occs_whole_374")

# CRS: African Equal Area
target_crs <- "ESRI:102022"


#==============================================================================#
#                         ---- 0. Helper Functions ----
#==============================================================================#

# Function to log messages
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste(timestamp, "-", message)
  # cat(full_message, "\n", file = log_file, append = TRUE)
  cat(full_message, "\n")  # Also print to console
}

# Function to check if a file exists
file_exists <- function(file_path) {
  return(file.exists(file_path))
}

# Function to load species data before stacking
load_species_data <- function(species) {
  file_path <- file.path(import_dir, "input", "species", paste0(species, ".RData"))
  log_message(paste("Loading data for", species, "from", file_path))
  
  if (file_exists(file_path)) {
    load(file_path)
    if (exists(species)) {
      data <- get(species)
      log_message(paste("Loaded occurrence data for", species, ". Dimensions:", paste(dim(data), collapse = "x")))
      return(data)
    } else {
      log_message(paste("Error: Data for", species, "not found in the loaded file"))
      return(NULL)
    }
  } else {
    log_message(paste("Error: File not found for", species))
    return(NULL)
  }
}

# Load template raster
template <- terra::rast(template_path)
cat("Template raster loaded successfully.\n")



#==============================================================================#
#                    ---- 1. Climate data processing ----
#==============================================================================#

# process_climate_file <- function(file, template) {
#   tryCatch({
#     log_message(paste("Processing climate file:", file))
# 
#     # Create output file name
#     file_name <- basename(file)
#     new_file_name <- paste0(tools::file_path_sans_ext(file_name), "_v2.tif")
#     output_path <- file.path(output_dir, "processed_climate", new_file_name)
# 
#     # Check if output file already exists
#     if (file_exists(output_path)) {
#       log_message(paste("Processed file already exists:", output_path))
#       return(output_path)
#     }
# 
#     # Load the climate raster
#     climate_raster <- terra::rast(file)
# 
#     # Ensure the climate raster has the same CRS as the template
#     if (terra::crs(climate_raster) != terra::crs(template)) {
#       climate_raster <- terra::project(climate_raster, terra::crs(template))
#     }
# 
#     # Resample to match the template resolution and extent
#     climate_resampled <- terra::resample(climate_raster, template)
# 
#     gc()
#     rm(climate_raster)
# 
#     # Mask using the template
#     climate_masked <- terra::mask(climate_resampled, template)
# 
#     gc()
#     rm(climate_resampled)
# 
#     # Write the processed raster
#     terra::writeRaster(climate_masked, output_path, overwrite = TRUE)
# 
#     log_message(paste("Climate file processed and saved:", output_path))
#     return(output_path)
#   }, error = function(e) {
#     log_message(paste("Error processing climate file:", file, "-", e$message))
#     return(NULL)
#   })
# 
#   gc()
# }
# 
# process_climate_data <- function(template, scenario, year_range) {
#   log_message(paste("Starting climate data processing for", scenario, year_range))
# 
#   # Determine the correct path based on scenario
#   if (scenario == "current") {
#     climate_path <- file.path(enviro_base_path,"current", year_range, "bioclim")
#   } else {
#     climate_path <- file.path(enviro_base_path, "future", year_range, "bioclim", scenario)
#   }
# 
#   # List climate files
#   climate_files <- list.files(enviro_base_path, pattern = "\\.tif$", full.names = TRUE)
# 
#   if (length(climate_files) == 0) {
#     log_message(paste("No climate files found for", year_range, scenario))
#     return(NULL)
#   }
# 
#   # Process files
#   processed_files <- lapply(climate_files, process_climate_file, template = template)
# 
#   # Remove NULL values (failed processing)
#   processed_files <- processed_files[!sapply(processed_files, is.null)]
# 
#   log_message(paste("Completed climate data processing for", year_range, scenario, ". Processed", length(processed_files), "files"))
# 
#   return(processed_files)
# }

#==============================================================================#
#                    ---- 2. Forest data processing ----
#==============================================================================#

# process_forest_file <- function(file, template, year_range) {
#   tryCatch({
#     log_message(paste("Processing forest file:", file))
# 
#     # Create output file name
#     file_name <- basename(file)
#     new_file_name <- paste0("processed_", year_range, ".tif")
#     output_path <- file.path(output_dir, "processed_forest", new_file_name)
# 
#     # Check if output file already exists
#     if (file.exists(output_path)) {
#       log_message(paste("Processed file already exists:", output_path))
#       return(output_path)
#     }
# 
#     r <- terra::rast(file)
# 
#     # Check if value 3 is present in the raster (there is only 1 of the 4 files that has this)
#     if (3 %in% unique(r)) {
# 
#       # More efficient processing
#       process_raster <- function(r) {
#         values <- terra::values(r)
# 
#         values[values %in% c(1, 2)] <- 0
#         values[values == 3] <- 1
#         values[is.na(values)] <- 0
# 
#         terra::values(r) <- values
#         return(r)
#       }
# 
#       r_processed <- process_raster(r)
#       gc()
#       rm(r)
# 
#       # Calculate aggregation factor
#       input_res <- terra::res(r_processed)[1]
#       target_res <- 1000
#       agg_factor <- ceiling(target_res / input_res)
# 
#       cat("Aggregating from", input_res, "m to", target_res, "m resolution.\n")
#       cat("Aggregation factor:", agg_factor, "\n")
# 
#       # Aggregate the raster using the above calculated factor using mean function
#       r_aggregated <- terra::aggregate(r_processed, fact=agg_factor, fun="mean")
#       gc()
#       rm(r_processed)
# 
#       # Resample to match template
#       r_resampled <- terra::resample(r_aggregated, template, method="bilinear")
#       gc()
#       rm(r_aggregated)
# 
#       # Final mask to template
#       r_final <- terra::mask(r_resampled, template)
#       gc()
#       rm(r_resampled)
# 
#       # Write processed raster
#       terra::writeRaster(r_final, output_path, overwrite = TRUE, filetype = "GTiff")
# 
#       cat("Raster processed and written to:", output_path, "\n")
#       return(list(status = "success", path = output_path, file = file))
#     } else {
#       # For other rasters, assume 1 is forest, everything else is non-forest
#       cat("Value 3 not present in raster. Skipping processing.\n")
#       return(list(status = "skipped", path = file, message = "Value 3 not present", file = file))
#     }
# 
#     log_message(paste("Forest file processed and saved:", output_path))
#     return(output_path)
#   }, error = function(e) {
#     log_message(paste("Error processing forest file:", file, "-", e$message))
#     return(NULL)
#   })
# }
# 
# process_forest_data <- function(template) {
#   log_message("Starting forest data processing")
# 
#   # List forest files
#   forest_files <- list.files(forest_path, pattern = "^fcc_.*_AFR_aea\\.tif$", full.names = TRUE)
# 
#   processed_files <- list()
# 
#   for (file in forest_files) {
#     # Extract year range from file name
#     year_range <- gsub(".*fcc_(\\d{4}-\\d{4})_AFR_aea\\.tif", "\\1", basename(file))
# 
#     # Process the file
#     processed_file <- process_forest_file(file, template, year_range)
# 
#     if (!is.null(processed_file)) {
#       processed_files[[year_range]] <- processed_file
#     }
#   }
# 
#   log_message(paste("Completed forest data processing. Processed", length(processed_files), "files"))
# 
#   return(processed_files)
# }
# 

#==============================================================================#
#                           ---- 3. Data stacking ----
#==============================================================================#


# new stacking function updated for environmental data already processed (12.01.2025)

tryCatch({
stack_environmental_data <- function(scenario, year_range, template) {
  log_message(paste("Stacking data for", scenario, year_range))
  
  # Determine paths based on scenario
  if (scenario == "current") {
    climate_path <- file.path(enviro_base_path, "current", year_range, "bioclim")
    forest_path <- file.path(enviro_base_path, "current", year_range, "forest")
  } else {
    climate_path <- file.path(enviro_base_path, "future", year_range, "bioclim", scenario)
    forest_path <- file.path(enviro_base_path, "future", year_range, "forest")
  }
  
  # List files
  climate_files <- list.files(climate_path, pattern = "\\.tif$", full.names = TRUE)
  forest_file <- list.files(forest_path, pattern = "^processed_.*\\.tif$", full.names = TRUE)[1]
  
  all_files <- c(climate_files, forest_file, aspect_file)
  
  if (length(all_files) == 0) {
    log_message(paste("No files found for", scenario, year_range))
    return(NULL)
  }
  
  # Create an empty SpatRaster with the same properties as the template
  stacked <- terra::rast(template)
  
  # Add each layer to the stack
  for (file in all_files) {
    layer <- terra::rast(file)
    # if (terra::crs(layer) != terra::crs(template)) {
    #   layer <- terra::project(layer, terra::crs(template))
    # }
    # layer <- terra::resample(layer, template)
    layer <- terra::mask(layer, template)
    stacked <- c(stacked, layer)
  }
  
  # Remove the first layer (which was empty)
  stacked <- stacked[[2:terra::nlyr(stacked)]]
  
  # Save stacked raster
  output_file <- file.path(output_dir, "stacked", paste0(scenario, "_", year_range, ".tif"))
  terra::writeRaster(stacked, output_file, overwrite = TRUE)
  
  log_message(paste("Stacked raster saved:", output_file))
  return(output_file)
}
})


#==============================================================================#
#                          ---- 4. Data extraction ----
#==============================================================================#


extract_environmental_data <- function(species_data, species_name, stacked_file, scenario, year_range) {
  log_message(paste("Extracting data for", species_name, "from", stacked_file))
  
  # Check if output file already exists
  output_file <- file.path(output_dir, "extracted", paste0(species_name, "_", scenario, "_", year_range, ".csv"))
  if (file_exists(output_file)) {
    log_message(paste("Extracted file already exists:", output_file))
    return(output_file)
  }
  
  # Load stacked raster
  stacked_env <- terra::rast(stacked_file)
  
  # Extract values
  extracted <- terra::extract(stacked_env, species_data[, c("x", "y")])
  
  # Combine occurrence data with extracted values
  result <- cbind(species_data, extracted)
  
  # Add scenario and year_range columns
  result$scenario <- scenario
  result$year_range <- year_range
  
  # Save results
  fwrite(result, output_file)
  
  log_message(paste("Extraction completed. Results saved to", output_file))
  return(output_file)
}


#==============================================================================#
#                          5. Main execution function ----
#==============================================================================#


main <- function() {
  log_message("Starting main workflow")
  
  # Create a list of common objects to export to all workers
  common_objects <- list(
    template = template,
    scenarios = scenarios,
    year_ranges = year_ranges,
    #process_climate_data = process_climate_data,
    #process_climate_data_env = environment(process_climate_data),
    #process_forest_data = process_forest_data,
    #process_forest_data_env = environment(process_forest_data),
    stack_environmental_data = stack_environmental_data,
    load_species_data = load_species_data,
    extract_environmental_data = extract_environmental_data,
    species = species
  )

  # # 1. Process climate data for each scenario and year range
  # log_message("Processing climate data")
  # climate_data <- foreach(scenario = scenarios, .packages = c("terra"),
  #                         .export = c(names(common_objects),
  #                                     ls(envir = common_objects$process_climate_data_env),
  #                                     ls(envir = common_objects$process_forest_data_env))) %dopar% {
  #   if (scenario == "current") {
  #     list(current = process_climate_data(template, scenario, year_ranges$current))
  #   } else {
  #     lapply(year_ranges$future, function(year_range) {
  #       process_climate_data(template, scenario, year_range)
  #     })
  #   }
  # }
  # names(climate_data) <- scenarios
  # 
  # 
  # # 2. Process forest data (only once, shared across all scenarios)
  # log_message("Processing forest data")
  # processed_forest <- process_forest_data(template)
  # 
  
  # Initialize lists to store results
  all_results <- list()
  
  # 3. Stack environmental data
  log_message("Stacking environmental data")
  stacked_files <- list()
  
  # Process each scenario
  for(scenario in scenarios) {
    if(scenario == "current") {
      log_message(paste("Processing current scenario:", scenario))
      stacked_files[[scenario]] <- list(
        current = stack_environmental_data(scenario, year_ranges$current, template)
      )
    } else {
      log_message(paste("Processing future scenario:", scenario))
      stacked_files[[scenario]] <- list()
      for(yr in year_ranges$future) {
        stacked_files[[scenario]][[yr]] <- stack_environmental_data(scenario, yr, template)
      }
    }
  }
  
  # Store stacked files in results
  all_results$stacked_files <- stacked_files
  
  # Print stacked files paths for debugging
  log_message("Stacked files created:")
  print(str(stacked_files))  # Add this to see the structure
  
  # 4. Extract data for species
  log_message("Extracting data for species")
  extracted_files <- list()
  
  # Loop through each species with explicit error handling
  for(sp in species) {
    tryCatch({
      log_message(paste("Processing species:", sp))
      species_data <- load_species_data(sp)
      
      if (!is.null(species_data)) {
        sp_extracted <- list()
        
        # Extract for current scenario with path verification
        current_path <- stacked_files[["current"]][["current"]]
        log_message(paste("Current scenario path:", current_path))
        
        if(!is.null(current_path) && file.exists(current_path)) {
          current_key <- paste(sp, "current", year_ranges$current, sep = "_")
          sp_extracted[[current_key]] <- extract_environmental_data(
            species_data, sp, current_path, "current", year_ranges$current
          )
        } else {
          log_message("ERROR: Current scenario stacked file not found")
        }
        
        # Extract for future scenarios with path verification
        for(scenario in scenarios[-1]) {
          for(year_range in year_ranges$future) {
            future_path <- stacked_files[[scenario]][[year_range]]
            log_message(paste("Future scenario path:", future_path))
            
            if(!is.null(future_path) && file.exists(future_path)) {
              future_key <- paste(sp, scenario, year_range, sep = "_")
              sp_extracted[[future_key]] <- extract_environmental_data(
                species_data, sp, future_path, scenario, year_range
              )
            } else {
              log_message(paste("ERROR: Future scenario stacked file not found for", scenario, year_range))
            }
          }
        }
        
        extracted_files <- c(extracted_files, sp_extracted)
      }
    }, error = function(e) {
      log_message(paste("Error processing species", sp, ":", e$message))
    })
  }
  
  # Store extracted files in results
  all_results$extracted_files <- unlist(extracted_files, recursive = FALSE)
  
  return(all_results)
}



#==============================================================================#
#                          ---- Run workflow ----
#==============================================================================#


# Set up parallel processing

# Memory management
#options(future.globals.maxSize = 8000 * 1024^2) # 8GB
terraOptions(memfrac=0.5, tempdir = terra_tmp) 

# Open cluster environment
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
log_message(paste("Parallel processing set up with", num_cores, "cores"))

# Run the workflow
result <- main()

# Stop parallel processing
stopImplicitCluster()
log_message("Parallel processing stopped")

# # Print out the paths of all extracted files
# log_message("Path to extracted files:")
# print(unlist(all_results))

# Log session info for reproducibility
log_message("Session info:")
print(sessionInfo())


# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#