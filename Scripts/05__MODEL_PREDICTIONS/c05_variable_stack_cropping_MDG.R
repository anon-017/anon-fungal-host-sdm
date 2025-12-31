# ---
# title: "c01_stack_cropping_MDG.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-27"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Crops processed environmental stacks using terra::crop to MDG extent 
# - Exports as raster .tif csv for use in 05__model_predictions.R script
# - Crops final stacks down to template A (MDG) extent for future predictions later in workflow
# 
# HOW
# - Runs code on cluster due to memory constraints locally
# 
# WHY
# - Current and future environmental data/predictor variables need to be cropped to the study area
# - initial code crops to AFR and MDG (templates A + B)
# - this script must be run prior to projections and predictions scripts
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           ----  Workspace set up ----
#==============================================================================#


# individual packages that require installing directly in R interactive session while on cluster
library(terra)
library(doParallel)
library(foreach)
library(future.apply)
library(data.table)


# File paths on cluster and mnt from transfer folder
# mnt = where the renamed, processed, environmental stacked rasters are saved
mnt_folder <- "/mnt/xxx/05__model_prediction/input/variable_stacks" # T drive mount on cluster where data is stored
terra_tmp <- "/import/xxx/05__model_predictionw/tmp"
templateA_path <- "/mnt/xxx/00__species_occurrence_cleaning/input/AOI/MDG/MDG_aea.tif"

# Create output directories
mnt_out <- file.path(mnt_folder, "MDG")
dir.create(mnt_out, showWarnings = FALSE, recursive = TRUE) # fully renamed stacked rasters of all env data

# Set up logging
log_file <- file.path(mnt_out, "process_log.txt")
cat("", file = log_file) # Initialize log file

## CRS: African Equal Area ----
crs <- "ESRI:102022"


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

# Function to crop processed environmental raster stacks to MDG / template A extent for predictions
crop_envdata_MDG <- function(stack_file, templateA) {
  # Read the stack
  stacked <- terra::rast(stack_file)
  
  # Read templateA
  templateA <- terra::rast(templateA_path)
  
  # Get current names
  current_names <- names(stacked)
  
  # Create a list to store cropped layers
  cropped_layers <- list()
  
  # Process each layer name
  for(i in seq_along(current_names)) {
    # Extract and crop each layer individually
    layer <- stacked[[i]]
    cropped_layers[[i]] <- terra::crop(layer, templateA)
  }
  
  # Combine all cropped layers into a new stack
  cropped_stack <- terra::rast(cropped_layers)
  
  # Save back with new names
  output_file <- gsub("\\.tif$", "_MDG.tif", stack_file)
  terra::writeRaster(cropped_stack, output_file, overwrite = TRUE)
  
  return(list(
    output_file = output_file
  ))
}


#==============================================================================#
#                   1. Main execution function----
#==============================================================================#


# Set up parallel processing
num_cores <- min(5, parallel::detectCores() - 1)  # Use 5 cores or less if not available
cl <- makeCluster(num_cores)
registerDoParallel(cl)
log_message(paste("Parallel processing set up with", num_cores, "cores"))

# Memory management - adjust for 5 cores
memory_per_core <- (128 * 1024) / num_cores  # Dividing 128GB by number of cores
options(future.globals.maxSize = memory_per_core * 1024^2)  # Set memory limit per core
terraOptions(memfrac = 0.8/num_cores, tempdir = terra_tmp)  # Adjust terra memory fraction


# The main workflow function for this script
main <- function() { # rename the processed raster stacks
  log_message("Starting cropping workflow")
  
  # List all stack files
  stack_files <- list.files(mnt_folder, pattern = "\\.tif$", recursive=T, full.names = TRUE)
  log_message(sprintf("Found %d stacked files to crop to MDG", length(stack_files)))
  
  # Process files in parallel giving each worker all functions and packages needed
  results <- foreach(
    stack_file = stack_files,
    .packages = c("terra"),
    .export = c("crop_envdata_MDG", "templateA_path", "mnt_out", "log_message", "log_file"),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch({
      # Crop the stack
      result <- crop_envdata_MDG(stack_file, templateA_path)
      
      # Move renamed file to final destination
      final_path <- file.path(mnt_out, basename(result$output_file))
      file.copy(result$output_file, final_path, overwrite = TRUE)
      file.remove(result$output_file)  # Remove the temporary file
      
      log_message(sprintf("Successfully cropped %s", basename(stack_file)))
      
      
      # Return results
      list(
        input_file = stack_file,
        output_file = final_path
      )
    }, error = function(e) {
      log_message(sprintf("Error processing %s: %s", stack_file, e$message))
      list(
        input_file = stack_file,
        error = e$message
      )
    })
  }
  
  return(results)
}


#==============================================================================#
#                          ---- Run workflow ----
#==============================================================================#


# Run the workflow
results <- main()  # Changed from result to results to match later references

# Stop parallel processing
stopCluster(cl)
stopImplicitCluster()
log_message("Parallel processing stopped")

# Print out summary of results
successful <- sum(sapply(results, function(x) !is.null(x) && !("error" %in% names(x))))
failed <- sum(sapply(results, function(x) !is.null(x) && ("error" %in% names(x))))
log_message(sprintf("Processing complete: %d files successful, %d files failed", successful, failed))

# Log session info for reproducibility
log_message("Session info:")
print(sessionInfo())

# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#