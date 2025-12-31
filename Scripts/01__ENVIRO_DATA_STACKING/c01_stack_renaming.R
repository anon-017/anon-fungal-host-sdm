# ---
# title: "c01_stack_renaming.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# update: "2025-12-31"
# ---


# WORK FLOW NOTES:

# WHAT
# - Renames processed environmental stacks using biolook up
# - Exports as csv for use in next script "02__variable_selection.R"
# 
# HOW
# - Runs code on cluster due to memory constraints locally
# 
# WHY
# - Current environmental data/predictor variables need to be consistently named for whole workflow
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
mnt_folder <- "/mnt/xxx/01__extract_species_enviro_data" # T drive mount on cluster where data is stored
import_dir <- "/import/xxx/01__extract_species_enviro_data" # cluster folder where AOI and template is stored (updated for ecoc9)
input_dir <- file.path(mnt_folder, "output", "stacked")
terra_tmp <- "/import/xxx/01__extract_species_enviro_data/tmp"

# Create output directories
import <- "/import/xxx/01__extract_species_enviro_data/output"
output_dir <- file.path(import, "18_02_2025_rename_stacks")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # fully renamed stacked rasters of all env data
 
# limit memory
memory.limit(size = 128 * 1024)  # Set memory limit to 128GB

## CRS: African Equal Area ----
crs <- "ESRI:102022"

## Environmental data look up ----
bio_lookup <- c(
  x = "x",                        # longitude (crs = "ESRI:102022")
  y = "y",                        # latitude (crs = "ESRI:102022")
  occ = "occ",                    # binary occurrence value (1 = presence, 0 = background point)
  ID = "ID",                      # cell ID of extracted environmental grid
  bio1 = "Ann Temp",              # man annual temperature in degrees Celsius
  bio2 = "Diurnal Range",         # mean diurnal range (mean monthly (max temp - min temp))
  bio3 = "Isothermality",         # isothermality (bio2/bio7) (* 100)
  bio4 = "Temp Seasonality",      # temperature seasonality (standard deviation *100)
  bio5 = "Max Temp Warmest",      # max temperature of warmest month
  bio6 = "Min Temp Coldest",      # min temperature of coldest month
  bio7 = "Ann Temp Range",        # temperature annual range (bio5-bio6)
  bio8 = "Mean Temp Wet Q",       # mean temperature of wettest quarter
  bio9 = "Mean Temp Dry Q",       # mean temperature of driest quarter
  bio10 = "Mean Temp Warm Q",     # mean temperature of warmest quarter
  bio11 = "Mean Temp Cold Q",     # mean temperature of coldest quarter
  bio12 = "Ann Precip",           # annual precipitation
  bio13 = "Precip Wet Month",     # precipitation of wettest month
  bio14 = "Precip Dry Month",     # precipitation of driest month
  bio15 = "Precip Seasonality",   # precipitation seasonality (coefficient of variation)
  bio16 = "Precip Wet Q",         # precipitation of wettest quarter
  bio17 = "Precip Dry Q",         # precipitation of driest quarter
  bio18 = "Precip Warm Q",        # precipitation of warmest quarter
  bio19 = "Precip Cold Q",        # precipitation of coldest quarter
  forest = "forest Cover",        # forest cover 2010
  aspect = "aspect",              # aspect (derived from elevation, as a proxy for microclimate, measured in degrees)
  scenario = "scenario",          # CHELSA shared socio-economic pathway (current, or future: ssp126, ssp370, ssp585)
  year_range = "future year"      # future year time span (2011-2040, 2041-2070, 2071-2100)
)


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

# Function to rename the columns (individual processed environmental raster names)
rename_stack_layers <- function(stack_file, bio_lookup) {
  # Read the stack
  stacked <- terra::rast(stack_file)
  
  # Get current names
  current_names <- names(stacked)
  
  # Process each layer name
  for(i in seq_along(current_names)) {
    bio_code <- extract_bio_number(current_names[i])
    if(!is.na(bio_code) && bio_code %in% names(bio_lookup)) {
      names(stacked)[i] <- bio_lookup[bio_code]
    } else {
      warning(sprintf("No lookup found for %s, keeping original name", current_names[i]))
    }
  }
  
  # Save back with new names
  output_file <- gsub("\\.tif$", "_renamed.tif", stack_file)
  terra::writeRaster(stacked, output_file, overwrite = TRUE)
  
  # Log the changes
  name_changes <- data.frame(
    original = current_names,
    new = names(stacked),
    stringsAsFactors = FALSE
  )
  
  return(list(
    output_file = output_file,
    name_changes = name_changes
  ))
}

# Function to extract bio number from environmental layer
extract_bio_number <- function(filename) {
  # For bioclim variables
  bio_num <- regexpr("bio[0-9]{1,2}", filename)
  if (bio_num != -1) {
    bio_code <- substr(filename, bio_num, bio_num + attr(bio_num, "match.length") - 1)
    return(bio_code)
  }
  # For forest cover
  if (grepl("fcc_123", filename)) {
    return("forest")
  }
  # For aspect
  if (grepl("aspect", filename)) {
    return("aspect")
  }
  return(NA_character_)
}


#==============================================================================#
#                       1. Main execution function ----
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
  log_message("Starting main workflow")
  
  # List all stack files
  stack_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
  log_message(sprintf("Found %d stacked files to process", length(stack_files)))
  
  # Process files in parallel giving each worker all functions and packages needed
  results <- foreach(
    stack_file = stack_files,
    .packages = c("terra"),
    .export = c("bio_lookup", "extract_bio_number", "rename_stack_layers"),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch({
      # Rename the stack
      result <- rename_stack_layers(stack_file, bio_lookup)
      
      # Move renamed file to final destination
      final_path <- file.path(output_dir, "renamed", basename(result$output_file))
      file.rename(result$output_file, final_path)
      
      # Return results
      list(
        input_file = stack_file,
        output_file = final_path,
        name_changes = result$name_changes
      )
    }, error = function(e) {
      log_message(sprintf("Error processing %s: %s", stack_file, e$message))
      NULL
    })
  }
  
  # Save summary of changes after processing
  if (!is.null(results)) {
    all_changes <- do.call(rbind, lapply(results[!sapply(results, is.null)], 
                                         function(x) x$name_changes))
    write.csv(all_changes, 
              file.path(output_dir, "name_changes_summary.csv"), 
              row.names = FALSE)
  }
  
  return(results)
}
  

#==============================================================================#
#                          ---- Run workflow ----
#==============================================================================#


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