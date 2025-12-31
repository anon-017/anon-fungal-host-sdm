# ---
# title: "01_single_env_stack_extraction.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-17"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Extracts environmental data for the combined presence-background points exported from scripts:
#    -- 00__species_occurrence_cleaning_bg.R (combined occurrence outputs)
#    -- c01__extract_species_enviro_data.R (only the "current" stacked env_vars output)
# - Renames/shortens environmental predictor variables for use within SDM process
# - Exports as csv for use in next script "02__variable_selection.R"
# 
# HOW
# - Takes "current" environmental stack (only) output from cluster stack output c01__script
# - Uses terra::extract for each species occurrence location (presence + background samples)
# - Saves csvs to file (in 01_ output and 02__input)
# 
# WHY
# - Current environmental data/predictor variables need to be tested for collinearity (next script)
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#

#### Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

#### Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

# Remove function that are not required for this script
rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, 
   align.forest.raster, align.raster, validate.processed.climate.files, calc_buffer_distance, calculate_metrics,
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, random_bg_whole_area)

#### CRS: African Equal Area ----
crs <- "ESRI:102022"

#### Environmental data look up ----
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
#                       1. Folders and directories ----
#==============================================================================#

#### Base folders for existing data ----
#datafolder <- file.path("//xxx")

# Data folder relating to script 00 (species occurrences)
datafolder00 <- file.path(datafolder, "00__species_occurrence_cleaning")

# Data folder relating to script 01 (env data extraction)
datafolder01 <- file.path(datafolder, "01__extract_species_enviro_data")

# Data folder relating to script 02 (collinearity)
datafolder02 <- file.path(datafolder, "02__variable_selection")


#### Create output folders ----
# List of output folders
out_dirs <- list(
  out_dir__01 = file.path(datafolder01, "output", "extracted", "current"),
  out_dir__02 = file.path(datafolder02, "input", "extracted_envdata", "current")
)

# current environmental data stack path
current_path <- file.path(datafolder01, "output", "stacked", "current_1981-2010.tif")


#==============================================================================#
#                     2. Post-process environmental stacks ----
#==============================================================================#

## Load existing environmental stacks

stack_files <- list.files(file.path(datafolder01, "output", "stacked"), pattern = "\\.tif$", full.names = TRUE)
for(stack_file in stack_files) {
  result <- rename_stack_layers(stack_file, bio_lookup)
  print(paste("Processed:", stack_file))
  print("Name changes:")
  print(result$name_changes)
}


#==============================================================================#
#                         3. Load species occurrences ----
#==============================================================================#

## combined presence-background occurrences
# 4 sets of combined species presence-background occurrences
# 2 per spp (ecors, 2far and x2 whole area random)
comb_occs_list <- list.files(file.path(datafolder00, "output", "combined_occs"), full.names=T)

# Function to load each species presence-background occurrence combo
load_species_data <- function(file_paths) {
  species_dfs <- lapply(file_paths, function(path) {
    e <- new.env()
    load(path, envir = e)
    e[[ls(e)[1]]]
  })
  names(species_dfs) <- basename(file_paths)
  return(species_dfs)
} # End function load_species_data

# Use "load_species_data" function to load species data for each set of combined presence-background point sets
species_data <- load_species_data(comb_occs_list)


#==============================================================================#
#                      4. Tidy environmental data  ----
#==============================================================================#

# Load "current" environmental data stack
current_env <- rast(current_path)

# Renaming and tidying of environmental data
# Function to rename environmental variables
rename_env_columns <- function(df) {
  old_names <- names(df)
  new_names <- old_names
  
  # Rename and shorten names of bioclimatic variables (column names)
  for(i in 1:19) {
    old_pattern <- paste0("temp_bio", i, "_1981-2010_historical_processed")
    new_pattern <- paste0("bio", i)
    new_names[old_names == old_pattern] <- new_pattern
  }
  
  # Rename other variables (non-CHELSA variables)
  other_patterns <- c(
    "fcc_123_AFR_aea" = "forest", # forest
    "aspect" = "aspect", # elevation/aspect
    "scenario" = "scenario", # climate scenario
    "year_range" = "year_range" # current or future temporal year range
  )
  
  # Use patterns to rename
  new_names[old_names %in% names(other_patterns)] <- 
    other_patterns[old_names[old_names %in% names(other_patterns)]]
  
  names(df) <- new_names
  return(df)
} # End function rename_env_columns


#==============================================================================#
#                   5. Extract current env data for each spp ----
#==============================================================================#

# Loop to extract current env data for each species
# Process each species occurrence file
for (sp_name in names(species_data)) {
    # Get current species dataframe
    sp_df <- species_data[[sp_name]]
    
    # Create spat vector from XY coordinates (in ESRI:102022 projection)
    sp_vect <- terra::vect(sp_df, geom= c("x", "y"), crs=crs)

    # Extract values from all layers for these points
    env_values <- terra::extract(current_env, sp_vect)
    
    # Combine with original coordinates
    extract_df <- cbind(sp_df[, c("x", "y", "occ")],
                       as.data.frame(env_values))
    
    # Remove any NAs
    extracted_df <- extract_df[complete.cases(extract_df), ]
  
    # Tidy environmental data based on "rename_env_columns" function
    extracted_df <- rename_env_columns(extracted_df)
    
    # Output file name based on input name
    output_name <- paste0("env_extract_", 
                              tools::file_path_sans_ext(sp_name),
                              ".csv")
    # Write to both output directories
    for (out_dir in out_dirs) {
      output_path <- file.path(out_dir, output_name)
      write.csv(extracted_df, file = output_path, row.names = FALSE)
    }
}

# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#