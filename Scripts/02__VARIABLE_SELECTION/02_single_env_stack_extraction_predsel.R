# ---
# title: "02_single_env_stack_extraction_predsel.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-17"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Extracts environmental data for the combined presence-background points exported from scripts:
# 
# HOW
# - Takes "current" environmental stack (only) output from cluster stack output c01__script
# - extracts for each selected predicted variable set from output of "02__variable_selection_by_species.R"
# - Uses terra::extract for each species occurrence location (presence + background samples)
# - Saves to file (in 02_ output and 03__input)
# 
# WHY
# - Current environmental data/predictor variables need to be tested for collinearity (next script)
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           ----  0. Workspace set up ----
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
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, random_bg_whole_area)

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
#                      ----  1. Folders and directories ----
#==============================================================================#

### Base folders for existing data ----
datafolder <- file.path("//xxx")

# Base data folder relating to script 00 (species occurrences)
datafolder00 <- file.path(datafolder, "00__species_occurrence_cleaning")
# Base data folder relating to script 01 (first env data extraction)
datafolder01 <- file.path(datafolder, "01__extract_species_enviro_data")
# Base data folder relating to script 02 (variable selection)
datafolder02 <- file.path(datafolder, "02__variable_selection")

### Load existing species and environmental data from previous scripts ----
# datafolder01 output folder: List csv of combined species/enviro occs (01_single_env_stack_extraction.R)
c_occs <- list.files(file.path(datafolder01, "output", "extracted", "current"), pattern = "calo", full.names = T)
vl_occs <- list.files(file.path(datafolder01, "output", "extracted", "current"), pattern = "vertlept", full.names = T)

# datafolder02 output folder: containing final predictor variable sets script (02__variable_selection_by_species.R)
csdm_vars_dir <- file.path(datafolder02, "output", "selected_vars", "sdm_variable_sets", "calo") 
vsdm_vars_dir <- file.path(datafolder02, "output", "selected_vars", "sdm_variable_sets", "vertlept") 

### Create output directories ----
# Save outputs in data folder relating to 03 (model building)
datafolder03 <- file.path(datafolder, "03__model_training", "input", "final_predictor_sets")
dir.create(datafolder03, recursive = TRUE, showWarnings = FALSE) 

### File path to "current" environmental raster stack ----
current_path <- file.path(datafolder01, "output", "stacked", "current_1981-2010_renamed.tif")


#==============================================================================#
#                      ----  2. Extraction function ----
#==============================================================================#


# Function to extract environmental data for each:
  # background method - 2 per spp:(CALO: ecor, whole; vertlept: 2far, whole)
  # correlation threshold, - 2 x 2 per spp (0.5, 0.7, per background method)
  # expert selection (Calo only)
  # selected predictor set - 2 x 2 x 4 per spp
  # (base = selected max vars (via thresholding and VIF scores)), +  forest, + aspect,  & forest &  aspect)


# Function to extract environmental data for each selected predictor variable set
extract_environmental_data <- function(occurrence_files, 
                                       env_raster_path,
                                       sdm_vars_dir,
                                       output_dir,
                                       species_name,
                                       crs) {
  
  # Validate species name - standardize case but maintain original for output
  species_display <- species_name  # Keep original for display/output
  species_name <- tolower(species_name)
  
  # Validate against known species
  valid_species <- c("calo", "vertlept")
  if(!species_name %in% valid_species) {
    stop("Invalid species name. Must be one of: ", paste(valid_species, collapse=", "))
  }
  
  # Define species-specific background methods
  species_bg_methods <- list(
    calo = c("whole", "ecor"),
    vertlept = c("whole", "2far")
  )
  
  # Define expected selection types
  selection_types <- c("thresh_0.5", "thresh_0.7", "expert")
  
  # Create output directory structure
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  log_dir <- file.path(output_dir, "logs")
  dir.create(log_dir, showWarnings = FALSE)
  
  # Load environmental raster stack
  cat("Loading environmental raster stack from:", env_raster_path, "\n")
  current_env <- rast(env_raster_path)
  
  # Filter occurrence files for this species
  species_files <- data.frame(
    filepath = occurrence_files,
    filename = basename(occurrence_files),
    stringsAsFactors = FALSE
  )
  
  # Filter for the specified species (case insensitive)
  species_files <- species_files[grepl(species_name, tolower(species_files$filename)), ]
  
  if(nrow(species_files) == 0) {
    stop("No occurrence files found for species: ", species_name)
  }
  
  # Extract background methods from filenames
  species_files$bg_method <- sapply(species_files$filename, function(x) {
    x_lower <- tolower(x)
    bg_methods <- species_bg_methods[[species_name]]
    
    for(method in bg_methods) {
      if(grepl(method, x_lower)) return(method)
    }
    return(NA)
  })
  
  # Remove any files where we couldn't identify the background method
  species_files <- species_files[!is.na(species_files$bg_method), ]
  
  if(nrow(species_files) == 0) {
    stop("No files with valid background methods found for species: ", species_name)
  }
  
  cat("Found background methods:", paste(species_files$bg_method, collapse = ", "), "\n")
  
  # Check for variable set files in species subdirectory
  cat("\nSearching for variable set files for species:", species_name, "\n")
  
  # Look for files in species subdirectory with lowercase name
  species_sdm_dir <- file.path(sdm_vars_dir, species_name)
  
  if(dir.exists(species_sdm_dir)) {
    cat("Looking in species directory:", species_sdm_dir, "\n")
    var_files <- list.files(species_sdm_dir, 
                            pattern = "\\.rds$", 
                            full.names = TRUE)
    
    # List contents to debug
    cat("Files found in directory:\n")
    all_files <- list.files(species_sdm_dir)
    if(length(all_files) > 0) {
      cat(paste(" -", all_files), sep = "\n")
    } else {
      cat("  (directory is empty)\n")
    }
  } else {
    cat("Species directory not found:", species_sdm_dir, "\n")
    
    # Try with capitalized name
    capitalized_name <- paste0(toupper(substr(species_name, 1, 1)), 
                               substr(species_name, 2, nchar(species_name)))
    capitalized_dir <- file.path(sdm_vars_dir, capitalized_name)
    
    if(dir.exists(capitalized_dir)) {
      cat("Trying capitalized directory:", capitalized_dir, "\n")
      var_files <- list.files(capitalized_dir, 
                              pattern = "\\.rds$", 
                              full.names = TRUE)
    } else {
      cat("Looking in main directory for pattern:", paste0(species_name, "_"), "\n")
      var_files <- list.files(sdm_vars_dir, 
                              pattern = paste0(species_name, "_.*\\.rds$"), 
                              full.names = TRUE)
    }
  }
  
  # Remove sdm_formula files and summary files
  var_files <- var_files[!grepl("(sdm_formula|summary)", var_files)]
  
  if(length(var_files) == 0) {
    stop("No variable set files found for species: ", species_name, 
         "\nChecked in: ", species_sdm_dir)
  }
  
  cat("Found", length(var_files), "variable set files\n")
  cat(paste(" -", basename(var_files)), sep = "\n")
  
  # Create tracking dataframe
  extraction_log <- data.frame(
    species = character(),
    bg_method = character(),
    threshold = character(),
    var_set = character(),
    n_points = numeric(),
    n_variables = numeric(),
    output_file = character(),
    stringsAsFactors = FALSE
  )
  
  # Create output directory for this species
  species_out_dir <- file.path(output_dir, species_name)
  dir.create(species_out_dir, showWarnings = FALSE)
  
  # Process each background method
  for(i in 1:nrow(species_files)) {
    bg_method <- species_files$bg_method[i]
    cat("\nProcessing background method:", bg_method, "\n")
    
    # Read CSV file
    cat("Loading:", species_files$filename[i], "\n")
    sp_df <- read.csv(species_files$filepath[i])
    
    # Create spatial vector
    sp_vect <- terra::vect(sp_df, geom = c("x", "y"), crs = crs)
    
    # Extract all environmental values once
    cat("Extracting environmental values...\n")
    env_values <- terra::extract(current_env, sp_vect)
    
    # Combine with coordinates and presence/absence
    # Ensure we have ID column (needed for proper matching)
    if(!"ID" %in% names(sp_df)) {
      sp_df$ID <- 1:nrow(sp_df)
    }
    
    # Create base dataframe with consistent column names
    base_df <- cbind(sp_df[, c("x", "y", "occ", "ID")],
                     as.data.frame(env_values))
    
    # Remove any NAs
    base_df <- base_df[complete.cases(base_df), ]
    
    # Standardize forest and aspect variable names
    names(base_df) <- gsub("forest[_ ]cover", "forest", names(base_df), ignore.case = TRUE)
    
    cat("\nAvailable variables in environmental data:\n")
    print(names(base_df))
    
    # Get all files for this background method - be more inclusive
    bg_pattern <- paste0("_", bg_method, "_")
    bg_var_files <- var_files[grepl(bg_pattern, basename(var_files))]
    
    # Debug: print all found files for this bg method
    cat("\nFiles found for background method", bg_method, ":\n")
    cat(paste(" -", basename(bg_var_files)), sep = "\n")
    
    # Process each variable set
    for(var_file in bg_var_files) {
      # Parse filename components
      var_filename <- basename(var_file)
      file_parts <- tools::file_path_sans_ext(var_filename)
      components <- strsplit(file_parts, "_")[[1]]
      
      # Debug: print components for each file
      cat("\nProcessing file:", var_filename, "\n")
      cat("File components:", paste(components, collapse=" | "), "\n")
      
      # Extract threshold/selection approach - modified to handle both patterns
      if(any(grepl("expert", components))) {
        threshold <- "expert"
        cat("Detected expert selection\n")
      } else if(any(grepl("thresh", components))) {
        thresh_idx <- which(grepl("thresh", components))[1]
        if(thresh_idx + 1 <= length(components)) {
          threshold <- paste(components[thresh_idx], components[thresh_idx + 1], sep="_")
          cat("Detected threshold selection:", threshold, "\n")
        } else {
          cat("Malformed threshold in filename:", var_filename, "\n")
          next
        }
      } else {
        cat("Cannot determine selection approach for file:", var_filename, "\n")
        next
      }
      
      # Determine variable set (base, forest, aspect, forest_aspect)
      # Search through all components
      has_forest <- any(grepl("forest", components))
      has_aspect <- any(grepl("aspect", components))
      has_forest_aspect <- any(grepl("forest_aspect", components))
      
      if(has_forest_aspect) {
        var_set <- "forest_aspect"
      } else if(has_forest && has_aspect) {
        var_set <- "forest_aspect"
      } else if(has_forest) {
        var_set <- "forest"
      } else if(has_aspect) {
        var_set <- "aspect"
      } else {
        var_set <- "base"
      }
      
      cat("Determined variable set:", var_set, "\n")
      
      # Load variables
      tryCatch({
        variables <- readRDS(var_file)
        
        # Check if we got the expected format
        if(!is.character(variables)) {
          cat("WARNING: Unexpected format in file:", var_filename, "\n")
          cat("Expected character vector, got:", class(variables), "\n")
          next
        }
        
        # Print loaded variables
        cat("Loaded variables:", paste(variables, collapse=", "), "\n")
        
        # Standardize variable names for forest and aspect
        variables <- gsub("forest[_ ]cover", "forest", variables, ignore.case = TRUE)
        
        cat("\nChecking variables for", var_set, ":\n")
        cat("Required:", paste(variables, collapse = ", "), "\n")
        cat("Available:", paste(names(base_df), collapse = ", "), "\n")
        
        # Check if all variables exist in the data
        missing_vars <- setdiff(variables, names(base_df))
        if(length(missing_vars) > 0) {
          warning("Missing variables for ", var_set, ": ", 
                  paste(missing_vars, collapse = ", "))
          next
        }
        
        cat("All required variables found for", var_set, "\n")
        
        # Create output dataframe with only needed variables
        output_cols <- c("x", "y", "occ", "ID", variables)
        output_df <- base_df[, output_cols]
        
        # Create output filename
        output_name <- paste(species_name, bg_method, threshold, var_set, "dataset.csv", 
                             sep = "_")
        output_path <- file.path(species_out_dir, output_name)
        
        # Save to CSV
        write.csv(output_df, file = output_path, row.names = FALSE)
        cat("Saved dataset to:", output_path, "\n")
        
        # Update tracking
        extraction_log <- rbind(extraction_log, data.frame(
          species = species_name,
          bg_method = bg_method,
          threshold = threshold,
          var_set = var_set,
          n_points = nrow(output_df),
          n_variables = length(variables),
          output_file = output_path,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        cat("ERROR processing file:", var_filename, "\n")
        cat("Error message:", e$message, "\n")
      })
    }
  }
  
  # Check if we got all expected combinations
  expected_bg_methods <- species_bg_methods[[species_name]]
  expected_thresholds <- c("thresh_0.5", "thresh_0.7", "expert")
  expected_var_sets <- c("base", "forest", "aspect", "forest_aspect")
  
  expected_combos <- expand.grid(
    bg_method = expected_bg_methods,
    threshold = expected_thresholds,
    var_set = expected_var_sets,
    stringsAsFactors = FALSE
  )
  
  # Check which combinations are missing
  actual_combos <- extraction_log[, c("bg_method", "threshold", "var_set")]
  missing_combos <- anti_join(expected_combos, actual_combos)
  
  if(nrow(missing_combos) > 0) {
    cat("\nWARNING: Some expected combinations were not processed:\n")
    for(i in 1:nrow(missing_combos)) {
      cat("- Missing:", 
          missing_combos$bg_method[i], 
          missing_combos$threshold[i], 
          missing_combos$var_set[i], "\n")
    }
  }
  
  # Count by selection type
  sel_counts <- table(extraction_log$threshold)
  cat("\nExtractions by selection type:\n")
  print(sel_counts)
  
  # Save extraction log with timestamp
  log_filename <- paste0(species_name, "_extraction_log_", 
                         format(Sys.time(), "%Y%m%d_%H%M"), ".rds")
  saveRDS(extraction_log, file.path(log_dir, log_filename))
  
  # Create summary report
  summary_filename <- paste0(species_name, "_extraction_summary_", 
                             format(Sys.time(), "%Y%m%d_%H%M"), ".txt")
  sink(file.path(log_dir, summary_filename))
  cat("Environmental Data Extraction Summary for", species_display, "\n")
  cat(paste(rep("=", 50), collapse=""), "\n\n")
  
  cat("Overall Statistics:\n")
  cat("-----------------\n")
  cat("Total extractions:", nrow(extraction_log), "\n")
  cat("Background methods found:", paste(unique(extraction_log$bg_method), collapse = ", "), "\n")
  cat("Thresholds processed:", paste(unique(extraction_log$threshold), collapse = ", "), "\n")
  cat("Variable sets processed:", paste(unique(extraction_log$var_set), collapse = ", "), "\n\n")
  
  # Per-background method summary
  for(bg in unique(extraction_log$bg_method)) {
    bg_data <- extraction_log[extraction_log$bg_method == bg, ]
    cat("\nBackground method:", bg, "\n")
    cat("  Variable sets processed:", length(unique(bg_data$var_set)), "\n")
    cat("  Thresholds processed:", paste(unique(bg_data$threshold), collapse = ", "), "\n")
    cat("  Number of points:", unique(bg_data$n_points), "\n")
    cat("  Average variables per set:", mean(bg_data$n_variables), "\n")
  }
  
  # Per-threshold summary
  cat("\n\nPer-threshold Statistics:\n")
  cat("----------------------\n")
  for(thresh in unique(extraction_log$threshold)) {
    thresh_data <- extraction_log[extraction_log$threshold == thresh, ]
    cat("\nThreshold:", thresh, "\n")
    cat("  Background methods:", paste(unique(thresh_data$bg_method), collapse = ", "), "\n")
    cat("  Variable sets:", paste(unique(thresh_data$var_set), collapse = ", "), "\n")
    cat("  Average variables:", mean(thresh_data$n_variables), "\n")
  }
  sink()
  
  return(extraction_log)
}


#==============================================================================#
#                 ---- 3. Run extraction on predictor variable sets ----
#==============================================================================#


# For calo
calo_extraction <- extract_environmental_data(
  occurrence_files = list.files(file.path(datafolder01, "output", "extracted", "current"), pattern = "calo", full.names = T),
  env_raster_path = current_path, # current environmental data stack
  sdm_vars_dir = csdm_vars_dir, # variable sets from 02__outputs
  output_dir = datafolder03, # output folder
  species_name = "calo", # naming to pick up only relevant files
  crs = "ESRI:102022" # crs for project
)

gc()

# For vertlept
vertlept_extraction <- extract_environmental_data(
  occurrence_files = list.files(file.path(datafolder01, "output", "extracted", "current"), pattern = "vertlept", full.names = T),
  env_raster_path = current_path, # current environmental data stack
  sdm_vars_dir = vsdm_vars_dir, # variable sets from 02__outputs
  output_dir = , datafolder03, # output folder
  species_name = "vertlept", # naming to pick up only relevant files
  crs = "ESRI:102022" # crs for project
)


# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#