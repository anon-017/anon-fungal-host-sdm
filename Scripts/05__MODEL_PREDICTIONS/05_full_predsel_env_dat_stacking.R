# ---
# title: "05_full_predsel_env_dat_stacking.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-27"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES: Creates variable-specific Raster Stacks for SDM Prediction

# WHAT
# - Creates raster stacks of individual final predictor variable sets used to train models
# - Creates for all time frames and climate scenarios
# - Creates 10 stacks per model variable set (1 current, 9 future (3x years x 3 ssps))
#
# HOW
# - Identifies required variables from model training output
# - Creates exact raster stacks with identically named layers
# 
# WHY
# - Exact matching of predictor variables is required for making predictions from trained models
# - Env_dat ready for making predictions in 05__model_predictions.R
# - this script must run before 05__model_prediction.R
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           ----  Workspace set up ----
#==============================================================================#

## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

## Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

## CRS: African Equal Area ----
crs <- "ESRI:102022"

#==============================================================================#
##                      ----  Folders and directories ----
#==============================================================================#

### Base folders for existing data ----
datafolder <- file.path("//xxx")

# Base data folders from various stages of the workflow
datafolder00 <- file.path(datafolder, "00__species_occurrence_cleaning")
datafolder01 <- file.path(datafolder, "01__extract_species_enviro_data")
datafolder02 <- file.path(datafolder, "02__variable_selection")
datafolder03 <- file.path(datafolder, "03__model_training")
datafolder04 <- file.path(datafolder, "04__model_eval_perf")

# Top performing models location
datafolder04_top <- file.path(datafolder04, "output", "top_performing_models")

# Input/output for this script
datafolder05 <- file.path(datafolder, "05__model_prediction")
data05_input <- file.path(datafolder05, "input")
data05_output <- file.path(datafolder05, "output")

# Create output directories
for (dir in c("variable_stacks")) {
  dir.create(file.path(data05_input, dir), recursive = TRUE, showWarnings = FALSE)
}

# Environmental data directories
env_data_path <- file.path(datafolder01, "output", "stacked")

# Variable lookup for column standardization
bio_lookup <- c(
  "bio1" = "Ann.Temp",
  "bio2" = "Diurnal.Range",
  "bio3" = "Isothermality",
  "bio4" = "Temp.Seasonality",
  "bio5" = "Max.Temp.Warmest",
  "bio6" = "Min.Temp.Coldest",
  "bio7" = "Ann.Temp.Range",
  "bio8" = "Mean.Temp.Wet.Q",
  "bio9" = "Mean.Temp.Dry.Q",
  "bio10" = "Mean.Temp.Warm.Q",
  "bio11" = "Mean.Temp.Cold.Q",
  "bio12" = "Ann.Precip",
  "bio13" = "Precip.Wet.Month",
  "bio14" = "Precip.Dry.Month",
  "bio15" = "Precip.Seasonality",
  "bio16" = "Precip.Wet.Q",
  "bio17" = "Precip.Dry.Q",
  "bio18" = "Precip.Warm.Q",
  "bio19" = "Precip.Cold.Q",
  "forest" = "forest",
  "aspect" = "aspect"
)

# Reverse lookup for converting model variable names back to raster names
reverse_lookup <- setNames(names(bio_lookup), bio_lookup)

# Define time periods and climate pathways
time_periods <- c("current", "2011-2040", "2041-2070", "2071-2100")
climate_pathways <- c("current", "ssp126", "ssp370", "ssp585")

# Species list
species_list <- c("calo", "vertlept")


#==============================================================================#
#                 ---- Function to create raster stacks ----
#==============================================================================#


# Function to create prediction stacks with flexible name handling
create_prediction_stacks <- function(species, top_models_file, model_dir, 
                                            env_data_path, output_dir) {
  cat("\n====== Creating prediction stacks for", species, "======\n")
  
  # Create a comprehensive lookup table for variable name mappings
  # This maps between all possible naming conventions for the same variables
  var_lookup <- data.frame(
    # Standard model names (as used in your models)
    model_name = c(
      "Ann.Temp", "Diurnal.Range", "Isothermality", "Temp.Seasonality", 
      "Max.Temp.Warmest", "Min.Temp.Coldest", "Ann.Temp.Range", 
      "Mean.Temp.Wet.Q", "Mean.Temp.Dry.Q", "Mean.Temp.Warm.Q", "Mean.Temp.Cold.Q",
      "Ann.Precip", "Precip.Wet.Month", "Precip.Dry.Month", "Precip.Seasonality",
      "Precip.Wet.Q", "Precip.Dry.Q", "Precip.Warm.Q", "Precip.Cold.Q",
      "forest", "aspect"
    ),
    
    # Standard bio numbers
    bio_number = c(
      "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
      "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15",
      "bio16", "bio17", "bio18", "bio19", "forest", "aspect"
    ),
    
    # Bio numbers with uppercase
    bio_upper = c(
      "Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7", 
      "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
      "Bio16", "Bio17", "Bio18", "Bio19", "Forest", "Aspect"
    ),
    
    # Full names with spaces
    full_name = c(
      "Annual Mean Temperature", "Mean Diurnal Range", "Isothermality", "Temperature Seasonality",
      "Max Temperature of Warmest Month", "Min Temperature of Coldest Month", "Temperature Annual Range",
      "Mean Temperature of Wettest Quarter", "Mean Temperature of Driest Quarter",
      "Mean Temperature of Warmest Quarter", "Mean Temperature of Coldest Quarter",
      "Annual Precipitation", "Precipitation of Wettest Month", "Precipitation of Driest Month",
      "Precipitation Seasonality", "Precipitation of Wettest Quarter", "Precipitation of Driest Quarter",
      "Precipitation of Warmest Quarter", "Precipitation of Coldest Quarter",
      "Forest Cover", "Aspect"
    ),
    
    # Space-separated names (matching your raster layer names)
    space_name = c(
      "Ann Temp", "Diurnal Range", "Isothermality", "Temp Seasonality", 
      "Max Temp Warmest", "Min Temp Coldest", "Ann Temp Range", 
      "Mean Temp Wet Q", "Mean Temp Dry Q", "Mean Temp Warm Q", "Mean Temp Cold Q",
      "Ann Precip", "Precip Wet Month", "Precip Dry Month", "Precip Seasonality",
      "Precip Wet Q", "Precip Dry Q", "Precip Warm Q", "Precip Cold Q",
      "forest Cover", "aspect"
    ),
    
    # Names for forest in different stacks (by year)
    forest_names = c(
      rep(NA, 19),  # NA for non-forest variables
      "fcc_2040_AFR_aea,fcc_2070_AFR_aea,fcc_2100_AFR_aea",  # Comma-separated list of possible forest names
      NA  # NA for aspect
    ),
    
    stringsAsFactors = FALSE
  )
  
  # Load top models
  if(!file.exists(top_models_file)) {
    cat("Top models file not found:", top_models_file, "\n")
    return(NULL)
  }
  
  top_models <- read.csv(top_models_file)
  cat("Found", nrow(top_models), "top models\n")
  
  # Get list of all available raster stacks
  all_stacks <- list.files(env_data_path, pattern="\\.tif$", full.names=TRUE)
  cat("Found", length(all_stacks), "environmental raster stacks\n")
  
  # Process each model
  all_results <- list()
  
  for(i in 1:nrow(top_models)) {
    model_name <- top_models$Model[i]
    cat("\n==== Processing model:", model_name, "====\n")
    
    # Extract model components
    model_parts <- strsplit(model_name, "_")[[1]]
    
    # Get background method, threshold, and variable set
    bg_method <- NULL
    threshold <- NULL
    var_set <- NULL
    
    # Background method
    for(method in c("whole", "ecor", "2far")) {
      if(any(grepl(method, model_parts))) {
        bg_method <- method
        break
      }
    }
    
    # Threshold
    thresh_idx <- which(grepl("thresh", model_parts))
    if(length(thresh_idx) > 0 && thresh_idx[1] + 1 <= length(model_parts)) {
      threshold <- paste(model_parts[thresh_idx[1]], model_parts[thresh_idx[1] + 1], sep="_")
    } else if(any(grepl("expert", model_parts))) {
      threshold <- "expert"
    }
    
    # Variable set
    if(any(grepl("forest_aspect", model_parts))) {
      var_set <- "forest_aspect"
    } else if(any(grepl("aspect", model_parts)) && !any(grepl("forest", model_parts))) {
      var_set <- "aspect"
    } else if(any(grepl("forest", model_parts)) && !any(grepl("aspect", model_parts))) {
      var_set <- "forest"
    } else {
      var_set <- "base"
    }
    
    cat("Background method:", bg_method, "\n")
    cat("Threshold:", threshold, "\n")
    cat("Variable set:", var_set, "\n")
    
    # Load the model to get variables
    model_path <- file.path(model_dir, species, "Algsets", model_name)
    
    if(!file.exists(model_path)) {
      cat("Model file not found:", model_path, "\n")
      next
    }
    
    # Load the model
    model_env <- new.env()
    load(model_path, envir = model_env)
    
    if(!"sdm_model" %in% ls(model_env)) {
      cat("No sdm_model object found in model file\n")
      # Clean up environment 
      rm(model_env)
      gc(verbose = FALSE)
      next
    }
    
    model_obj <- model_env$sdm_model
    
    # Get variables from model object
    variables <- NULL
    
    # Try multiple paths to extract variables
    if(isS4(model_obj) && "models" %in% slotNames(model_obj) && 
       !is.null(model_obj@models[["occ"]]) && 
       !is.null(model_obj@models[["occ"]][["rf"]]) && 
       !is.null(model_obj@models[["occ"]][["rf"]][["1"]]@varImportance) &&
       !is.null(model_obj@models[["occ"]][["rf"]][["1"]]@varImportance[["training"]]@variables)) {
      
      variables <- model_obj@models[["occ"]][["rf"]][["1"]]@varImportance[["training"]]@variables
      cat("Found variables from model varImportance:", paste(variables, collapse=", "), "\n")
    } else if(isS4(model_obj) && "data" %in% slotNames(model_obj)) {
      # Try to extract from data
      data_cols <- names(model_obj@data)
      # Remove non-predictor columns
      exclude_cols <- c("x", "y", "ID", "occ", "id", "X", "Y", "Occ")
      variables <- setdiff(data_cols, exclude_cols)
      cat("Found variables from model data:", paste(variables, collapse=", "), "\n")
    } else if(isS4(model_obj) && "formula" %in% slotNames(model_obj)) {
      # Try to extract from formula
      formula_text <- as.character(model_obj@formula)[3]  # Right side of formula
      # Split by operators and clean up
      formula_parts <- unlist(strsplit(formula_text, "[+:\\-*]"))
      variables <- trimws(formula_parts)
      cat("Found variables from formula:", paste(variables, collapse=", "), "\n")
    } else {
      cat("Could not extract variables from model\n")
      next
    }
    
    if(length(variables) == 0) {
      cat("No variables found for this model\n")
      next
    }
    
    # Process each environmental stack
    for(stack_file in all_stacks) {
      stack_name <- basename(stack_file)
      
      # Determine climate scenario and time period
      climate <- NA
      period <- NA
      
      if(grepl("current", stack_name)) {
        climate <- "current"
        period <- "current"
      } else {
        # Extract climate scenario
        for(scenario in c("ssp126", "ssp370", "ssp585")) {
          if(grepl(scenario, stack_name)) {
            climate <- scenario
            break
          }
        }
        
        # Extract time period
        for(time_period in c("2011-2040", "2041-2070", "2071-2100")) {
          if(grepl(time_period, stack_name)) {
            period <- time_period
            break
          }
        }
      }
      
      if(is.na(climate) || is.na(period)) {
        cat("Could not determine climate scenario or time period for:", stack_name, "\n")
        next
      }
      
      cat("\nProcessing stack:", stack_name, "\n")
      cat("Climate scenario:", climate, "\n")
      cat("Time period:", period, "\n")
      
      # Create output directory
      out_dir <- file.path(output_dir, "variable_stacks", species)
      dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
      
      # Create output filename
      out_name <- paste(species, bg_method, threshold, var_set, 
                        climate, period, ".tif", sep="_")
      out_file <- file.path(out_dir, out_name)
      
      # Load the full stack
      cat("Loading full stack:", stack_file, "\n")
      full_stack <- terra::rast(stack_file)
      
      # Get layer names from the stack
      stack_layers <- names(full_stack)
      cat("Stack layers:", paste(stack_layers, collapse=", "), "\n")
      
      # Create a list to store the layers we'll use
      selected_layers <- list()
      
      # For each model variable, find the corresponding layer in the stack
      for(var in variables) {
        cat("Looking for variable:", var, "\n")
        
        # Find this variable in our lookup table
        var_row <- which(var_lookup$model_name == var)
        
        if(length(var_row) == 0) {
          cat("WARNING: Variable", var, "not found in lookup table\n")
          # Special handling for 'forest' and 'aspect'
          if(var == "forest" || var == "aspect") {
            # Will handle separately
            next
          } else {
            cat("Skipping unknown variable:", var, "\n")
            next
          }
        }
        
        # Get all possible names for this variable
        possible_names <- as.character(var_lookup[var_row, ])
        cat("Possible names:", paste(possible_names, collapse=", "), "\n")
        
        # Check if any of these names exist in the stack
        matching_layers <- stack_layers[stack_layers %in% possible_names]
        
        # For forest, also check the additional forest naming patterns
        if(var == "forest" && length(matching_layers) == 0) {
          # Get the forest variants from the lookup table
          forest_variants <- var_lookup$forest_names[which(var_lookup$model_name == "forest")]
          if(!is.na(forest_variants)) {
            # Split the comma-separated list
            forest_names <- unlist(strsplit(forest_variants, ","))
            # Check for matches
            forest_matches <- stack_layers[stack_layers %in% forest_names]
            if(length(forest_matches) > 0) {
              matching_layers <- forest_matches
            }
          }
        }
        
        if(length(matching_layers) > 0) {
          # Found at least one match
          layer_name <- matching_layers[1]
          cat("Found matching layer:", layer_name, "\n")
          
          # Add this layer to our selection
          selected_layers[[var]] <- full_stack[[layer_name]]
        } else {
          cat("No matching layer found for variable:", var, "\n")
          cat("ERROR: Required variable", var, "not found in stack and no alternative available\n")
          # This is a critical variable, so we can't proceed with this stack
          selected_layers <- list()
          break
        }
      }
      
      # Release the full stack to save memory after selections are made
      full_stack <- NULL
      gc(verbose = FALSE)
      
      # Check if we found all critical variables
      if(length(selected_layers) == 0 || length(selected_layers) < length(variables)) {
        cat("ERROR: Could not find all required variables in the stack. Skipping.\n")
        rm(selected_layers)
        gc(verbose = FALSE)
        next
      }
      
      # Create raster stack from selected layers
      if(length(selected_layers) > 0) {
        # Convert list to raster stack
        model_stack <- terra::rast(selected_layers)
        
        # Final check - ensure all variable names match what the model expects
        final_names <- names(model_stack)
        model_names <- variables[variables %in% names(var_lookup$model_name)]
        
        # Rename layers to match model variables
        new_names <- final_names
        for(i in 1:length(final_names)) {
          current_name <- final_names[i]
          
          # Check if this name is in any column of the lookup table
          for(col in 2:ncol(var_lookup)) {
            if(col == which(colnames(var_lookup) == "forest_names")) {
              # Skip the forest_names column for renaming
              next
            }
            
            matches <- which(var_lookup[,col] == current_name)
            if(length(matches) > 0) {
              # Found a match - get the model name
              new_names[i] <- var_lookup$model_name[matches[1]]
              break
            }
          }
          
          # Handle forest names specially
          if(grepl("fcc_.*_AFR_aea", current_name)) {
            forest_idx <- which(var_lookup$model_name == "forest")
            if(length(forest_idx) > 0) {
              new_names[i] <- "forest"
            }
          }
        }
        
        # Apply new names
        names(model_stack) <- new_names
        cat("Final layer names:", paste(names(model_stack), collapse=", "), "\n")
        
        # Verify that all required variables are present
        missing_vars <- setdiff(variables, names(model_stack))
        if(length(missing_vars) > 0) {
          cat("WARNING: Still missing variables:", paste(missing_vars, collapse=", "), "\n")
          cat("Missing critical variables, skipping this stack\n")
          # Clean up
          rm(model_stack)
          gc(verbose = FALSE)
          next
        }
        
        # Save the stack
        terra::writeRaster(model_stack, out_file, overwrite=TRUE)
        cat("Saved stack to:", out_file, "\n")
        
        # Add to results
        all_results[[length(all_results) + 1]] <- list(
          species = species,
          model_name = model_name,
          bg_method = bg_method,
          threshold = threshold,
          var_set = var_set,
          climate = climate,
          period = period,
          variables = variables,
          found_vars = names(model_stack),
          stack_file = out_file
        )
        
        # Clean up
        rm(model_stack)
        gc(verbose = FALSE)
      }
    }
  }
  
  # Create summary dataframe
  if(length(all_results) == 0) {
    cat("No stacks were created\n")
    return(NULL)
  }
  
  summary_df <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      species = x$species,
      model_name = x$model_name,
      bg_method = x$bg_method,
      threshold = x$threshold,
      var_set = x$var_set,
      climate = x$climate,
      period = x$period,
      num_variables = length(x$variables),
      found_variables = length(x$found_vars),
      stack_file = x$stack_file,
      stringsAsFactors = FALSE
    )
  }))
  
  # Save summary
  summary_file <- file.path(output_dir, "variable_stacks", paste0(species, "_stack_summary.csv"))
  write.csv(summary_df, summary_file, row.names=FALSE)
  cat("\nSaved summary to:", summary_file, "\n")
  
  # Final cleanup
  gc(verbose = FALSE)
  
  # Return summary
  return(summary_df)
}


#==============================================================================#
#                 ---- Function (wrapper - process both species) ----
#==============================================================================#


# Create a wrapper function to process multiple species
process_all_species <- function(species_list, datafolder04_top, datafolder03, 
                                env_data_path, output_dir) {
  all_summaries <- list()
  
  for(species in species_list) {
    cat("\n\n========== PROCESSING SPECIES:", species, "==========\n")
    
    # Process this species
    summary <- create_prediction_stacks(
      species = species,
      top_models_file = file.path(datafolder04_top, paste0(species, "_top3_models.csv")),
      model_dir = file.path(datafolder03, "output"),
      env_data_path = file.path(env_data_path),
      output_dir = output_dir
    )
    
    if(!is.null(summary)) {
      all_summaries[[species]] <- summary
    }
    
    # Major garbage collection between species
    gc(verbose = TRUE)
  }
  
  # Combine summaries
  if(length(all_summaries) > 0) {
    master_summary <- do.call(rbind, all_summaries)
    master_file <- file.path(output_dir, "variable_stacks", "master_stack_summary.csv")
    write.csv(master_summary, master_file, row.names=FALSE)
    cat("\nSaved master summary to:", master_file, "\n")
    return(master_summary)
  } else {
    cat("No summaries were created for any species\n")
    return(NULL)
  }
}


#==============================================================================#
#                           ---- Do the stacking now ----
#==============================================================================#


species_list <- c("calo", "vertlept")
master_summary <- process_all_species(
  species_list = species_list,
  datafolder04_top = file.path(datafolder, "04__model_eval_perf", "output", "top_performing_models"),
  datafolder03 = file.path(datafolder, "03__model_training"),
  env_data_path = file.path(datafolder, "01__extract_species_enviro_data", "output", "stacked"),
  output_dir = file.path(datafolder, "05__model_prediction", "input")
)


#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#

# Clean up
gc()
cat("\nScript completed successfully!\n")