# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-07"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Takes final sdm_variable_sets from 02__variable_selection_by_species.R outputs and build a set of models with 3 ML algorithms
# - Current that is 16 sets of variables per species, see: "03__model_training\input\final_predictor_sets\logs" for summary
# - This set of code uses species occurrences (00), and background points (00), and current environmental predictor variables  for collinearity between predictor variables
# - It outputs a set of SDM package model files
# 
# HOW
# - SDM package is used to build a set of SDM models for two species, and testing different bg methods, correlation thresholds, and predictors
# - this code has been built to run locally (not via Linux cluster)
# 
# WHY
# - Model training and building required before predictions can be made
# - To test which variables are the most important at predicting species presence# 
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (zzELUtemp) and run from there or saved locally to machine TEMP folder but should be deleted
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
   random_bg_whole_area, rename_stack_layers, select07_cv)


## Set up base folders ----

# Base folders
datafolder <- file.path("//xxx")

# Inputs: combined, extracted environmental occurrence data (csv) from script: "02__variable_selection" folders saved to 03__input
datafolder03 <- file.path(datafolder, "03__model_training", "input", "final_predictor_sets")

# Inputs: predictor variable sets (character names)
# This folder contains RDS files of character names of environmental predictors and model formulas
sdm_vars_dir <- file.path(datafolder, "02__variable_selection", "output", "selected_vars", "sdm_variable_sets")

# Outputs: Set of trained models, model performance measures
# 03__ base output model training and performance outputs
data__03_out <- file.path(datafolder, "03__model_training", "output")
dir.create(data__03_out, recursive = TRUE, showWarnings = FALSE)


#==============================================================================#
#                       1. Model building function ----
#==============================================================================#

# For loop that goes round all possible SDM formulas:
# - species
# - bg method
# - threshold
# - base (automatic variable selection based on correlation matrix thresholds and VIF scores (no aspect, no forest))
# - with forest
# - with aspect
# - with forest + aspect

## nested for loop that uses two different algorithm sets:
#   alg3 = "rf", "brt", "maxent"

## All models use a 80:20 training/testing split and weights based on 
# In total this function will produce 32 models per spp. Taking 16 input variable sets and producing models
# with two sets of algorithms per set.

# Function to build SDMs 
# SDM training script 
train_single_sdm <- function(env_file, output_dir, species_name, force_retrain = FALSE) {
  # Define algorithm combinations
  alg_sets <- list(
    alg3 = c("rf", "brt", "maxent")
  )
  
  # Create output directory structure
  species_out_dir <- file.path(output_dir, species_name)
  dir.create(species_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create "Algsets" sub-directory for different algorithm sets
  dir.create(file.path(species_out_dir, "Algsets"), recursive = TRUE, showWarnings = FALSE)
  
  
  # Create results tracking dataframe
  results_tracker <- data.frame(
    species = character(),
    bg_method = character(),
    threshold = character(),
    var_set = character(),
    alg_set = character(),
    model_file = character(),
    n_train_presence = numeric(),
    n_train_background = numeric(),
    n_test_presence = numeric(),
    n_test_background = numeric(),
    variables = character(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  # Extract components from filename
  file_parts <- tools::file_path_sans_ext(basename(env_file))
  parts <- strsplit(file_parts, "_")[[1]]
  
  # Parse filename components - handling different naming patterns
  # Determine species from first part
  curr_species <- tolower(parts[1])
  
  # Extract background method
  bg_method <- parts[2]
  
  # Extract threshold - might be expert or thresh_0.X
  if("expert" %in% parts) {
    threshold <- "expert"
  } else if(any(grepl("thresh", parts))) {
    # Find the thresh part index
    thresh_idx <- which(grepl("thresh", parts))[1]
    # Combine thresh with its value (usually 0.5 or 0.7)
    threshold <- paste(parts[thresh_idx], parts[thresh_idx + 1], sep="_")
  } else {
    # Fallback
    threshold <- paste(parts[3], parts[4], sep="_")
  }
  
  # Determine variable set (base, forest, aspect, forest_aspect)
  if("forest_aspect" %in% parts || (all(c("forest", "aspect") %in% parts) && !("base" %in% parts))) {
    var_set <- "forest_aspect"
  } else if("forest" %in% parts) {
    var_set <- "forest"
  } else if("aspect" %in% parts) {
    var_set <- "aspect"
  } else {
    var_set <- "base"
  }
  
  # Check the correct files / naming
  cat("\nProcessing file:", basename(env_file),
      "\nParsed components:",
      "\n  species:", curr_species,
      "\n  bg_method:", bg_method,
      "\n  threshold:", threshold,
      "\n  var_set:", var_set, "\n")
  
  # Load environmental data
  env_data <- read.csv(env_file)
  cat("\nData loaded:", nrow(env_data), "rows,", ncol(env_data), "columns\n")
  
  # Split presence and background points
  presence_data <- env_data[env_data$occ == 1, ]
  background_data <- env_data[env_data$occ == 0, ]
  
  # Create 80-20 split indices for both presence and background
  set.seed(123)  # for reproducibility
  n_presence <- nrow(presence_data)
  n_background <- nrow(background_data)
  
  train_presence_idx <- sample(n_presence, size = round(0.8 * n_presence))
  train_background_idx <- sample(n_background, size = round(0.8 * n_background))
  
  # Create training and testing datasets
  train_data <- rbind(
    presence_data[train_presence_idx, ],
    background_data[train_background_idx, ]
  )
  test_data <- rbind(
    presence_data[-train_presence_idx, ],
    background_data[-train_background_idx, ]
  )
  
  # Calculate weights for background points
  n_train_presence <- sum(train_data$occ == 1)
  n_train_background <- sum(train_data$occ == 0)
  
  # Weights: make background points sum to same weight as presence points
  presence_weight <- 1
  background_weight <- n_train_presence / n_train_background
  train_weights <- ifelse(train_data$occ == 1, presence_weight, background_weight)
  
  # Get predictor variables - exclude metadata columns
  metadata_cols <- c("x", "y", "occ", "ID")
  variables <- setdiff(names(env_data), metadata_cols)
  cat("Variables:", paste(variables, collapse = ", "), "\n")
  
  # Print sample sizes
  cat("\nTraining set:",
      "\n  Presence points:", sum(train_data$occ == 1),
      "\n  Background points:", sum(train_data$occ == 0),
      "\n  Background weight:", round(background_weight, 4),
      "\nTesting set:",
      "\n  Presence points:", sum(test_data$occ == 1),
      "\n  Background points:", sum(test_data$occ == 0), "\n")
  
  # Save train/test split for later evaluation
  split_info <- list(
    train_data = train_data,
    test_data = test_data,
    train_weights = train_weights
  )
  
  # Only save split info if we don't already have it
  split_info_file <- file.path(species_out_dir, paste0(file_parts, "_split_info.rds"))
  if(!file.exists(split_info_file)) {
    saveRDS(split_info, split_info_file)
    cat("Saved train/test split information to:", split_info_file, "\n")
  } else {
    cat("Using existing train/test split information from:", split_info_file, "\n")
  }
  
  # Process each algorithm set
  for(alg_name in names(alg_sets)) {
    tryCatch({
      # Create unique model identifier
      model_id <- paste(species_name, bg_method, threshold, var_set, alg_name, sep="_")
      
      # Check if model already exists
      existing_model <- check_existing_models(species_out_dir, model_id)
      
      # Skip alg3 models if they already exist and we're not forcing a retrain
      if(alg_name == "alg3" && !force_retrain && existing_model != FALSE) {
        cat("\nSkipping existing alg3 model:", existing_model, "\n")
        
        # Add to results tracker with 'existing' status
        results_tracker <- rbind(results_tracker, data.frame(
          species = species_name,
          bg_method = bg_method,
          threshold = threshold,
          var_set = var_set,
          alg_set = alg_name,
          model_file = existing_model,
          n_train_presence = sum(train_data$occ == 1),
          n_train_background = sum(train_data$occ == 0),
          n_test_presence = sum(test_data$occ == 1),
          n_test_background = sum(test_data$occ == 0),
          variables = paste(variables, collapse = ", "),
          status = "existing",
          stringsAsFactors = FALSE
        ))
        
        next  # Skip to next iteration
      }
      
      # If we reach here, we need to train the model
      cat("\nFitting model:", model_id, "\n")
      
      # Create standard linear formula for ML algorithms
      sdm_formula_linear <- as.formula(paste("occ ~", paste(variables, collapse = " + ")))
      
      # Create SDM data object for training
      d <- sdmData(formula = sdm_formula_linear, train = train_data)
      
        # For alg3 (pure ML algorithms), follow Barbet-Massin (2012) recommendation:
        # - Use equal weights for presence/absence
        # - Replicate 10 times with bootstrapping
        
        # Check if number of presences is very different from number of absences
        # If yes, consider subsampling background points to match presences
        if(n_train_presence * 3 < n_train_background) {
          cat("Large imbalance detected - subsetting background points for ML algorithms\n")
          
          # Create balanced subset for ML methods as per Barbet-Massin
          ml_train_data <- train_data
          bg_indices <- which(ml_train_data$occ == 0)
          selected_bg <- sample(bg_indices, size = n_train_presence, replace = FALSE)
          
          # Keep only selected background points and all presence points
          pres_indices <- which(ml_train_data$occ == 1)
          ml_train_data <- ml_train_data[c(pres_indices, selected_bg), ]
          
          # Create new SDM data object with balanced data
          d <- sdmData(formula = sdm_formula_linear, train = ml_train_data)
          
          # Equal weights since data is now balanced
          ml_weights <- rep(1, nrow(ml_train_data))
        } else {
          # Use the calculated weights if not extremely imbalanced
          ml_weights <- train_weights
        }
        
        # Fit ML models with 10 bootstrap replicates
        sdm_model <- sdm(
          formula = sdm_formula_linear,
          data = d,
          methods = alg_sets[[alg_name]],
          replication = "boot",
          n = 10,
          weights = ml_weights
        )
      
      # Save model in algorithm-specific directory
      model_file <- file.path(species_out_dir, paste0(model_id, ".sdm"))
      save(sdm_model, file = model_file)
      cat("Model saved to:", model_file, "\n")
      
      # Update results tracker
      results_tracker <- rbind(results_tracker, data.frame(
        species = species_name,
        bg_method = bg_method,
        threshold = threshold,
        var_set = var_set,
        alg_set = alg_name,
        model_file = model_file,
        n_train_presence = sum(train_data$occ == 1),
        n_train_background = sum(train_data$occ == 0),
        n_test_presence = sum(test_data$occ == 1),
        n_test_background = sum(test_data$occ == 0),
        variables = paste(variables, collapse = ", "),
        status = "success",
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("\nError fitting model:", model_id, "\n")
      cat("Error message:", e$message, "\n")
      
      results_tracker <- rbind(results_tracker, data.frame(
        species = species_name,
        bg_method = bg_method,
        threshold = threshold,
        var_set = var_set,
        alg_set = alg_name,
        model_file = NA,
        n_train_presence = sum(train_data$occ == 1),
        n_train_background = sum(train_data$occ == 0),
        n_test_presence = sum(test_data$occ == 1),
        n_test_background = sum(test_data$occ == 0),
        variables = paste(variables, collapse = ", "),
        status = paste("error:", e$message),
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Save results log
  results_file <- file.path(species_out_dir, "Algsets", paste0(species_name, "_model_results.rds"))
  saveRDS(results_tracker, results_file)
  
  return(results_tracker)
}


# Wrapper function to process all files for a species
train_all_sdm <- function(env_tables_dir, output_dir, species_name, force_retrain = FALSE) {
  # Validate species
  valid_species <- c("calo", "vertlept")
  if(!tolower(species_name) %in% valid_species) {
    stop("Invalid species name. Must be one of: ", paste(valid_species, collapse=", "))
  }
  
  # Expected background methods per species
  bg_methods <- list(
    calo = c("ecor", "whole"),
    vertlept = c("whole", "2far")
  )
  
  # Get all environmental tables for this species
  env_files <- list.files(
    path = file.path(env_tables_dir, species_name),
    pattern = "_dataset\\.csv$",  # Changed pattern to match our naming
    full.names = TRUE
  )
  
  if(length(env_files) == 0) {
    cat("No files found with _dataset.csv pattern, checking for _env.csv pattern...\n")
    env_files <- list.files(
      path = file.path(env_tables_dir, species_name),
      pattern = "_env\\.csv$",
      full.names = TRUE
    )
    
    if(length(env_files) == 0) {
      stop("No environmental tables found for species: ", species_name)
    }
  }
  
  # Validate that we have expected background methods
  expected_bg_methods <- bg_methods[[tolower(species_name)]]
  file_bg_methods <- sapply(env_files, function(f) {
    parts <- strsplit(basename(f), "_")[[1]]
    parts[2]  # Background method is typically the second part
  })
  
  if(!all(file_bg_methods %in% expected_bg_methods)) {
    unexpected <- setdiff(file_bg_methods, expected_bg_methods)
    warning("Found files with unexpected background methods: ", 
            paste(unexpected, collapse=", "), 
            ". Expected: ", paste(expected_bg_methods, collapse=", "))
  }
  
  # Print summary of files to be processed
  cat("Found", length(env_files), "files to process:\n")
  cat(paste("-", basename(env_files)), sep = "\n")
  
  # Create list to store results
  all_results <- list()
  
  # Process each file
  for(i in seq_along(env_files)) {
    cat("\n\nProcessing file", i, "of", length(env_files), "\n")
    cat("File:", basename(env_files[i]), "\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    
    # Run single SDM training with force_retrain parameter
    results <- train_single_sdm(
      env_file = env_files[i],
      output_dir = output_dir,
      species_name = species_name,
      force_retrain = force_retrain
    )
    
    # Store results
    all_results[[basename(env_files[i])]] <- results
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  
  # Save combined results
  combined_file <- file.path(output_dir, species_name,
                             paste0(species_name, "_all_model_results.rds"))
  saveRDS(combined_results, combined_file)
  
  # Print summary
  cat("\n\nProcessing Complete\n")
  cat("===================\n")
  cat("Total files processed:", length(env_files), "\n")
  cat("Total models built/found:", nrow(combined_results), "\n")
  cat("Successful new models:", sum(combined_results$status == "success"), "\n")
  cat("Existing models skipped:", sum(combined_results$status == "existing"), "\n")
  cat("Failed models:", sum(!combined_results$status %in% c("success", "existing")), "\n")
  
  # Print breakdown by background method and threshold
  cat("\nResults by background method and threshold:\n")
  bg_thresh_counts <- table(combined_results$bg_method, combined_results$threshold)
  print(bg_thresh_counts)
  
  # Print breakdown by variable set
  cat("\nResults by variable set:\n")
  var_set_counts <- table(combined_results$var_set)
  print(var_set_counts)
  
  return(combined_results)
}

# Helper function to check if models already exist
check_existing_models <- function(species_out_dir, model_id) {
  # Check if model file exists
  model_file <- file.path(species_out_dir, paste0(model_id, ".sdm"))
  exists <- file.exists(model_file)
  
  # Return path if exists, FALSE otherwise
  if(exists) {
    return(model_file)
  } else {
    return(FALSE)
  }
}


#==============================================================================#
#                       2. Model building runs ----
#==============================================================================#

# standard runs - all algorithms
# # Run the model building function for Calophyllum
# calo_results <- train_all_sdm(
#   env_tables_dir = datafolder03,
#   output_dir = data__03_out,
#   species_name = "calo"
# )
# 
# gc()
# 
# # Run the model building function for Verticillium/Leptographium
# vertlept_results <- train_all_sdm(
#   env_tables_dir = datafolder03,
#   output_dir = data__03_out,
#   species_name = "vertlept"
# )

calo_results <- train_all_sdm(
  env_tables_dir = datafolder03,
  output_dir = data__03_out,
  species_name = "calo",
  force_retrain = FALSE  # Default is FALSE (only run for missing algs - i.e. GAM)
)
gc()
vertlept_results <- train_all_sdm(
  env_tables_dir = datafolder03,
  output_dir = data__03_out,
  species_name = "vertlept",
  force_retrain = FALSE  # Default is FALSE
)
gc()


# # # Clean up
# gc()
# rm(list = ls())

#==============================================================================#
#                        ----  End of workflow  ----
#==============================================================================#