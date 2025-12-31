# ---
# title: "02__variable_selection_by_species.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-24"
# update: "2025-12-31"
# notes:
# ---

## WHAT
# - Tests the 21 environmental predictor variables for collinearity between one another across two study extents,
# and multiple background point methods (calo: ecoregions, whole area; vertlept: 2far (pseudo), and whole area)
# - It outputs a .csv and .RData file of the names of the non-correlated, selected variables
# - calo can take forward 8 variables, vertlept can take forward 4 (forest and aspect removed from wilt)
# 
# HOW
# - corrr package is used to test for highly correlated variables (0.7, 0.5 thresholds)
# - VIF scores are computed, predictors are chosen based on a lower VIF score
# 
# WHY
# - To test which predictor variables are spatially autocorrelated / multicollinear
# - To test which variables are the most important at predicting species presence
# - Taking too many predictors into model building can lead to overfitting
# - Need robust models as the outputs of the SDMs become the inputs of the dynamic landscapes of RangeShiftR (large data chapter 6)
# 
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           ---- 0. Workspace set up ----
#==============================================================================#

## Setwd before running source() ----
# First set working directory so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

## Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, 
   align.forest.raster, align.raster, validate.processed.climate.files, calc_buffer_distance, calculate_metrics,
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, random_bg_whole_area)


## Set up base folders ----
# Base folders
datafolder <- file.path("//xxx")

# Data folder relating to script 00 (species occurrences)
#datafolder00 <- file.path(datafolder, "00__00__species_occurrence_cleaning")

# Data folder relating to script 01 (env data extraction)
#datafolder01 <- file.path(datafolder, "01__extract_species_enviro_data")

# Data folder relating to script 02 (collinearity)
datafolder02 <- file.path(datafolder, "02__variable_selection")

# Set output dir (for selected environmental variables)
out_dir__02 <- file.path(datafolder02, "output", "selected_vars")
dir.create(out_dir__02, recursive = TRUE, showWarnings = FALSE) 

# Create correlation plots directory
corr_dir <- file.path(out_dir__02, "corr_plots")
dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE) # sub folder for the correlation matrix plots

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
  forest = "Forest Cover",        # forest cover 2010
  aspect = "Aspect",              # aspect (derived from elevation, as a proxy for microclimate, measured in degrees)
  scenario = "Scenario",          # CHELSA shared socio-economic pathway (current, or future: ssp126, ssp370, ssp585)
  year_range = "Future year"      # future year time span (2011-2040, 2041-2070, 2071-2100)
)


#==============================================================================#
#                           ----  Helper functions  ----
#==============================================================================#

###  Get bg_set_name from file name ----
# Function to extract background method from filename
get_bg_method <- function(filename) {
  cat("\nChecking background method for file:", filename, "\n")
  # Extract method from filename
  if(grepl("ecor", filename, ignore.case = TRUE)) {
    cat("Found background method: ecor\n")
    return("ecor")
  } else if(grepl("whole", filename, ignore.case = TRUE)) {
    cat("Found background method: whole\n")
    return("whole")
  } else if(grepl("2far", filename, ignore.case = TRUE)) {
    cat("Found background method: 2far\n")
    return("2far")
  } else {
    warning("Unknown background method in filename: ", filename)
    return("unknown")
  }
}


### Variable selection for different species and thresholds ----
# Function to process species data with different thresholds


process_species_data <- function(extracted_file, species_name, max_vars,
                                 remove_forest = FALSE, remove_aspect = FALSE,
                                 thresholds = c(0.5, 0.7)) {
  
  cat("\nProcessing file:", extracted_file, "\n")
  cat("Species:", species_name, "\n")
  
  # Get background method from filename
  bg_method <- get_bg_method(extracted_file)
  cat("Background method identified:", bg_method, "\n")
  
  # Read and process data
  env_values <- fread(extracted_file)
  env_values_df <- as.data.frame(env_values)
  names(env_values_df) <- bio_lookup[names(env_values_df)]
  
  # Remove forest if specified
  if(remove_forest) {
    env_values_df <- subset(env_values_df, select = -c(`Forest Cover`))
  }
  
  # Remove aspect if specified
  if(remove_aspect) {
    env_values_df <- subset(env_values_df, select = -c(`Aspect`))
  }
  
  # Separate metadata and environmental columns
  meta_cols <- c("ID", "x", "y", "occ")
  env_cols <- setdiff(names(env_values_df), meta_cols)
  
  results_list <- list()
  
  # Process for each correlation threshold
  for(thresh in thresholds) {
    
    # Correlation analysis
    cor_matrix <- cor(env_values_df[, env_cols])
    high_cor <- findCorrelation(cor_matrix, cutoff = thresh)
    selected_vars <- env_cols[-high_cor]
    
    # Calculate VIF scores for remaining variables
    env_subset <- env_values_df[, selected_vars]
    vif_scores <- usdm::vif(env_subset)
    
    # Order by VIF and select top variables based on species limit
    vif_ordered <- vif_scores[order(vif_scores$VIF), ]
    final_vars <- head(vif_ordered$Variables, max_vars)
    
    # Create final datasets
    final_data <- cbind(
      env_values_df[, meta_cols],
      env_values_df[, final_vars]
    )
    
    # Store results
    results_list[[paste0("thresh_", thresh)]] <- list(
      data = final_data,
      selected_vars = final_vars,
      cor_matrix = cor_matrix,
      vif_scores = vif_ordered
    )
    
    # Generate plots with bg_set identifier
    # 1. Correlation plot
    png(file.path(corr_dir, 
                  paste0(species_name, "_", bg_method, "_corrplot_", thresh, ".png")),
        width = 1200, height = 1200, res = 150)
    corrplot(cor_matrix, 
             method = "color", 
             type = "upper",
             tl.col = "black", 
             tl.srt = 45,
             number.cex = 0.5,        # Reduce text size inside squares
             tl.cex = 0.7,           # Reduce variable names text size
             addCoef.col = "black",  
             digits = 1)             # Reduce to 1 decimal place
    dev.off()
    
    # 2. VIF scores plot
    ggplot(vif_ordered, aes(x = reorder(Variables, -VIF), y = VIF)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(species_name, bg_method, "VIF Scores (corr threshold:", thresh, ")"),
           x = "Variables", y = "VIF Score")
    ggsave(file.path(corr_dir, 
                     paste0(species_name, "_", bg_method, "_vif_plot_", thresh, ".png")))
  }
  
  return(results_list)
}


### Save outputs to file ----
# Save outputs for each species and threshold
# Includes bg_set_name to name according to how background points were generated
save_species_results <- function(results, species_name, bg_method, out_dir) {
  cat("\nSaving results for:", species_name, "\n")
  cat("Background method:", bg_method, "\n")
  
  for(thresh in names(results)) {
    cat("\nProcessing threshold:", thresh, "\n")
    
    # Create filenames
    base_filename <- paste0(species_name, "_", bg_method, "_", thresh)
    cat("Creating files with base name:", base_filename, "\n")
    
    # Save full dataset with clear naming
    rds_file <- file.path(out_dir, paste0(base_filename, "_dataset.RDS"))
    cat("Saving dataset to:", rds_file, "\n")
    saveRDS(results[[thresh]]$data, rds_file)
    
    # Save CSV version
    csv_file <- file.path(out_dir, paste0(base_filename, "_dataset.csv"))
    cat("Saving CSV to:", csv_file, "\n")
    write.csv(results[[thresh]]$data, csv_file, row.names = FALSE)
    
    # Save variable names vector
    vars_file <- file.path(out_dir, paste0(base_filename, "_predsel_only.RDS"))
    cat("Saving variables to:", vars_file, "\n")
    saveRDS(results[[thresh]]$selected_vars, vars_file)
  }
}


#==============================================================================#
#                           1. Main processing  ----
#==============================================================================#


### Main processing function (all combined)
# Main processing function for all background sets of a species
process_all_species_sets <- function(extracted_files, species_name, max_vars, remove_forest = FALSE,
                                     remove_aspect = FALSE) {
  for(file in extracted_files) {
    # Extract background method from filename
    bg_method <- get_bg_method(file)
    
    results <- process_species_data(
      extracted_file = file,
      species_name = species_name,
      max_vars = max_vars,
      remove_forest = remove_forest,
      remove_aspect = remove_aspect
    )
    
    save_species_results(results, species_name, bg_method, out_dir__02)
  }
}

# List the combined species-enviro dfs for each species
calo_files <- list.files(file.path(datafolder02, "input", "extracted_envdata", "current"), pattern="calo", full.names = T)
vertlept_files <- list.files(file.path(datafolder02, "input", "extracted_envdata", "current"), pattern="vertlept", full.names = T)


##### Process Calophyllum environmental predictor variables ----
# Process Calo (max 8 variables, base is without forest and aspect)
process_all_species_sets(
  extracted_files = calo_files,
  species_name = "calo",
  max_vars = 8,
  remove_forest = TRUE,
  remove_aspect = TRUE
)

gc()

##### Process Verticillium/Leptographium environmental predictor variables ----
# Process vertlept (max 4 variables, base is without forest and aspect)
process_all_species_sets(
  extracted_files = vertlept_files,
  species_name = "vertlept",
  max_vars = 4,
  remove_forest = TRUE,
  remove_aspect = TRUE
)


#==============================================================================#
#                       ----  2. Selection analysis  ----
#==============================================================================#

### Per species analysis ----

# Function per species summary/analysis
analyze_variables_per_spp <- function(directory_path) {
  # Get list of files
  files <- list.files(path = directory_path, 
                      pattern = "_predsel_only.RDS$", 
                      full.names = TRUE,
                      recursive = FALSE)
  
  if(length(files) == 0) {
    stop("No valid files found for analysis")
  }
  
  # Create results table
  results <- data.frame(
    species = character(),
    bg_method = character(),
    threshold = character(),
    n_variables = numeric(),
    variables = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each file
  for(file in files) {
    filename <- basename(file)
    cat("Processing file:", filename, "\n")
    
    # Parse filename components
      parts <- strsplit(filename, "_")[[1]]
      species <- parts[1]
      bg_method <- parts[2]
      threshold <- paste0(parts[3], "_", parts[4])
  
    
    # Load selected variables
    selected_vars <- readRDS(file)
    
    # Add to results
    results <- rbind(results, data.frame(
      species = species,
      bg_method = bg_method,
      threshold = threshold,
      n_variables = length(selected_vars),
      variables = paste(selected_vars, collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort results
  results <- results[order(results$species, results$bg_method, results$threshold), ]
  
  # Initialize list to store species-specific results
  species_results <- list()
  
  # Process each species separately
  for(sp in unique(results$species)) {
    sp_data <- results[results$species == sp, ]
    
    # 1. Create comparison plot for this species
    var_comp_plot <- ggplot(sp_data, 
                            aes(x = threshold, y = n_variables, fill = bg_method)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = sprintf("Variable Selection Comparison for %s", sp),
           x = "Threshold",
           y = "Number of Variables",
           fill = "Background Method") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 2. Create variable occurrence matrix
    vars_list <- strsplit(sp_data$variables, ", ")
    unique_vars <- unique(unlist(vars_list))
    var_matrix <- matrix(0, 
                         nrow = length(unique_vars),
                         ncol = nrow(sp_data),
                         dimnames = list(unique_vars,
                                         paste(sp_data$bg_method, sp_data$threshold)))
    
    for(i in 1:nrow(sp_data)) {
      these_vars <- unlist(strsplit(sp_data$variables[i], ", "))
      var_matrix[these_vars, i] <- 1
    }
    
    # Convert to data frame for plotting
    var_matrix_df <- as.data.frame(var_matrix) %>%
      tibble::rownames_to_column("variable")
    
    var_matrix_long <- tidyr::pivot_longer(var_matrix_df,
                                           cols = -variable,
                                           names_to = "combination",
                                           values_to = "present")
    
    # 3. Create heatmap of variable presence
    var_heatmap <- ggplot(var_matrix_long,
                          aes(x = combination, y = variable, fill = factor(present))) +
      geom_tile() +
      scale_fill_manual(values = c("0" = "white", "1" = "steelblue"),
                        labels = c("Absent", "Selected")) +
      theme_minimal() +
      labs(title = sprintf("Variable Selection Pattern for %s", sp),
           x = "Method and Threshold Combination",
           y = "Environmental Variable",
           fill = "Status") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save species-specific plots
    tryCatch({
      ggsave(file.path(directory_path, 
                       sprintf("%s_variable_comparison.png", sp)), 
             var_comp_plot,
             width = 8, height = 6)
      ggsave(file.path(directory_path, 
                       sprintf("%s_variable_heatmap.png", sp)), 
             var_heatmap,
             width = 10, height = 8)
    }, error = function(e) {
      warning(sprintf("Error saving plots for %s: %s", sp, e$message))
    })
    
    # Create species-specific summary text file
    sink(file.path(directory_path, sprintf("%s_variable_summary.txt", sp)))
    cat(sprintf("Variable Selection Analysis for %s\n", sp))
    cat("=====================================\n\n")
    
    # Compare variables across methods and thresholds
    for(th in unique(sp_data$threshold)) {
      cat(sprintf("\nThreshold: %s\n", th))
      th_data <- sp_data[sp_data$threshold == th, ]
      
      cat("\nVariables by background method:\n")
      for(bg in unique(sp_data$bg_method)) {
        bg_vars <- th_data$variables[th_data$bg_method == bg]
        if(length(bg_vars) > 0) {
          cat(sprintf("  %s: %s\n", bg, bg_vars))
        }
      }
      
      # If we have multiple methods, show common and unique variables
      bg_methods <- unique(sp_data$bg_method)
      if(length(bg_methods) > 1) {
        vars_by_method <- lapply(bg_methods, function(bg) {
          vars <- th_data$variables[th_data$bg_method == bg]
          if(length(vars) > 0) {
            unlist(strsplit(vars, ", "))
          } else {
            character(0)
          }
        })
        
        if(length(vars_by_method) >= 2 && all(sapply(vars_by_method, length) > 0)) {
          common_vars <- Reduce(intersect, vars_by_method)
          cat("\nCommon variables across methods:\n  ", 
              paste(common_vars, collapse = ", "), "\n")
          
          # Unique variables for each method
          for(i in seq_along(bg_methods)) {
            unique_vars <- setdiff(vars_by_method[[i]], 
                                   unlist(vars_by_method[-i]))
            if(length(unique_vars) > 0) {
              cat(sprintf("\nUnique to %s:\n  %s\n", 
                          bg_methods[i], 
                          paste(unique_vars, collapse = ", ")))
            }
          }
        }
      }
    }
    sink()
    
    # Store species-specific results
    species_results[[sp]] <- list(
      data = sp_data,
      comparison_plot = var_comp_plot,
      heatmap = var_heatmap,
      variable_matrix = var_matrix
    )
  }
  
  # Print summary of what was processed
  cat(sprintf("\nProcessed %d files\n", nrow(results)))
  cat(sprintf("Found %d unique species\n", length(unique(results$species))))
  cat(sprintf("Found %d unique background methods\n", length(unique(results$bg_method))))
  
  return(list(
    results_table = results,
    species_results = species_results
  ))
}

# Run the second analysis function per spp
analysis_results_per_spp <- analyze_variables_per_spp(out_dir__02)

# Access results
analysis_results_per_spp$results_table  # Overall summary table
analysis_results_per_spp$species_results  # List containing species-specific analyses


#==============================================================================#
#                ---- 2a. Calo specific: manual pred vars  ----
#==============================================================================#


# Create manual predictor set for Calo based on expert knowledge
create_manual_predictor_set <- function(extracted_files, out_dir) {
  # Define expert-selected variables
  expert_vars <- c(
    "bio3",   # Isothermality
    "bio4",   # Temperature seasonality
    "bio5",   # Max temp warmest month
    "bio6",   # Min temp coldest month
    "bio13",  # Precip wettest month
    "bio14",  # Precip driest month
    "bio15"   # Precip seasonality
  )
  
  # Required metadata columns
  required_cols <- c("x", "y", "occ", "ID")
  
  # Process each background method file
  for(file in extracted_files) {
    # Get background method from filename
    bg_method <- get_bg_method(file)
    cat(sprintf("\nProcessing %s background method from file: %s\n", bg_method, basename(file)))
    
    # Read environmental data
    env_values <- fread(file)
    
    # Verify all required columns exist
    missing_required <- setdiff(required_cols, names(env_values))
    if(length(missing_required) > 0) {
      stop(sprintf("Missing required columns: %s", 
                   paste(missing_required, collapse = ", ")))
    }
    
    missing_expert <- setdiff(expert_vars, names(env_values))
    if(length(missing_expert) > 0) {
      stop(sprintf("Missing expert variables: %s", 
                   paste(missing_expert, collapse = ", ")))
    }
    
    # Create variants of variable sets
    variants <- list(
      base = list(
        vars = expert_vars,
        suffix = "expert"
      ),
      forest = list(
        vars = c(expert_vars, "forest"),
        suffix = "expert_forest"
      ),
      aspect = list(
        vars = c(expert_vars, "aspect"),
        suffix = "expert_aspect"
      ),
      forest_aspect = list(
        vars = c(expert_vars, "forest", "aspect"),
        suffix = "expert_forest_aspect"
      )
    )
    
    # Process each variant
    for(variant_name in names(variants)) {
      variant <- variants[[variant_name]]
      current_vars <- variant$vars
      
      # Create filenames matching original pattern
      base_filename <- paste0("calo_", bg_method, "_", variant$suffix)
      
      # Create and save the complete dataset for this variant
      if(all(current_vars %in% names(env_values))) {
        final_data <- data.frame(
          env_values[, ..required_cols],  # metadata
          env_values[, ..current_vars]    # selected variables
        )
        
        # Save the variable list (matching original _predsel_only.RDS pattern)
        saveRDS(current_vars, 
                file.path(out_dir, paste0(base_filename, "_predsel_only.RDS")))
        
        # Save complete dataset
        saveRDS(final_data,
                file.path(out_dir, paste0(base_filename, "_dataset.RDS")))
        write.csv(final_data,
                  file.path(out_dir, paste0(base_filename, "_dataset.csv")),
                  row.names = FALSE)
        
        # Create SDM formula
        sdm_formula <- paste("presence ~", paste(current_vars, collapse = " + "))
        writeLines(sdm_formula, 
                   file.path(out_dir, paste0(base_filename, "_sdm_formula.txt")))
      } else {
        warning(sprintf("Could not find all variables for %s variant", variant_name))
      }
    }
    
    # Run correlation analysis for base variables
    cor_matrix <- cor(env_values[, ..expert_vars])
    
    # Save correlation matrix in corr_dir (assuming it exists from original code)
    write.csv(cor_matrix,
              file.path(corr_dir, paste0("calo_", bg_method, "_expert_corrplot.csv")))
    
    # Create correlation plot
    png(file.path(corr_dir, paste0("calo_", bg_method, "_expert_corrplot.png")),
        width = 1200, height = 1000, res = 150)
    corrplot(cor_matrix,
             method = "color",
             type = "upper",
             tl.col = "black",
             tl.srt = 45,
             addCoef.col = "black",
             col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
             title = paste("Correlation Matrix of Expert-Selected Variables -", bg_method),
             mar = c(0,0,1,0))
    dev.off()
  }
  
  # Create documentation file
  sink(file.path(out_dir, "expert_selection_documentation.txt"))
  cat("Expert-Selected Variables for Calophyllum\n")
  cat("=======================================\n\n")
  
  cat("Selected variables:\n")
  for(var in expert_vars) {
    cat(sprintf("%s - %s\n", var, bio_lookup[var]))
  }
  
  cat("\nBackground methods processed:\n")
  bg_methods <- sapply(extracted_files, get_bg_method)
  cat(paste(unique(bg_methods), collapse = "\n"))
  
  cat("\n\nFiles created for each background method:\n")
  for(variant_name in names(variants)) {
    variant <- variants[[variant_name]]
    cat(sprintf("\n%s variant:\n", variant_name))
    cat("Variables included:", paste(variant$vars, collapse = ", "), "\n")
  }
  sink()
  
  return(list(
    expert_vars = expert_vars,
    variants = variants,
    background_methods = unique(bg_methods)
  ))
}

# Run the expert selected variables (temp and precip seasonality and extremes) for Calophyllum
expert_results <- create_manual_predictor_set(
  extracted_files = calo_files[grep("(ecor|whole)", calo_files)],
  out_dir = out_dir__02
)


#==============================================================================#
#                  ----  3. Create final subsets for models  ----
#==============================================================================#


# Function to create expanded variable sets and format for SDM package
create_sdm_variable_sets <- function(analysis_results, out_dir) {
  
  # Define expected selection types
  selection_types <- c("thresh_0.5", "thresh_0.7", "expert")
  
  # Create directory for SDM variable sets if it doesn't exist
  sdm_vars_dir <- file.path(out_dir, "sdm_variable_sets")
  dir.create(sdm_vars_dir, showWarnings = FALSE)
  
  # Get all relevant RDS files
  all_files <- list.files(out_dir, pattern = "_predsel_only.RDS$", full.names = TRUE)
  
  # Process each species
  for(sp in unique(c(analysis_results$results_table$species, "calo"))) {
    cat(sprintf("\nProcessing species: %s\n", sp))
    cat("===============================\n")
    
    # Create directory for species
    sp_dir <- file.path(sdm_vars_dir, sp)
    dir.create(sp_dir, showWarnings = FALSE)
    
    # Get files for this species
    sp_files <- grep(paste0(sp, "_"), all_files, value = TRUE)
    
    # Track what we find for summary
    found_combinations <- list()
    
    # Process each file
    for(file in sp_files) {
      filename <- basename(file)
      parts <- strsplit(filename, "_")[[1]]
      bg_method <- parts[2]
      
      # Determine selection type
      if(grepl("expert", filename)) {
        selection <- "expert"
      } else if(grepl("0.5", filename)) {
        selection <- "thresh_0.5"
      } else if(grepl("0.7", filename)) {
        selection <- "thresh_0.7"
      } else {
        warning(sprintf("Unknown selection type in file: %s", filename))
        next
      }
      
      cat(sprintf("\nProcessing %s - Background: %s, Selection: %s\n", 
                  sp, bg_method, selection))
      
      # Get base variables
      base_vars <- readRDS(file)
      if(!is.character(base_vars)) {
        warning(sprintf("Unexpected format in file %s", filename))
        next
      }
      
      # Create the four variants
      var_sets <- list(
        base = base_vars,
        forest = c(base_vars, "forest"),
        aspect = c(base_vars, "aspect"),
        forest_aspect = c(base_vars, "forest", "aspect")
      )
      
      # Save each variant
      for(set_name in names(var_sets)) {
        # Create descriptive filename
        variant_filename <- paste0(sp, "_", bg_method, "_", selection, "_", set_name, ".rds")
        
        cat(sprintf("  Creating variant: %s\n", set_name))
        
        # Save as RDS
        saveRDS(var_sets[[set_name]], 
                file.path(sp_dir, variant_filename))
        
        # Create SDM-formatted variable string
        sdm_vars <- paste(var_sets[[set_name]], collapse = " + ")
        
        # Save formula for SDM package
        formula_file <- sub(".rds", "_sdm_formula.txt", variant_filename)
        sdm_formula <- paste("presence ~", sdm_vars)
        writeLines(sdm_formula, 
                   file.path(sp_dir, formula_file))
        
        cat(sprintf("    - Saved: %s\n", variant_filename))
        cat(sprintf("    - Created formula: %s\n", formula_file))
      }
      
      # Track what we found
      combo_key <- paste(bg_method, selection)
      found_combinations[[combo_key]] <- list(
        bg_method = bg_method,
        selection = selection,
        variables = base_vars
      )
    }
    
    # Create summary file for this species
    sink(file.path(sp_dir, paste0(sp, "_variable_sets_summary.txt")))
    cat(sprintf("Variable Sets Summary for %s\n", sp))
    cat("==============================\n\n")
    
    # Summarize what was found vs expected
    cat("Selection Methods Found:\n")
    cat("------------------------\n")
    bg_methods <- unique(sapply(found_combinations, function(x) x$bg_method))
    
    for(bg in bg_methods) {
      cat(sprintf("\nBackground Method: %s\n", bg))
      for(sel in selection_types) {
        combo_key <- paste(bg, sel)
        if(combo_key %in% names(found_combinations)) {
          cat(sprintf("  %s: Found (%d variables)\n", 
                      sel, 
                      length(found_combinations[[combo_key]]$variables)))
        } else {
          cat(sprintf("  %s: Not found\n", sel))
        }
      }
    }
    
    cat("\nDetailed Variable Sets:\n")
    cat("---------------------\n")
    
    # Process each combination found
    for(combo in names(found_combinations)) {
      info <- found_combinations[[combo]]
      
      cat(sprintf("\nCombination: %s (%s)\n", combo, info$selection))
      cat("------------------\n")
      
      base_vars <- info$variables
      
      cat("\n1. Base variables:\n")
      cat(paste(base_vars, collapse = ", "), "\n")
      
      cat("\n2. Base + Forest:\n")
      cat(paste(c(base_vars, "forest"), collapse = ", "), "\n")
      
      cat("\n3. Base + Aspect:\n")
      cat(paste(c(base_vars, "aspect"), collapse = ", "), "\n")
      
      cat("\n4. Base + Forest + Aspect:\n")
      cat(paste(c(base_vars, "forest", "aspect"), collapse = ", "), "\n\n")
    }
    sink()
    
    cat("\nCompleted processing for", sp, "\n")
    cat("----------------------------------------\n")
  }
  
  return(list(
    directory = sdm_vars_dir,
    found_combinations = found_combinations
  ))
}

# Updated load function to handle all selection types
load_sdm_variables <- function(species, bg_method, selection, variant, sdm_vars_dir) {
  # Validate selection type
  valid_selections <- c("thresh_0.5", "thresh_0.7", "expert")
  if(!selection %in% valid_selections) {
    stop(sprintf("Invalid selection type. Must be one of: %s", 
                 paste(valid_selections, collapse = ", ")))
  }
  
  # Construct file path
  sp_dir <- file.path(sdm_vars_dir, species)
  filename <- paste0(species, "_", bg_method, "_", selection, "_", variant, ".rds")
  file_path <- file.path(sp_dir, filename)
  
  # Load and return variables
  if(file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    stop(sprintf("Variable set file not found: %s", filename))
  }
}

# 1. Create all variable sets
sdm_vars_dir <- create_sdm_variable_sets(analysis_results_per_spp, out_dir__02)


# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                        ----  End of workflow  ----
#==============================================================================#