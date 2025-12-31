# ---
# title: "05__ensemble_projection_feasibility.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-28"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Runs feasibility checks on top 3 ensemble sdms per spp - projection on whole MDG study area (template A)
# - Creates raster habitat suitability maps using current climate
# 
# HOW
# - Runs code locally sequentially due to memory and cluster errors
# 
# WHY
# - Will only take the most feasible model foward to make full ensemble future predictions per species
# 
# WHERE
# - Latest scripts always saved to Github, with copies saved locally and on servers if required
# - Data for code always taken from server (xxx) and run from there or saved locally to machine TEMP folder but should be deleted
# once output computed/backed up


#==============================================================================#
#                           0. Workspace set up ----
#==============================================================================#


## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")



##### Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, 
   align.forest.raster, align.raster, validate.processed.climate.files, calc_buffer_distance, calculate_metrics,
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, 
   random_bg_whole_area, rename_stack_layers, select07_cv, 
   analyze_elev_distribution, generate_buf_bg_pts, generate_dispersal_bg_pts, generate_elev_bg_pts, plot_bg_comparison)


##### CRS: African Equal Area ----
crs <- "ESRI:102022"


##### Base folders for existing data ----
datafolder <- file.path("//xxx")


# Load required libraries
library(terra)
library(data.table)
library(sdm)
installAll()
# Additional packages required for SDM models to work properly
library(dismo)
library(raster)
library(sp)
library(gbm)
library(mgcv)
library(randomForest)
library(rpart)
library(kernlab)

# Input data folders
datafolder03 <- file.path(datafolder, "03__model_training", "output")
datafolder05_in <- file.path(datafolder, "05__model_prediction", "input", "variable_stacks", "MDG", "top1")

# Output data folder for this script
datafolder05_out <- file.path(datafolder, "05__model_prediction", "output", "current_projections", "weighted_mean_ensembles", "extra_settings")
dir.create(file.path(datafolder05_out), recursive = TRUE, showWarnings = FALSE)


#==============================================================================#
#                 ---- 1. Find & load enviro data stacks ----
#==============================================================================#

# Locate and load stacks relating to the top 3 models per spp

# These rasters directly related to the top 3 models produced in 04__model_perf_eval.R output:
# "top_models_main_table.csv"
# found in: \04__model_eval_perf\output\top_performing_models\publication_tables

#### Calo load 'top 3' raster stacks ----

# 1 : calo_whole_thresh_0.7_forest_alg3.sdm
# 2 : calo_whole_thresh_0.7_base_alg3.sdm
# 3 : calo_whole_thresh_0.7_aspect_alg3.sdm

# 1 : current env_data_stack for "calo_whole_thresh_0.7_forest_alg3.sdm"
calo_whole_thresh_0.7_forest_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                      pattern="calo_whole_thresh_0.7_forest_current_current__MDG.tif", full.names = T))

# 2 : current env_data_stack for "calo_whole_thresh_0.7_base_alg3.sdm"
calo_whole_thresh_0.7_base_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                      pattern="calo_whole_thresh_0.7_base_current_current__MDG.tif", full.names = T))

# 3 : current env_data_stack for "calo_whole_thresh_0.7_aspect_alg3.sdm"
calo_whole_thresh_0.7_forest_aspect_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                      pattern="calo_whole_thresh_0.7_aspect_current_current__MDG.tif", full.names = T))


#### vertlept load 'top 3' raster stacks ----

# 1 : vertlept_2far_thresh_0.7_base_alg3.sdm
# 2 : vertlept_2far_thresh_0.7_aspect_alg3.sdm
# 3 : vertlept_2far_thresh_0.7_forest_alg3.sdm

# 1 : current env_data_stack for "vertlept_2far_thresh_0.7_base_alg3.sdm"
vertlept_2far_thresh_0.7_base_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                        pattern="vertlept_2far_thresh_0.7_base_current_current__MDG.tif", full.names = T))

# 2 : current env_data_stack for "vertlept_2far_thresh_0.7_aspect_alg3.sdm"
vertlept_2far_thresh_0.7_aspect_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                        pattern="vertlept_2far_thresh_0.7_aspect_current_current__MDG.tif", full.names = T))

# 3:  current env_data_stack for "vertlept_2far_thresh_0.7_forest_alg3.sdm" 
vertlept_2far_thresh_0.7_forest_curr_stack <- terra::rast(list.files(file.path(datafolder05_in),
                                        pattern="vertlept_2far_thresh_0.7_forest_current_current__MDG.tif", full.names = T))


#==============================================================================#
#                 ---- 2. Find and load top 3 models per spp ----
#==============================================================================#

#### Calo load top 3 models ----
# 1 : calo_whole_thresh_0.7_forest_alg3.sdm
# 2 : calo_whole_thresh_0.7_base_alg3.sdm
# 3 : calo_whole_thresh_0.7_aspect_alg3.sdm

load(file.path(datafolder03, "calo", "calo_whole_thresh_0.7_forest_alg3.sdm"))
calo_whole_thresh_0.7_forest_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)

load(file.path(datafolder03, "calo", "calo_whole_thresh_0.7_base_alg3.sdm"))
calo_whole_thresh_0.7_base_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)

load(file.path(datafolder03, "calo", "calo_whole_thresh_0.7_aspect_alg3.sdm"))
calo_whole_thresh_0.7_forest_aspect_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)


#### vertlept load top 3 models ----
# 1 : vertlept_2far_thresh_0.7_base_alg3.sdm
# 2 : vertlept_2far_thresh_0.7_aspect_alg3.sdm
# 3 : vertlept_2far_thresh_0.7_forest_alg3.sdm

load(file.path(datafolder03, "vertlept", "vertlept_2far_thresh_0.7_base_alg3.sdm"))
vertlept_2far_thresh_0.7_base_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)

load(file.path(datafolder03, "vertlept", "vertlept_2far_thresh_0.7_aspect_alg3.sdm"))
vertlept_2far_thresh_0.7_aspect_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)

load(file.path(datafolder03, "vertlept", "vertlept_2far_thresh_0.7_forest_alg3.sdm"))
vertlept_2far_thresh_0.7_forest_alg3 <- sdm_model # rename to match model base/file name

rm(sdm_model)


#==============================================================================#
#                 ---- 3. Make continuous current projections ----
#==============================================================================#

# Weighted ensemble projections across current climate of MDG - top 3 models
# Weighting chosen to optimise AUC performance, ensemble optimised using max(se+sp) - max TSS
# Weighting with the performance calculation to rank the top 3 models in "04__top_perf_models.R" script

# 1 : calo_whole_thresh_0.7_forest_alg3.sdm
# 2 : calo_whole_thresh_0.7_base_alg3.sdm
# 3 : calo_whole_thresh_0.7_aspect_alg3.sdm

#### Calo ensemble projections ----
### 1 : calo_whole_thresh_0.7_forest_alg3.sdm
#### Ensemble weighted AUC with opt 2 = "max(se+sp)" (maximum sum of sensitivity and specificity, which equals maximizing TSS).
ensname_proj_calo1_AUCopt2 <- file.path(datafolder05_out, "ensproj_AUCopt2_calo_whole_thresh_0.7_forest_alg3.tif")
ens_projcalo1_AUCopt2 <- ensemble(calo_whole_thresh_0.7_forest_alg3, calo_whole_thresh_0.7_forest_curr_stack, ensname_proj_calo1_AUCopt2,
                          setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                       stat='AUC',opt=2, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
gc()
#### 2 : calo_whole_thresh_0.7_base_alg3.sdm
ensname_proj_calo2_AUCopt2 <- file.path(datafolder05_out, "ensproj_AUCopt2_calo_whole_thresh_0.7_base_alg3.tif")
ens_projcalo2_AUCopt2 <- ensemble(calo_whole_thresh_0.7_base_alg3, calo_whole_thresh_0.7_base_curr_stack, ensname_proj_calo2_AUCopt2,
                          setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                       stat='AUC',opt=2, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
gc()
#### 3 : calo_whole_thresh_0.7_aspect_alg3.sdm
ensname_proj_calo3_AUCopt2 <- file.path(datafolder05_out,  "ensproj_AUCopt2_calo_whole_thresh_0.7_aspect_alg3.tif")
ens_projcalo3_AUCopt2 <- ensemble(calo_whole_thresh_0.7_aspect_alg3, calo_whole_thresh_0.7_aspect_curr_stack,
                          ensname_proj_calo3_AUCopt2,
                          setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                       stat='AUC',opt=2, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
gc()

#### vertlept ensemble projections ----
# Weighting chosen to optimise MCC (Matthews Correlation Coefficient) performance, 
# ensemble optimised using "minROCdist" (opt=4)

# 1 : vertlept_2far_thresh_0.7_base_alg3.sdm
# 2 : vertlept_2far_thresh_0.7_aspect_alg3.sdm
# 3 : vertlept_2far_thresh_0.7_forest_alg3.sdm

#### 1 : vertlept_2far_thresh_0.7_base_alg3.sdm
##### Ensemble projections
ensname_proj_vertlept1_MCCopt4 <- file.path(datafolder05_out, "ensproj_MCCopt4_vertlept_2far_thresh_0.7_base_alg3.tif")
ensproj_vertlept1_MCCopt4 <- ensemble(vertlept_2far_thresh_0.7_base_alg3, vertlept_2far_thresh_0.7_base_curr_stack, 
                              ensname_proj_vertlept1_MCCopt4, 
                              setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                           stat='MCC',opt=4, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
writeRaster(ensname_proj_vertlept1_MCCopt4, )
gc()
#### 2 : vertlept_2far_thresh_0.7_aspect_alg3.sdm
ensname_proj_vertlept2_MCCopt4  <- file.path(datafolder05_out, "ensproj_MCCopt4_vertlept_2far_thresh_0.7_aspect_alg3.tif")
ensproj_vertlept2_MCCopt4 <- ensemble(vertlept_2far_thresh_0.7_aspect_alg3, vertlept_2far_thresh_0.7_aspect_curr_stack, 
                              ensname_proj_vertlept2_MCCopt4, 
                              setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                           stat='MCC',opt=4, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
gc()

#### 3 : vertlept_2far_thresh_0.7_forest_alg3.sdm
ensname_proj_vertlept3_MCCopt4  <- file.path(datafolder05_out, "ensproj_MCCopt4_vertlept_2far_thresh_0.7_forest_alg3.tif")
ensproj_vertlept3_MCCopt4 <- ensemble(vertlept_2far_thresh_0.7_forest_alg3, vertlept_2far_thresh_0.7_forest_curr_stack, 
                              ensname_proj_vertlept3_MCCopt4, 
                              setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                           stat='MCC',opt=4, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)


# Write a copy of the final models to file to easily findable in 05__model_prediction.R script
# calo
#write.sdm(ensname_proj_calo1_AUCopt2, filename = file.path(datafolder03, "calo_top1_AUCopt2.sdm"))
# vert
#write.sdm(ensproj_vertlept1_MCCopt4, filename = file.path(datafolder03, "vertlept_top1_MCCopt4.sdm"))


#==============================================================================#
#             ---- 4. Make binary thresholds from trained models ----
#==============================================================================#

# opt can be specified to select one of the criteria for optimising the threshold.
# The value can be between 1 to 15 for "sp=se", "max(se+sp)", "min(cost)",
# "minROCdist", "max(kappa)", "max(ppv+npv)", "ppv=npv", "max(NMI)", "max(ccr)",
# "prevalence", "max(MCC)", "P10", "P5", "P1", "P0" criteria.
# P10, P5, P1 refer to 10, 5, and 1 percentile of presence records in the evaluation dataset,
# for which the suitability value is used as the threshold.
# By choosing P0, the minimum suitability value across presence records is selected as the threshold.


# 1 : calo_whole_thresh_0.7_forest_alg3.sdm

# 1 : vertlept_2far_thresh_0.7_base_alg3.sdm

# Load each model gui manually to check thresholding value
sdm::gui(calo_whole_thresh_0.7_forest_alg3) # threshold opt 2 = max(se+sp) (maximizing the sum of sensitivity and specificity)
sdm::gui(vertlept_2far_thresh_0.7_base_alg3) # threshold opt 2 = ?

## Get individual thresholds for chosen models ----
# (top 1 model from the initial ranking based on different optimisations of ensemble performance metrics during thresholding)
#### Calo top threshold using opt2 "max(se+sp)" ----
# Jiménez-Valverde & Lobo (2007) showed works well for native range models
thresh_calo_opt2 <- sdm::threshold(calo_whole_thresh_0.7_forest_alg3,id="ensemble",opt=2)
gc()

#### Vertlept top threshold using opt 4 "minROCdist" (minimum distance to ROC curve) ----
# Liu et al. (2016) showed this works better for transferability across regions
thresh_vertlept_opt4 <- sdm::threshold(vertlept_2far_thresh_0.7_base_alg3,id="ensemble",opt=4)
gc()

# Secondary thresholding to account for small occurrence data
### Thresholds using opt 12 - "10th Percentile Training Presence Thresholds / P10" ----
# opt = P0, the minimum suitability value across presence records
thresh_calo_opt12 <- sdm::threshold(calo_whole_thresh_0.7_forest_alg3,id="ensemble",opt=12)
gc()
thresh_vertlept_opt12 <- sdm::threshold(vertlept_2far_thresh_0.7_base_alg3,id="ensemble",opt=12)
gc()

# Make intermediate output folder
thres_outdir <- file.path(datafolder, "05__model_prediction", "intermediate", "thresholds")
dir.create(file.path(thres_outdir), recursive = TRUE, showWarnings = FALSE)

### Write threshold values to file for binary predictions ----

# Save top thresholds opt 2 & 4 to file to use to make binary predictions for future climate scenarios
saveRDS(thresh_calo_opt2, file=file.path(thres_outdir, "thresh_calo_opt2.rds"))
saveRDS(thresh_vertlept_opt4, file=file.path(thres_outdir, "thresh_vertlept_opt4.rds"))

# Save secondary thresholds opt 12 to file which account for small species occurrences
saveRDS(thresh_calo_opt12, file=file.path(thres_outdir, "thresh_calo_opt12.rds"))
saveRDS(thresh_vertlept_opt12, file=file.path(thres_outdir, "thresh_vertlept_opt12.rds"))


#==============================================================================#
#                     ---- 5. Calculate suitable cells per thresh ----
#==============================================================================#

#  N.B. Need the models to be loaded for this section:
# Extract values from prediction at test locations
# Extract coordinates from the sdm object

  # Function to count suitable pixels
  count_suitable_cells <- function(rast_layer, threshold) {
    # Create binary raster where TRUE = suitable
    binary_rast <- rast_layer >= threshold
    
    # Get frequency table
    freq_table <- freq(binary_rast)
    
    # Find and sum the TRUE/1 cells
    suitable_count <- 0
    true_row <- which(freq_table$value == 1 | freq_table$value == TRUE)
    if(length(true_row) > 0) {
      suitable_count <- freq_table$count[true_row]
    }
    
    return(suitable_count)
  }
  
  # Function to calculate classification accuracy
  calculate_accuracy <- function(rast_file, threshold, sdm_obj, layer_index = 1) {
    # Load raster if it's a file path
    if(is.character(rast_file)) {
      rast <- rast(rast_file)
    } else {
      rast <- rast_file
    }
    
    # Select the desired layer if multi-layer raster
    if(nlyr(rast) > 1) {
      rast <- rast[[layer_index]]
    }
    
    # Extract presence/absence data from SDM object
    species_data <- sdm_obj@data@species$occ
    
    # Create a numeric vector of 0s (absence)
    pa_vector <- numeric(length(sdm_obj@data@info@coords[,1]))
    
    # Set presence indices to 1
    pa_vector[species_data@presence] <- 1
    
    # Get coordinates from SDM object
    coords <- sdm_obj@data@info@coords[, 2:3]
    
    # Extract values at those coordinates
    extracted <- terra::extract(rast, coords)
    
    # Get the prediction values
    if(is.data.frame(extracted)) {
      # For data frames, take the first column that's not ID
      pred_column <- if(ncol(extracted) > 1) 2 else 1
      pred_values <- extracted[[pred_column]]
    } else {
      # For vectors, use directly
      pred_values <- as.numeric(extracted)
    }
    
    # Convert to binary predictions based on threshold
    binary_preds <- ifelse(pred_values >= threshold, 1, 0)
    
    # Calculate correct classification percentage
    valid_indices <- !is.na(binary_preds) & !is.na(pa_vector)
    if(sum(valid_indices) == 0) {
      return(NA)  # No valid comparison points
    }
    
    percent_correct <- sum(binary_preds[valid_indices] == pa_vector[valid_indices]) / 
      sum(valid_indices) * 100
    
    return(percent_correct)
  }
  
  # Define threshold values for each species and method
  thresholds <- list(
    Calophyllum = list(
      "max(se+sp)" = thresh_calo_opt2,
      "minROCdist" = NA,  # Add this with NA value since it's not available
      "10th percentile" = thresh_calo_opt12
    ),
    Verticillium = list(
      "max(se+sp)" = NA,  # Add this with NA value since it's not available
      "minROCdist" = thresh_vertlept_opt4,
      "10th percentile" = thresh_vertlept_opt12
    )
  )
  
  
  # Define prediction raster data for each species
  rasters <- list(
    Calophyllum = ensname_proj_calo1_AUCopt2, # calo opt 2 ensemble
    Verticillium = ensname_proj_vertlept1_MCCopt4  # verleptop4 ensemble
  )
  
  # Define SDM objects for each species
  sdm_objects <- list(
    Calophyllum = calo_whole_thresh_0.7_forest_alg3, # calo sdm
    Verticillium = vertlept_2far_thresh_0.7_base_alg3  # vertlept sdm
  )
  
  # Create storage for results
  results <- data.frame(
    Optimization_Measure = character(),
    Metric = character(),
    Calophyllum = character(),
    Verticillium = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each threshold method
  methods <- c("max(se+sp)", "minROCdist", "10th percentile")
  method_labels <- c(
    "Max sum of sensitivity and specificity (\"max(se+sp)\") (opt 2)",
    "Minimum distance to ROC curve criterion (\"minROCdist\") (opt 4)",
    "10th percentile training presence (opt 12)"
  )
  
  for(i in 1:length(methods)) {
    method <- methods[i]
    method_label <- method_labels[i]
    
    # For each species, calculate stats for this threshold method
    species_results <- list()
    
    for(species in c("Calophyllum", "Verticillium")) {
      # Get the threshold value for this species and method
      threshold_value <- thresholds[[species]][[method]]
      
      # Initialize with NA values
      species_results[[species]] <- list(
        threshold = threshold_value,
        suitable_pixels = NA,
        percent_correct = NA
      )
      
      # Skip if threshold is NA
      if(is.na(threshold_value)) {
        next
      }
      
      # Get first layer of raster
      rast_layer1 <- tryCatch({
        rasters[[species]][[1]]
      }, error = function(e) {
        cat("Error getting raster layer for", species, ":", e$message, "\n")
        return(NULL)
      })
      
      # Skip if raster is NULL
      if(is.null(rast_layer1)) {
        next
      }
      
      # Calculate suitable pixels
      suitable_pixels <- tryCatch({
        count_suitable_cells(rast_layer1, threshold_value)
      }, error = function(e) {
        cat("Error counting suitable pixels for", species, "with method", method, ":", e$message, "\n")
        return(NA)
      })
      
      # Calculate classification accuracy
      percent_correct <- tryCatch({
        calculate_accuracy(
          rast_layer1,
          threshold_value,
          sdm_objects[[species]],
          1
        )
      }, error = function(e) {
        cat("Error calculating accuracy for", species, "with method", method, ":", e$message, "\n")
        return(NA)
      })
      
      # Update results
      species_results[[species]]$suitable_pixels <- suitable_pixels
      species_results[[species]]$percent_correct <- percent_correct
    }
    
    # Add threshold row
    results <- rbind(results, data.frame(
      Optimization_Measure = method_label,
      Metric = "Threshold value",
      Calophyllum = ifelse(is.na(species_results$Calophyllum$threshold), "n/a", 
                           format(round(species_results$Calophyllum$threshold, 5), nsmall = 5)),
      Verticillium = ifelse(is.na(species_results$Verticillium$threshold), "n/a", 
                            format(round(species_results$Verticillium$threshold, 5), nsmall = 5)),
      stringsAsFactors = FALSE
    ))
    
    # Add suitable pixels row
    results <- rbind(results, data.frame(
      Optimization_Measure = "",
      Metric = "No. of suitable pixels (1km)",
      Calophyllum = ifelse(is.na(species_results$Calophyllum$suitable_pixels), "n/a", 
                           format(species_results$Calophyllum$suitable_pixels, big.mark = ",")),
      Verticillium = ifelse(is.na(species_results$Verticillium$suitable_pixels), "n/a", 
                            format(species_results$Verticillium$suitable_pixels, big.mark = ",")),
      stringsAsFactors = FALSE
    ))
    
    # Add percent correct classification row
    results <- rbind(results, data.frame(
      Optimization_Measure = "",
      Metric = "% correct classification",
      Calophyllum = ifelse(is.na(species_results$Calophyllum$percent_correct), "n/a", 
                           paste0(format(round(species_results$Calophyllum$percent_correct, 2)), " %")),
      Verticillium = ifelse(is.na(species_results$Verticillium$percent_correct), "n/a", 
                            paste0(format(round(species_results$Verticillium$percent_correct, 2)), " %")),
      stringsAsFactors = FALSE
    ))
  }
  
  # Print the formatted table
  print(results)
  
  # output path
  out_path <- file.path(datafolder05_out, "threshold_classification")
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save the table to CSV
  write.csv(results, file.path(out_path, paste0("threshold_optimization_classification_results.csv")), row.names = FALSE)



#==============================================================================#
#                             6. MESS Analysis ----
#==============================================================================#

# MESS (Multivariate Environmental Similarity Surface) 
# (Zurell, Elith, and Schroeder 2012, Elith et al., 2010)
# Novel environments are conditions not realised in the sampled data but 
# are realised in the projection data. 
# If the entire niche of the species is encompassed by data, 
# then the model does not need to extrapolate 
# even if the projection data contain some novel environments. 
# Novel environments prove problematic if the niche is truncated in the sampled data 


# Add back in the ensemble projections if not already loaded in environment
## Load back in Calo ensemble projections ----
ensname_proj_calo1_AUCopt2 <- terra::rast(file.path(datafolder05_out, "ensproj_AUCopt2_calo_whole_thresh_0.7_forest_alg3.tif"))
## Load back in vertlept ensemble projections ----
ensname_proj_vertlept1_MCCopt4 <- terra::rast(file.path(datafolder05_out, "ensproj_MCCopt4_vertlept_2far_thresh_0.7_base_alg3.tif"))

# Function to extract reference environmental data from a model
extract_reference_env <- function(model) {
  # Get data from the sdm model
  species_data <- model@data
  
  # Extract feature data (environmental variables)
  feature_data <- species_data@features
  
  # Convert to data frame and return
  return(as.data.frame(feature_data))
}

# Function to calculate MESS using terra
calculate_mess <- function(ref_data, proj_rast) {
  # Create output SpatRaster with same dimensions as proj_rast
  mess_rast <- terra::rast(proj_rast[[1]])
  mess_rast <- setNames(mess_rast, "mess")
  
  # Convert projection raster to points for processing
  proj_points <- terra::as.points(proj_rast, na.rm=TRUE)
  proj_values <- terra::extract(proj_rast, proj_points)
  
  # Initialize results
  n_points <- nrow(proj_values)
  n_vars <- ncol(proj_values) - 1  # Subtract 1 for ID column
  
  # Initialize results matrices
  mess_values <- numeric(n_points)
  var_contributions <- matrix(0, nrow=n_points, ncol=n_vars)
  most_dissimilar <- integer(n_points)
  
  # Calculate MESS for each variable and point
  for (i in 1:n_points) {
    var_mess <- numeric(n_vars)
    
    for (j in 1:n_vars) {
      var_idx <- j + 1  # +1 because first column is ID
      var_name <- names(proj_values)[var_idx]
      
      # Get the value for this point
      point_value <- proj_values[i, var_idx]
      
      # Get reference distribution for this variable
      ref_values <- ref_data[, var_name]
      
      # Calculate the percentile of point_value in ref_values
      percent_below <- sum(ref_values < point_value) / length(ref_values) * 100
      
      # MESS calculation
      if (point_value < min(ref_values)) {
        # Below minimum
        var_mess[j] <- (point_value - min(ref_values)) / (max(ref_values) - min(ref_values)) * 100
      } else if (point_value > max(ref_values)) {
        # Above maximum
        var_mess[j] <- (point_value - max(ref_values)) / (max(ref_values) - min(ref_values)) * 100
      } else {
        # Within range
        var_mess[j] <- min(percent_below, 100 - percent_below) * 2
      }
      
      var_contributions[i, j] <- var_mess[j]
    }
    
    # Overall MESS is the minimum of the individual variable MESS values
    mess_values[i] <- min(var_mess)
    
    # Most dissimilar variable is the one with the lowest MESS value
    most_dissimilar[i] <- which.min(var_mess)
  }
  
  # Create data frame with coordinates and results
  coords <- terra::crds(proj_points)
  results_df <- data.frame(
    x = coords[,1],
    y = coords[,2],
    mess = mess_values
  )
  
  # Add variable contributions
  for (j in 1:n_vars) {
    var_name <- names(proj_values)[j+1]
    results_df[[paste0("mess_", var_name)]] <- var_contributions[, j]
  }
  
  # Add most dissimilar variable
  results_df$mod <- most_dissimilar
  
  # Rasterize results
  mess_result <- terra::rasterize(
    x = cbind(results_df$x, results_df$y),
    y = mess_rast,
    values = results_df$mess,
    fun = "mean"
  )
  
  return(mess_result)
}

# Function to match variable names between reference and projection data
match_variables <- function(ref_vars, proj_vars) {
  # Try exact matches first
  exact_matches <- intersect(ref_vars, proj_vars)
  
  if (length(exact_matches) > 0) {
    return(data.frame(
      ref_var = exact_matches,
      proj_var = exact_matches,
      stringsAsFactors = FALSE
    ))
  }
  
  # If no exact matches, try fuzzy matching
  matches <- data.frame(
    ref_var = character(0),
    proj_var = character(0),
    stringsAsFactors = FALSE
  )
  
  for (ref_var in ref_vars) {
    # Create clean versions of variable names for matching
    clean_ref <- tolower(gsub("[^[:alnum:]]", "", ref_var))
    
    for (proj_var in proj_vars) {
      clean_proj <- tolower(gsub("[^[:alnum:]]", "", proj_var))
      
      # Check if one name contains the other
      if (grepl(clean_ref, clean_proj, fixed=TRUE) || 
          grepl(clean_proj, clean_ref, fixed=TRUE)) {
        matches <- rbind(matches, 
                         data.frame(ref_var = ref_var, 
                                    proj_var = proj_var,
                                    stringsAsFactors = FALSE))
      }
    }
  }
  
  return(matches)
}

# Version 10: MESS visualisation function - carefully reviewed for syntax errors
create_mess_visualisations <- function(species, mess_raster, output_dir) {
  # Create output directory if it doesn't exist
  viz_dir <- file.path(output_dir, "visualisation")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract overall MESS values
  mess_overall <- mess_raster
  
  # Create binary novel environment map (MESS < 0 indicates extrapolation)
  novel_env <- mess_overall < 0
  novel_env <- setNames(novel_env, "novel_environment")
  
  # Calculate statistics on novel environments
  total_cells <- terra::global(novel_env, fun="notNA", na.rm=TRUE)$notNA
  novel_cells <- terra::global(novel_env, fun="sum", na.rm=TRUE)$sum
  pct_novel <- (novel_cells / total_cells) * 100
  
  # Create statistics data frame
  stats <- data.frame(
    Species = species,
    Total_Cells = total_cells,
    Novel_Env_Cells = novel_cells,
    Pct_Novel = pct_novel
  )
  
  # Save rasters
  terra::writeRaster(
    mess_overall,
    filename = file.path(output_dir, paste0(species, "_mess_analysis.tif")),
    overwrite = TRUE
  )
  
  terra::writeRaster(
    novel_env,
    filename = file.path(output_dir, paste0(species, "_novel_environment_binary.tif")),
    overwrite = TRUE
  )
  
  # Create PNG visualisation for MESS analysis
  png_file <- file.path(viz_dir, paste0(species, "_mess_plots.png"))
  png(png_file, width=800, height=400, res=100)
  
  # Set up plot parameters
  par(mfrow=c(1,2), mar=c(0,0,2,2), oma=c(0,0,0,0), xaxs="i", yaxs="i")
  
  # Plot 1: MESS values
  terra::plot(mess_overall, 
              main="MESS Analysis",
              axes=FALSE, box=FALSE,
              col=rev(terrain.colors(100)))
  
  # Plot 2: Binary novel environments
  terra::plot(novel_env, 
              main="Novel Environments (MESS < 0)",
              axes=FALSE, box=FALSE,
              col=c("grey", "red"),
              legend=TRUE)
  
  # Add MESS percentage text below the plot
  mtext(paste0("Novel environment: ", format(pct_novel, digits=2), "%"), 
        side=1, line=-1, cex=0.8)
  
  dev.off()
  
  # Try to create prediction visualisations
  create_prediction_viz <- function() {
    # Determine correct prediction file/object based on species
    pred_obj_name <- ifelse(species == "calo", "ensproj_calo1_AUCopt2", "ensproj_vertlept1_MCCopt4")
    pred_file_name <- ifelse(species == "calo", "ensname_proj_calo1_AUCopt2", "ensname_proj_vertlept1_MCCopt4")
    
    # Try to get the prediction object from memory first
    if (exists(pred_obj_name)) {
      cat("Using prediction object from memory:", pred_obj_name, "\n")
      pred_rast <- get(pred_obj_name)
      if (!inherits(pred_rast, "SpatRaster")) {
        pred_rast <- terra::rast(pred_rast)
      }
    } 
    # Then try to get the file path and load it
    else if (exists(pred_file_name)) {
      pred_path <- get(pred_file_name)
      
      # Handle S4 object issue
      if (is(pred_path, "S4")) {
        cat("Prediction path is an S4 object, attempting to convert to raster directly\n")
        pred_rast <- try(terra::rast(pred_path), silent = TRUE)
        if (inherits(pred_rast, "try-error")) {
          cat("Could not convert S4 object to raster\n")
          return(FALSE)
        }
      } else if (is.character(pred_path)) {
        cat("Attempting to load prediction from file:", pred_path, "\n")
        if (file.exists(pred_path)) {
          pred_rast <- terra::rast(pred_path)
        } else {
          cat("Prediction file does not exist:", pred_path, "\n")
          return(FALSE)
        }
      } else {
        cat("Prediction path is neither a character string nor an S4 object\n")
        return(FALSE)
      }
    } else {
      cat("No prediction object or file path found for", species, "\n")
      return(FALSE)
    }
    
    # Create masked prediction
    known_env_mask <- mess_overall >= 0  # Areas with non-negative MESS values
    masked_pred <- terra::mask(pred_rast, known_env_mask)
    
    # Save masked prediction
    mask_file <- file.path(output_dir, paste0(species, "_prediction_novel_env_masked.tif"))
    terra::writeRaster(masked_pred, filename=mask_file, overwrite=TRUE)
    
    # Create PNG for prediction visualisation
    png_pred_file <- file.path(viz_dir, paste0(species, "_prediction_plots.png"))
    png(png_pred_file, width=800, height=400, res=100)
    
    # Set up plot parameters
    par(mfrow=c(1,2), mar=c(0,0,2,2), oma=c(0,0,0,0), xaxs="i", yaxs="i")
    
    # Plot original prediction
    terra::plot(pred_rast,
                main="Original SDM Prediction",
                axes=FALSE, box=FALSE,
                col=rev(terrain.colors(100)))
    
    # Plot masked prediction
    terra::plot(masked_pred,
                main="Prediction (Known Environments Only)",
                axes=FALSE, box=FALSE,
                col=rev(terrain.colors(100)))
    
    dev.off()
    cat("Created prediction visualisations for", species, "\n")
    return(TRUE)
  }
  
  # Try to create prediction visualisations, with error handling
  tryCatch({
    result <- create_prediction_viz()
    if (!result) {
      cat("Could not create prediction visualisations for", species, "\n")
    }
  }, error = function(e) {
    cat("Error creating prediction visualisations for", species, ":", e$message, "\n")
  })
  
  # Return the statistics
  return(stats)
}

# Main MESS analysis function
run_mess_analysis <- function(model, env_stack, species, output_dir) {
  cat("\n==== Running MESS analysis for", species, "====\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract reference environmental data from model
  reference <- extract_reference_env(model)
  ref_var_names <- names(reference)
  
  cat("Reference variables:\n")
  print(ref_var_names)
  
  # Get projection variables
  proj_var_names <- names(env_stack)
  
  cat("Projection variables:\n")
  print(proj_var_names)
  
  # Match variables
  matched_vars <- match_variables(ref_var_names, proj_var_names)
  
  cat("Matched variables:\n")
  print(matched_vars)
  
  if (nrow(matched_vars) == 0) {
    stop("No matching variables found for ", species)
  }
  
  # Subset reference data and projection raster to matched variables
  ref_subset <- reference[, matched_vars$ref_var, drop=FALSE]
  proj_subset <- env_stack[[matched_vars$proj_var]]
  
  # Rename projection raster layers to match reference data
  names(proj_subset) <- matched_vars$ref_var
  
  # Calculate MESS
  cat("Calculating MESS...\n")
  mess_result <- calculate_mess(ref_subset, proj_subset)
  
  # Create visualisations and get statistics
  cat("Creating visualisations...\n")
  stats <- create_mess_visualisations(species, mess_result, output_dir)
  
  cat("MESS analysis completed for", species, "\n")
  cat("Percentage of novel environment:", round(stats$Pct_Novel, 2), "%\n\n")
  
  return(stats)
}

# Set output directory
mess_outdir <- file.path(datafolder, "05__model_prediction", "output", "mess_analysis")

# Run MESS analysis for both species
stats_list <- list()

# Calo
tryCatch({
  calo_stats <- run_mess_analysis(
    model = calo_whole_thresh_0.7_forest_alg3,
    env_stack = calo_whole_thresh_0.7_forest_curr_stack,
    species = "calo",
    output_dir = mess_outdir
  )
  stats_list[["calo"]] <- calo_stats
}, error = function(e) {
  cat("Error in Calo MESS analysis:", e$message, "\n")
})

# Vertlept
tryCatch({
  vertlept_stats <- run_mess_analysis(
    model = vertlept_2far_thresh_0.7_base_alg3,
    env_stack = vertlept_2far_thresh_0.7_base_curr_stack,
    species = "vertlept",
    output_dir = mess_outdir
  )
  stats_list[["vertlept"]] <- vertlept_stats
}, error = function(e) {
  cat("Error in Vertlept MESS analysis:", e$message, "\n")
})

# Combine statistics and write to CSV
if (length(stats_list) > 0) {
  novel_env_summary <- do.call(rbind, stats_list)
  write.csv(novel_env_summary, file=file.path(mess_outdir, "novel_environment_summary.csv"), row.names=FALSE)
  cat("\nNovel environment summary:\n")
  print(novel_env_summary)
} else {
  cat("\nNo MESS analyses were completed successfully.\n")
}

# Final message
cat("\nMESS analysis workflow completed.\n")
cat("Results are saved in:", mess_outdir, "\n")

# Tidy up
rm(list = ls())
gc()

#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#