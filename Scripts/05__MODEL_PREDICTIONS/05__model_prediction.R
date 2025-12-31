# ---
# title: "05__model_predictions"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-26"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# - Takes top trained models (model files in 03__model_building to project onto current climate, larger geography (incl. MDG)
#  then predicts through time under future climate scenarios between 2011-2040, 2041-2070, and 2071-2100
# - Uses top 3 performing models per species, see: "04__model_eval_perf\output\top_performing_models" for summary plots and full table
# - This set of code uses combined species occurrences (00), current environmental predictor variables, and future environment stacks
# - It outputs a set of SDM package prediction files
#
# HOW
# - sdm package "predict" function
# - Using first predictions of current environmental stack to look at feasibility i.e., is the wilt predicted to be in RNP?
# - ... and does Calo have the right distribution patterns/fits the presence and survey data.
# 
# WHY
# - Running and studying the predictions for current and future environments (climate scenarios, temporal ranges) to look at feasibility
# 
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

## Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

# Remove function that are not required for this script
rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, 
   align.forest.raster, align.raster, validate.processed.climate.files, calc_buffer_distance, calculate_metrics,
   create_sdm_ensemble, generate_ecoregion_bg_pts, plot_bg_comparison_calo, plot_bg_comparison_vert, process_2degfar_bg_pts, 
   random_bg_whole_area, rename_stack_layers, select07_cv, 
   analyze_elev_distribution, generate_buf_bg_pts, generate_dispersal_bg_pts, generate_elev_bg_pts, plot_bg_comparison)

library(sdm) 
installAll()

## Set up base folders ----

# Required derived data for predictions:
  ## model files (top 1 per spp)
  ## future environmental data cropped to MDG

# Base folders
datafolder <- file.path("//xxx")

# Inputs: Occurrence data
# Data folder relating to script 00 (species occurrences)
datafolder00 <- file.path(datafolder, "00__species_occurrence_cleaning")

# Inputs: Environmental data
# cropped, process environmental stacks to match final model predictors (top1) - manually copied into a separate folder due to naming clash "forest" and "forest_aspect"
env_dat_path <- file.path(datafolder, "05__model_prediction", "input", "variable_stacks", "top1")
calo_env_dat <- list.files(env_dat_path, pattern="calo", full.names = T)
vert_env_dat <- list.files(env_dat_path, pattern="vertlept", full.names = T)

# Inputs: final trained model per spp, model performance measures
# 03__ base output model training and performance outputs
datafolder03 <- file.path(datafolder, "03__model_training", "output")
final_models <- list.files(datafolder03, pattern=".sdm", full.names=T)

# Outputs: Prediction outputs
data__05_out <- file.path(datafolder, "05__model_prediction", "output")
# Create output directories
for(dir in c("model_evaluation", "current_projections/weighted_mean_ensembles/extra_settings", 
             "future_predictions/weighted_mean_ensembles/extra_settings", "figures/weighted_mean_ensembles/extra_settings")) {
  dir.create(file.path(data__05_out, dir), recursive = TRUE, showWarnings = FALSE)
}

#==============================================================================#
#                        1. Load in required data ----
#==============================================================================#


# Species occurrences
species_list <- c("calo", "vertlept")

# Set up time spans
# 2040, 2070, 2100 - 3 time stamps to predict
# Define future time periods
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

# Set up the climate pathway versions
# SSP126, SSP370, SSP585
climate_pathways <- c("ssp126", "ssp370", "ssp585")


#==============================================================================#
#                        2. Helper functions ----
#==============================================================================#


# Function to extract model information from file path
extract_model_info <- function(model_path) {
  model_name <- basename(model_path)
  species_name <- ifelse(grepl("calo", model_name), "calo", "vertlept")
  return(list(path = model_path, species = species_name))
}

# Function to extract scenario information from raster path
extract_scenario_info <- function(raster_path) {
  file_name <- basename(raster_path)
  
  # Extract climate scenario information
  is_current <- grepl("current_current", file_name)
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
    scenario = scenario,
    time_period = time_period,
    is_current = is_current
  ))
}


#==============================================================================#
#                        3. Run predictions ----
#==============================================================================#

# Loop round for each species, make the prediction for each environmental stack, 
# plot and save evaluation metric and prediction file.

# Process each species
for (species in species_list) {
  cat("\n\n=== Processing", species, "===\n")
  
  # Get the right environmental data and model
  env_files <- if (species == "calo") calo_env_dat else vert_env_dat
  
  # Find the final model for this species
  model_file <- final_models[grep(species, final_models)]
  
  if (length(model_file) == 0) {
    warning(paste("No model found for species:", species))
    next
  }
  
  if (length(model_file) > 1) {
    warning(paste("Multiple models found for species:", species, ". Using the first one."))
    model_file <- model_file[1]
  }
  
  cat("Loading model:", basename(model_file), "\n")
  load(model_file)
  model <- sdm_model 
  
  # Process each environmental stack
  for (env_file in env_files) {
    # Extract scenario information
    scenario_info <- extract_scenario_info(env_file)
    
    # Skip if we couldn't identify the scenario or time period
    if (is.na(scenario_info$scenario) || is.na(scenario_info$time_period)) {
      warning(paste("Could not identify scenario or time period for:", basename(env_file)))
      next
    }
    
    cat("\nPredicting for:", 
        "\n  Species:", species,
        "\n  Scenario:", scenario_info$scenario, 
        "\n  Time period:", scenario_info$time_period, "\n")
    
    # Load environmental data
    env_stack <- stack(env_file)
    
    # Create output filename
    output_dir <- ifelse(scenario_info$is_current, 
                      file.path(data__05_out, "current_projections", "weighted_mean_ensembles", "extra_settings"),
                      file.path(data__05_out, "future_predictions", "extra_settings"))
    dir.create(file.path(output_dir), recursive = TRUE, showWarnings = FALSE)
    
    
    output_filename <- paste0(species, "_", scenario_info$scenario, "_", 
                              scenario_info$time_period, "_prediction.tif")
    output_path <- file.path(output_dir, output_filename)
    
    # Make prediction using ensemble
    # ensemble using weighted averaging based on AUC statistic, collecting all uncertainty statistics
    # only includes models in the weighting that = auc > 0.7 & tss > 0.5
    # and optimum threshold criterion 2 (i.e., Max(spe+sen)) (maximises TSS) :  
    cat("Running ensemble prediction...\n")
    prediction <- ensemble(model, newdata = env_stack, filename = output_path,
                           setting=list(method=c('weighted',"uncertainty", "cv", "stdev", "ci"),
                                        stat='AUC',opt=2, expr='auc > 0.7 & tss > 0.5'), overwrite = TRUE)
    
    # Create a quick plot for visual inspection
    figure_path <- file.path(data__05_out, "figures", 
                             paste0(species, "_", scenario_info$scenario, "_", 
                                    scenario_info$time_period, "_predict.png"))
    
    png(figure_path, width = 800, height = 600)
    plot(prediction, main = paste(species, scenario_info$scenario, scenario_info$time_period))
    dev.off()
    
    cat("Prediction saved to:", output_path, "\n")
    cat("Figure saved to:", figure_path, "\n")
  }
  
  # Save model evaluation metrics to a file
  eval_path <- file.path(data__05_out, "model_evaluation", paste0(species, "_eval.csv"))
  if (file.exists(eval_path)) {
    cat("Model evaluation file already exists:", eval_path, "\n")
  } else {
    cat("Generating model evaluation metrics...\n")
    model_eval <- getEvaluation(model, stat = c("AUC", "TSS", "Kappa", "COR", "Deviance", "NMI", "phi", "ppv", "npv"))
    write.csv(model_eval, eval_path, row.names = FALSE)
    cat("Evaluation metrics saved to:", eval_path, "\n")
  }
}

cat("\n=== Predictios completed ===\n")

# # Clean up
gc()
rm(list = ls())

#==============================================================================#
#                        ----  End of workflow ----
#==============================================================================#