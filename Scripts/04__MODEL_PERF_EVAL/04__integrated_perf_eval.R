# ---
# title: "04__integrated_perf_eval.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-11"
# update: "2025-12-31"
# ---

# WORK FLOW NOTES:

# WHAT
# This script integrates SDM model evaluation, performance analysis, and 
# visualisation across individual models and ensembles.
# 
# Workflow:
# 1. Set up environment and paths
# 2. Define helper functions for evaluation
# 3. Process individual models and calculate metrics
# 4. Evaluate ensemble models 
# 5. Rank species ensembles against each other
# 6. Create visualisation outputs
# 7. Export results for publication

# Libraries
library(sdm)        # For SDM model handling
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping
library(ggplot2)    # For visualisation
library(viridis)    # For colourblind-friendly palettes
library(grid)       # For inset plots
library(plotly)     # For interactive plots
library(htmlwidgets) # For saving web-based visualisations
library(tibble)     # For improved data frames


#==============================================================================#
#                           ----  1. Workspace set up ----
#==============================================================================#

## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

# Define paths
datafolder <- file.path("//xxx")

# Input data folders
datafolder03 <- file.path(datafolder, "03__model_training", "output")
datafolder02 <- file.path(datafolder, "02__variable_selection", "output", "selected_vars", "sdm_variable_sets")

# Output directories
data__04_out <- file.path(datafolder, "04__model_eval_perf", "output")
outdir__04_top <- file.path(data__04_out, "top_performing_models")

# Create output directories
for (dir in c(
  file.path(data__04_out, "plots/pdp"), 
  file.path(data__04_out, "plots/pdp_3d"), 
  file.path(data__04_out, "plots/importance"), 
  file.path(data__04_out, "metrics"),
  file.path(outdir__04_top),
  file.path(outdir__04_top, "species_ensemble_rankings"),
  file.path(outdir__04_top, "explained_deviance")
)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Define the species list
species_list <- c("calo", "vertlept")
species_display_names <- c(
  "calo" = "Calophyllum paniculatum",
  "vertlept" = "Verticillium/Leptographium"
)

# Define model algorithm sets
alg_sets <- list(
  alg3 = c("rf", "brt", "maxent"),
  alg4 = c("rf", "brt", "maxent", "gam")
)

# Environmental data lookup with all unique variables
bio_lookup <- c(
  # Temperature variables
  "Mean.Temp.Wet.Q" = "Mean Temperature of Wettest Quarter",
  "Mean.Temp.Dry.Q" = "Mean Temperature of Driest Quarter",
  "Mean.Temp.Cold.Q" = "Mean Temperature of Coldest Quarter",
  "Min.Temp.Coldest" = "Minimum Temperature of Coldest Month",
  "Ann.Temp.Range" = "Annual Temperature Range",
  "Temp.Seasonality" = "Temperature Seasonality",
  "Diurnal.Range" = "Mean Diurnal Range",
  
  # Precipitation variables
  "Precip.Wet.Month" = "Precipitation of Wettest Month",
  "Precip.Dry.Month" = "Precipitation of Driest Month",
  "Precip.Cold.Q" = "Precipitation of Coldest Quarter",
  "Precip.Seasonality" = "Precipitation Seasonality",
  
  # Land cover variables
  "Forest.Cover" = "Forest Cover",
  
  # Additional variables that might be in file names
  "aspect" = "Aspect",
  "forest" = "Forest Cover",
  "forest_aspect" = "Forest and Aspect"
)

# Evaluation metrics lookup and descriptions
eval_metrics_desc <- c(
  "AUC" = "Area Under ROC Curve (0-1, >0.8 excellent)",
  "TSS" = "True Skill Statistic (-1 to +1, >0.6 good)",
  "COR" = "Point-biserial correlation",
  "Deviance" = "Residual deviance",
  "MCC" = "Matthews Correlation Coefficient",
  "F1" = "F1 Score"
)



#==============================================================================#
#                     ---- 2. Helper Functions ----                            
#==============================================================================#


# Memory management function to be called at appropriate points
clear_memory <- function() {
  gc(verbose = FALSE, full = TRUE)
}

# Function to clean and format variable labels
clean_labels <- function(text) {
  # Initialize output vector
  clean_text <- character(length(text))
  
  # Process each element
  for(i in seq_along(text)) {
    if(text[i] %in% names(bio_lookup)) {
      clean_text[i] <- bio_lookup[text[i]]
    } else {
      clean_text[i] <- gsub("[_.]", " ", text[i])
    }
  }
  
  return(clean_text)
}

# Function to create 3D partial dependency plots
create_3d_pdp <- function(response_curves, var1, var2) {
  # Extract response data
  data1 <- response_curves@response[[var1]]
  data2 <- response_curves@response[[var2]]
  
  # Calculate mean responses for each variable
  response_cols1 <- grep("ID-", names(data1), value = TRUE)
  response_cols2 <- grep("ID-", names(data2), value = TRUE)
  
  # Create grid of values
  x <- data1[[var1]]  # Variable 1 values
  y <- data2[[var2]]  # Variable 2 values
  
  # Calculate mean response surface
  z <- matrix(NA, nrow = length(x), ncol = length(y))
  for(i in seq_along(x)) {
    for(j in seq_along(y)) {
      resp1 <- mean(as.numeric(data1[i, response_cols1]))
      resp2 <- mean(as.numeric(data2[j, response_cols2]))
      z[i,j] <- (resp1 + resp2)/2
    }
  }
  
  # Create the 3D plot
  fig <- plotly::plot_ly() %>%
    add_surface(
      x = ~x,
      y = ~y,
      z = ~z,
      colourscale = "Viridis"
    ) %>%
    layout(
      scene = list(
        xaxis = list(title = clean_labels(var1)),
        yaxis = list(title = clean_labels(var2)),
        zaxis = list(title = "Response")
      ),
      title = list(
        text = paste("3D Response Surface:", 
                     clean_labels(var1), 
                     "vs", 
                     clean_labels(var2))
      )
    )
  
  return(fig)
}

# Function to parse model filename to extract components
parse_model_name <- function(model_name) {
  # Example: calo_ecor_thresh_0.5_forest_alg3.sdm
  parts <- strsplit(tools::file_path_sans_ext(model_name), "_")[[1]]
  
  # Extract species, bg_method, threshold, var_set, alg_set
  species <- parts[1]
  bg_method <- parts[2]
  
  # Determine threshold
  thresh_idx <- which(parts == "thresh")
  if(length(thresh_idx) > 0) {
    threshold <- paste("thresh", parts[thresh_idx + 1], sep = "_")
  } else if("expert" %in% parts) {
    threshold <- "expert"
  } else {
    threshold <- NA
  }
  
  # Determine variable set (base, forest, aspect, forest_aspect)
  var_set <- "base"  # Default
  if("forest_aspect" %in% parts || (all(c("forest", "aspect") %in% parts) && !("base" %in% parts))) {
    var_set <- "forest_aspect"
  } else if("forest" %in% parts) {
    var_set <- "forest"
  } else if("aspect" %in% parts) {
    var_set <- "aspect"
  }
  
  # Extract algorithm set
  alg_idx <- which(grepl("alg", parts))
  alg_set <- parts[alg_idx]
  
  return(list(
    species = species,
    bg_method = bg_method,
    threshold = threshold,
    var_set = var_set,
    alg_set = alg_set
  ))
}

# Function to load all .sdm models from a directory
load_all_sdm_models <- function(species_directory) {
  # Get the path to the directory
  algsets_dir <- file.path(species_directory)
  
  # List all .sdm files (non-recursive - only get top-level files)
  sdm_files <- list.files(algsets_dir, pattern = "\\.sdm$", 
                          recursive = FALSE, full.names = TRUE)
  
  # Initialize empty list to store models
  sdm_models_list <- list()
  
  # Load each model and store in the list
  for (file_path in sdm_files) {
    # Extract the model name from the file path (without extension)
    model_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Load the model
    message("Loading model: ", model_name)
    load(file_path)
    
    # Store in list with the model name
    sdm_models_list[[model_name]] <- sdm_model
    
    # Clean up memory after processing each ensemble
    rm(sdm_model)
    clear_memory()
  }
  
  # Combine all results into a single data frame
  combined_results <- do.call(rbind, all_ensemble_results)
  
  # Calculate summary statistics with explicit naming for easier access
  summary_stats <- combined_results %>%
    group_by(algorithm, ensemble_name) %>%
    summarize(
      # For AUC
      AUC_mean = mean(AUC, na.rm = TRUE),
      AUC_sd = sd(AUC, na.rm = TRUE),
      AUC_min = min(AUC, na.rm = TRUE),
      AUC_max = max(AUC, na.rm = TRUE),
      
      # For TSS
      TSS_mean = mean(TSS, na.rm = TRUE),
      TSS_sd = sd(TSS, na.rm = TRUE),
      TSS_min = min(TSS, na.rm = TRUE),
      TSS_max = max(TSS, na.rm = TRUE),
      
      # For explained deviance
      explained_deviance_mean = mean(explained_deviance, na.rm = TRUE),
      explained_deviance_sd = sd(explained_deviance, na.rm = TRUE),
      explained_deviance_min = min(explained_deviance, na.rm = TRUE),
      explained_deviance_max = max(explained_deviance, na.rm = TRUE),
      
      # For Kappa
      Kappa_mean = mean(Kappa, na.rm = TRUE),
      Kappa_sd = sd(Kappa, na.rm = TRUE),
      Kappa_min = min(Kappa, na.rm = TRUE),
      Kappa_max = max(Kappa, na.rm = TRUE),
      
      # Calculate combined scores using the specified formula
      Score_AUC_TSS = 0.6*AUC_mean + 0.4*TSS_mean - 0.5*(0.6*AUC_sd + 0.4*TSS_sd),
      
      # Number of models in this group
      n_models = n(),
      
      .groups = "drop"
    )
  
  # Add species column to summary
  summary_stats$species <- species_name
  
  # Save summary statistics
  if (!is.null(output_dir)) {
    summary_path <- file.path(output_dir, paste0(safe_species_name, "_sdm_summary_explicit.csv"))
    write.csv(summary_stats, summary_path, row.names = FALSE)
    message("Detailed summary statistics with scores saved to: ", summary_path)
    
    # Also save a ranked version based on the combined score
    ranked_summary <- summary_stats %>%
      arrange(desc(Score_AUC_TSS))
    
    ranked_path <- file.path(output_dir, paste0(safe_species_name, "_sdm_summary_ranked.csv"))
    write.csv(ranked_summary, ranked_path, row.names = FALSE)
    message("Ranked summary statistics saved to: ", ranked_path)
  }
  
  # Clean up memory
  clear_memory()
  
  # Return all results plus the enhanced summary
  return(list(
    individual_results = all_ensemble_results,
    combined_results = combined_results,
    summary_stats = summary_stats,
    ranked_summary = summary_stats %>% arrange(desc(Score_AUC_TSS))
  ))
}
#rm(sdm_model)


# Function to calculate null deviance (for explained deviance)
calculate_null_deviance <- function(sdm_model, response = "occ", algorithm, model_id) {
  # Extract the specific model from the ensemble with proper nesting
  model_data <- sdm_model@models[[response]][[algorithm]][[as.character(model_id)]]
  
  # Extract observed values (presence/absence) from the training data
  observed <- model_data@evaluation$training@observed
  
  # Calculate prevalence (proportion of presences)
  prevalence <- mean(observed)
  
  # For a binary response, the null model predicts the prevalence for all observations
  predicted_null <- rep(prevalence, length(observed))
  
  # Calculate null deviance using binomial deviance formula                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
  # To avoid log(0) issues, add a small constant
  epsilon <- 1e-10
  predicted_null_adj <- pmax(pmin(predicted_null, 1-epsilon), epsilon)
  
  # Calculate null deviance
  null_deviance <- -2 * sum(
    observed * log(predicted_null_adj) + 
      (1 - observed) * log(1 - predicted_null_adj)
  )
  
  return(null_deviance)
}

# Check if a model has already been processed
check_model_processed <- function(model_path, target_alg = "alg4") {
  model_name <- basename(model_path)
  model_components <- parse_model_name(model_name)
  
  # Check if this is an alg4 model (if we're targeting those)
  if (!is.null(target_alg) && model_components$alg_set != target_alg) {
    return(TRUE)  # Skip non-target models
  }
  
  # Check if results file already exists
  output_dir <- file.path(data__04_out, model_components$species, "metrics")
  result_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(model_name), "_results.rds"))
  
  # Return TRUE if the file exists (already processed)
  return(file.exists(result_path))
}

#==============================================================================#
#                 ---- 3. Individual Model Processing ----                     
#==============================================================================#


# Process a single model to extract metrics, create plots, etc.
process_one_model <- function(model_path, skip_3d_pdp = TRUE) {
  # Load model and extract name
  load(model_path)
  model_name <- basename(model_path)
  cat("Processing:", model_name, "\n")
  
  # Parse model name components
  model_components <- parse_model_name(model_name)
  
  # 1. Get evaluation metrics
  eval_metrics <- getEvaluation(sdm_model)
  metrics_summary <- data.frame(
    Model = model_name,
    Metric = names(eval_metrics_desc),
    Mean = sapply(names(eval_metrics_desc), function(m) mean(eval_metrics[[m]], na.rm = TRUE)),
    SD = sapply(names(eval_metrics_desc), function(m) sd(eval_metrics[[m]], na.rm = TRUE)),
    Description = unname(eval_metrics_desc),
    Species = model_components$species,
    BG_Method = model_components$bg_method,
    Threshold = model_components$threshold,
    Var_Set = model_components$var_set,
    Alg_Set = model_components$alg_set
  )
  
  # 2. Get variable importance - modified to handle the structure correctly
  var_imp <- getVarImp(sdm_model)
  # Safe handling for variable importance in case the structure varies
  var_imp_df <- NULL
  tryCatch({
    var_imp_df <- data.frame(
      variable = var_imp@varImportanceMean$corTest$variables,
      importance = var_imp@varImportanceMean$corTest$corTest,
      lower = var_imp@varImportanceMean$corTest$lower,
      upper = var_imp@varImportanceMean$corTest$upper
    )
  }, error = function(e) {
    warning(paste("Error extracting variable importance:", e$message))
    # Fallback method to extract variable importance
    var_imp_df <- data.frame(
      variable = names(var_imp@varImportanceMean),
      importance = unlist(var_imp@varImportanceMean),
      lower = NA,
      upper = NA
    )
  })
  
  # 3. Get response curves with explicit model IDs
  model_ids <- getModelId(sdm_model)
  response_curves <- getResponseCurve(sdm_model, id = model_ids)
  
  # 4. Create plots
  # a. Variable Importance Plot
  if(!is.null(var_imp_df)) {
    p_varimp <- ggplot(var_imp_df, aes(x = reorder(variable, importance), y = importance)) +
      geom_col(fill = "steelblue") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)
      ) +
      labs(
        title = "Variable Importance",
        subtitle = paste("Model:", gsub(".sdm", "", model_name)),
        x = "Environmental Variables",
        y = "Importance Score",
        caption = paste("AUC =", round(metrics_summary$Mean[metrics_summary$Metric == "AUC"], 3),
                        "±", round(metrics_summary$SD[metrics_summary$Metric == "AUC"], 3))
      ) +
      scale_x_discrete(labels = clean_labels)
    
    # Save variable importance plot
    output_dir <- file.path(data__04_out, model_components$species)
    ggsave(
      file.path(output_dir, "plots", "importance",
                paste0(tools::file_path_sans_ext(model_name), "_varImp.png")),
      p_varimp, dpi = 300, width = 8, height = 6
    )
  }
  
  # b. Response Curves (2D)
  for(var in response_curves@variables) {
    var_data <- response_curves@response[[var]]
    response_cols <- grep("ID-", names(var_data), value = TRUE)
    var_data$mean_response <- rowMeans(var_data[, response_cols])
    var_data$sd_response <- apply(var_data[, response_cols], 1, sd)
    
    p_response <- ggplot(var_data, aes_string(x = var)) +
      geom_ribbon(
        aes(ymin = mean_response - 1.96*sd_response,
            ymax = mean_response + 1.96*sd_response),
        alpha = 0.1, fill = "grey50"
      ) +
      geom_ribbon(
        aes(ymin = mean_response - sd_response,
            ymax = mean_response + sd_response),
        alpha = 0.2, fill = "steelblue"
      ) +
      geom_line(aes(y = mean_response), linewidth = 1, colour = "steelblue") +
      theme_minimal() +
      labs(
        title = clean_labels(var),
        subtitle = paste("Model:", gsub(".sdm", "", model_name)),
        x = clean_labels(var),
        y = "Response",
        caption = "Dark band: ±1 SE, Light band: 95% CI"
      )
    
    output_dir <- file.path(data__04_out, model_components$species)
    ggsave(
      file.path(output_dir, "plots", "pdp",
                paste0(tools::file_path_sans_ext(model_name), "_", var, "_response.png")),
      p_response, dpi = 300, width = 8, height = 6
    )
  }
  
  # 5. Save results
  results <- list(
    metrics = metrics_summary
  )
  
  output_dir <- file.path(data__04_out, model_components$species)
  saveRDS(results, 
          file.path(output_dir, "metrics",
                    paste0(tools::file_path_sans_ext(model_name), "_results.rds")))
  
  # Clean up memory
  clear_memory()
  
  return(metrics_summary)
}


# Process all individual models
process_all_models <- function(target_alg = "alg4", force_reprocess = FALSE) {
  # Initialize results list
  all_results <- list()
  
  # Process each species
  for(species in species_list) {
    cat("\nProcessing", species, "models...\n")
    
    # Create species-specific results directory
    species_dir <- file.path(data__04_out, species)
    for(dir in c("plots/pdp", "plots/importance", "metrics")) {
      dir.create(file.path(species_dir, dir), recursive = TRUE, showWarnings = FALSE)
    }
    
    # Get all model files for the species (non-recursive to avoid old files)
    model_files <- list.files(
      path = file.path(datafolder03, species),
      pattern = "\\.sdm$",
      recursive = FALSE,
      full.names = TRUE
    )
    
    if(length(model_files) == 0) {
      warning(paste("No model files found for species:", species))
      next
    }
    
    # Filter for target algorithm if specified
    if (!is.null(target_alg)) {
      target_files <- grep(target_alg, model_files, value = TRUE)
      if (length(target_files) > 0) {
        model_files <- target_files
        cat("Filtering for", target_alg, "models:", length(model_files), "found\n")
      } else {
        cat("No", target_alg, "models found for species:", species, "\n")
      }
    }
    
    # Check which models need processing
    models_to_process <- character(0)
    for(model_path in model_files) {
      if(force_reprocess || !check_model_processed(model_path, target_alg)) {
        models_to_process <- c(models_to_process, model_path)
      }
    }
    
    cat("Found", length(models_to_process), "models to process (out of", length(model_files), "total)\n")
    
    # Process each model that needs it
    species_results <- list()
    for(model_path in models_to_process) {
      model_name <- basename(model_path)
      cat("  Processing:", model_name, "\n")
      
      tryCatch({
        result <- process_one_model(model_path, skip_3d_pdp = TRUE)
        species_results[[model_name]] <- result
      }, error = function(e) {
        warning(paste("Error processing", model_name, ":", e$message))
      })
    }
    
    # Load previously processed results for this species
    metrics_files <- list.files(
      path = file.path(species_dir, "metrics"),
      pattern = "_results\\.rds$",
      full.names = TRUE
    )
    
    for(metrics_file in metrics_files) {
      model_name <- gsub("_results.rds$", ".sdm", basename(metrics_file))
      if(!(model_name %in% names(species_results))) {
        # Only load if we haven't already processed this model
        tryCatch({
          result <- readRDS(metrics_file)$metrics
          species_results[[model_name]] <- result
        }, error = function(e) {
          warning(paste("Error loading", metrics_file, ":", e$message))
        })
      }
    }
    
    # Combine results for this species if any were successfully processed
    if(length(species_results) > 0) {
      combined_metrics <- do.call(rbind, species_results)
      
      # Save combined results
      saveRDS(combined_metrics, 
              file.path(species_dir, "metrics", "all_model_metrics.rds"))
      write.csv(combined_metrics,
                file.path(species_dir, "metrics", "all_model_metrics.csv"),
                row.names = FALSE)
      
      # Create performance comparison plot for AUC metric
      auc_metrics <- combined_metrics[combined_metrics$Metric == "AUC", ]
      
      p_comparison <- ggplot(auc_metrics, 
                             aes(x = Var_Set, y = Mean, colour = BG_Method)) +
        geom_point(position = position_dodge(width = 0.3)) +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                      position = position_dodge(width = 0.3), width = 0.2) +
        facet_grid(Alg_Set ~ Threshold) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste(species, "Model Performance Comparison"),
             x = "Habitat Type",
             y = "AUC Score",
             colour = "Background Method")
      
      ggsave(file.path(species_dir, "plots", paste0(species, "_performance_comparison.png")),
             p_comparison, dpi = 300, width = 12, height = 8)
      
      all_results[[species]] <- combined_metrics
    }
    
    # Clean up memory after processing each species
    clear_memory()
  }
  
  # Create an overall summary across all species
  if(length(all_results) > 0) {
    all_combined <- do.call(rbind, all_results)
    saveRDS(all_combined, file.path(data__04_out, "all_species_metrics.rds"))
    write.csv(all_combined, file.path(data__04_out, "all_species_metrics.csv"), row.names = FALSE)
    
    # Create overall performance comparison
    auc_all <- all_combined[all_combined$Metric == "AUC", ]
    p_overall <- ggplot(auc_all, 
                        aes(x = Var_Set, y = Mean, colour = Species)) +
      geom_point(position = position_dodge(width = 0.3)) +
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                    position = position_dodge(width = 0.3), width = 0.2) +
      facet_grid(Alg_Set ~ Threshold) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Overall Model Performance Comparison",
           x = "Habitat Type",
           y = "AUC Score",
           colour = "Species")
    
    ggsave(file.path(data__04_out, "overall_performance_comparison.png"),
           p_overall, dpi = 300, width = 12, height = 8)
  }
  
  # Clean up memory
  clear_memory()
  
  return(all_results)
}


#==============================================================================#
#                    ---- 4. Ensemble Evaluation ----                           
#==============================================================================#


# Create evaluation table for SDM ensembles
create_sdm_evaluation_table <- function(sdm_model, response = "occ") {
  # Get all available algorithms
  algorithms <- names(sdm_model@models[[response]])
  
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop through each algorithm
  for (algo in algorithms) {
    # Get all model IDs for this algorithm
    model_ids <- names(sdm_model@models[[response]][[algo]])
    
    # Loop through each model ID
    for (mid in model_ids) {
      # Wrap in try-catch to handle potential errors
      tryCatch({
        # Extract the model
        model_data <- sdm_model@models[[response]][[algo]][[mid]]
        
        # Calculate null deviance and explained deviance
        null_dev <- calculate_null_deviance(sdm_model, response, algo, mid)
        residual_dev <- model_data@evaluation$training@statistics$Deviance
        explained_dev <- 1 - (residual_dev / null_dev)
        
        # Extract other standard metrics from the model
        auc <- model_data@evaluation$training@statistics$AUC
        correlation <- model_data@evaluation$training@statistics$COR[1]  # First element is correlation value
        prevalence <- model_data@evaluation$training@statistics$Prevalence
        
        # Get threshold-based metrics (using the prevalence threshold)
        thresh_idx <- which(model_data@evaluation$training@threshold_based$criteria == "prevalence")
        if (length(thresh_idx) > 0) {
          thresh_metrics <- model_data@evaluation$training@threshold_based[thresh_idx, ]
          sensitivity <- thresh_metrics$sensitivity
          specificity <- thresh_metrics$specificity
          tss <- thresh_metrics$TSS
          kappa <- thresh_metrics$Kappa
        } else {
          # Default values if prevalence threshold isn't available
          sensitivity <- NA
          specificity <- NA
          tss <- NA
          kappa <- NA
        }
        
        # Store all results in a list
        all_results[[length(all_results) + 1]] <- list(
          response = response,
          algorithm = algo,
          model_id = as.numeric(mid),
          null_deviance = null_dev,
          residual_deviance = residual_dev,
          explained_deviance = explained_dev,
          AUC = auc,
          correlation = correlation,
          prevalence = prevalence,
          sensitivity = sensitivity,
          specificity = specificity,
          TSS = tss,
          Kappa = kappa
        )
      }, error = function(e) {
        warning(paste("Error processing model ID", mid, "in algorithm", algo, ":", e$message))
      })
    }
  }
  
  # Check if any results were generated
  if (length(all_results) == 0) {
    stop("No results could be generated. Check if the ensemble models have the expected structure.")
  }
  
  # Convert list of results to a data frame - use tibble instead of data.frame
  results_df <- do.call(rbind, lapply(all_results, function(x) {
    tibble::tibble(
      response = x$response,
      algorithm = x$algorithm,
      model_id = x$model_id,
      null_deviance = x$null_deviance,
      residual_deviance = x$residual_deviance,
      explained_deviance = x$explained_deviance,
      AUC = x$AUC,
      correlation = x$correlation,
      prevalence = x$prevalence,
      sensitivity = x$sensitivity,
      specificity = x$specificity,
      TSS = x$TSS,
      Kappa = x$Kappa
    )
  }))
  
  # Sort by algorithm and model_id - using dplyr for consistency
  results_df <- results_df %>%
    arrange(algorithm, model_id)
  
  return(results_df)
}

# Check if ensemble results already exist
check_ensemble_processed <- function(ensemble_name, species_name, output_dir) {
  # Create a safe filename version of the species name
  safe_species_name <- gsub("[^a-zA-Z0-9]", "_", species_name)
  
  # Check if results file already exists
  results_path <- file.path(output_dir, paste0(safe_species_name, "_", ensemble_name, ".csv"))
  
  return(file.exists(results_path))
}

# Evaluate ensemble models
evaluate_sdm_ensembles <- function(ensemble_models, output_dir = NULL, species_name = NULL, force_reprocess = FALSE) {
  # Ensure output directory exists if specified
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Set default species name if not provided
  if (is.null(species_name)) {
    species_name <- "species"
  }
  
  # Create a safe filename version of the species name
  safe_species_name <- gsub("[^a-zA-Z0-9]", "_", species_name)
  
  # Initialize list to store results from all ensembles
  all_ensemble_results <- list()
  
  # Process each ensemble model
  for (i in seq_along(ensemble_models)) {
    # Get ensemble name
    if (!is.null(names(ensemble_models)) && names(ensemble_models)[i] != "") {
      ensemble_name <- names(ensemble_models)[i]
    } else {
      ensemble_name <- paste0("ensemble_", i)
    }
    
    # Check if this ensemble has already been processed
    if (!force_reprocess && !is.null(output_dir) && 
        check_ensemble_processed(ensemble_name, species_name, output_dir)) {
      message("Ensemble ", ensemble_name, " for species ", species_name, " already processed. Skipping.")
      
      # Load existing results
      file_path <- file.path(output_dir, paste0(safe_species_name, "_", ensemble_name, ".csv"))
      eval_table <- read.csv(file_path)
      all_ensemble_results[[ensemble_name]] <- eval_table
      
      next
    }
    
    message("Processing ensemble model: ", ensemble_name, " for species: ", species_name)
    
    # Extract the model
    sdm_model <- ensemble_models[[i]]
    
    # Get evaluation table for this ensemble
    eval_table <- create_sdm_evaluation_table(sdm_model)
    
    # Add species and ensemble name columns
    eval_table$species <- species_name
    eval_table$ensemble_name <- ensemble_name
    
    # Store in results list
    all_ensemble_results[[ensemble_name]] <- eval_table
    
    # Save individual results if output directory is specified
    if (!is.null(output_dir)) {
      file_path <- file.path(output_dir, paste0(safe_species_name, "_", ensemble_name, ".csv"))
      write.csv(eval_table, file_path, row.names = FALSE)
      message("Results saved to: ", file_path)
    }
  }}



#==============================================================================#
#                 ---- 5. Visualisation and Outputs ----                      
#==============================================================================#
  
# Process model metrics for individual and ensemble visualisation
process_models_metrics <- function(force_reprocess = FALSE) {
    # Get model metrics from saved outputs
    model_metrics_path <- file.path(data__04_out, "all_species_metrics_wide_clean.csv")
    
    # Check if file already exists and we're not forcing reprocessing
    if (!force_reprocess && file.exists(model_metrics_path)) {
      message("Loading model metrics from existing file:", model_metrics_path)
      model_metrics <- read.csv(model_metrics_path)
    } else {
      # If not, we'll need to create it from the individual species metrics
      cat("Creating combined metrics file...\n")
      
      # Get individual species metrics
      calo_metrics_path <- file.path(data__04_out, "calo", "metrics", "all_model_metrics_wide.csv")
      vertlept_metrics_path <- file.path(data__04_out, "vertlept", "metrics", "all_model_metrics_wide.csv")
      
      if (!file.exists(calo_metrics_path) || !file.exists(vertlept_metrics_path)) {
        warning("Individual species metrics not found. Please run the individual model processing first.")
        return(NULL)
      }
      
      # Load individual metrics
      calo_subset <- read.csv(calo_metrics_path)
      vl_subset <- read.csv(vertlept_metrics_path)
      
      # Combine into one table
      all_species_metrics_wide <- rbind(calo_subset, vl_subset)
      
      # Clean up column names for presentation
      all_species_metrics_wide <- all_species_metrics_wide %>%
        mutate(
          # Clean Species names
          Species = case_when(
            Species == "calo" ~ "Calophyllum paniculatum",
            Species == "vertlept" ~ "Verticillium/Leptographium",
            TRUE ~ Species
          ),
          # Clean BG_Method
          BG_Method = case_when(
            BG_Method == "ecor" ~ "Ecoregion",
            BG_Method == "whole" ~ "Whole",
            TRUE ~ BG_Method
          ),
          # Extract only the numeric part from Threshold
          Threshold = case_when(
            grepl("^thresh_", Threshold) ~ gsub("thresh_", "", Threshold),
            TRUE ~ Threshold
          ),
          Var_Set = case_when(
            Var_Set == "forest_aspect" ~ "Forest + Aspect",
            Var_Set == "aspect" ~ "Aspect",
            Var_Set == "base" ~ "Base",
            TRUE ~ Var_Set
          ),
          Alg_Set = case_when(
            grepl("^alg", Alg_Set) ~ paste("Alg", gsub("alg", "", Alg_Set)),
            TRUE ~ Alg_Set
          )
        )
      
      # Write the processed file
      write.csv(all_species_metrics_wide, file = model_metrics_path, row.names = FALSE)
      model_metrics <- all_species_metrics_wide
    }
    
    # Clean up the data for visualisation
    model_data <- model_metrics %>%
      # Create factors for better plot ordering
      mutate(
        Species = factor(Species),
        BG_Method = factor(BG_Method, levels = c("Whole", "Ecoregion", "2far")),
        Threshold = factor(Threshold, levels = c("0.5", "0.7")),
        Var_Set = factor(Var_Set, levels = c("Base", "Forest", "Aspect", "Forest + Aspect")),
        Alg_Set = factor(Alg_Set)
      )
    
    # Calculate a global SD for each metric 
    global_sd <- model_data %>%
      summarize(
        SD_AUC = sd(AUC, na.rm = TRUE),
        SD_TSS = sd(TSS, na.rm = TRUE)
      )
    
    # Use SD to create multi-metric scoring/ranking
    model_wide <- model_data %>%
      mutate(
        SD_AUC = global_sd$SD_AUC,
        SD_TSS = global_sd$SD_TSS,
        Score = 0.6*AUC + 0.4*TSS - 0.5*(0.6*SD_AUC + 0.4*SD_TSS)
      )
    
    # Check if output directory exists
    if (!dir.exists(outdir__04_top)) {
      dir.create(outdir__04_top, recursive = TRUE)
    }
    
    # Save the scored/ranked list of all models to file (both species combined into one table)
    write.csv(model_wide, file=(file.path(outdir__04_top, "model_full_rankings_AUC_TSS_combined_score.csv")))
    
    # Save only top 3 models per species
    species_list <- unique(model_wide$Species)
    
    for (species in species_list) {
      # Step 1: Filter rows per species
      subset <- model_wide[model_wide$Species == species, ]
      
      # Step 2: Order by "Score" in descending order
      sorted <- subset[order(subset$Score, decreasing = TRUE), ]
      
      # Step 3: Keep only the top 3 rows
      top_3 <- sorted[1:3, ]
      safe_species_name <- gsub("[^a-zA-Z0-9]", "_", species)
      write.csv(top_3, file=file.path(outdir__04_top, paste0(safe_species_name,"_top3_models.csv")))
    }
    
    # Clean up memory
    clear_memory()
    
    return(model_wide)
  }
  
# Create visualisation plots for model performance
create_performance_plots <- function(model_wide, force_reprocess = FALSE) {
    # Check if outputs already exist
    output_exists <- file.exists(file.path(outdir__04_top, "auc_tss_plot.png"))
    
    if (!force_reprocess && output_exists) {
      message("Performance plots already exist. Skipping visualisation creation.")
      return(NULL)
    }
    
    ## 1. Create the scatterplot "auc_tss_plot" - "Model performance: AUC vs TSS
    auc_tss_plot <- ggplot(model_wide, aes(x = AUC, y = TSS, colour = BG_Method, shape = Threshold)) +
      geom_point(size = 3, alpha = 0.8) +
      facet_wrap(~Species) +
      scale_colour_viridis(discrete = TRUE) +
      labs(
        title = "Model Performance: AUC vs TSS",
        x = "AUC",
        y = "TSS",
        colour = "Background Method",
        shape = "Threshold"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ## 2. Heatmap of model performance by parameter combinations
    # Create configuration names individually to handle potential duplicates
    heatmap_data <- model_wide %>%
      mutate(
        Model_Config = paste(BG_Method, Threshold, Var_Set, Alg_Set, sep = "_")
      ) %>%
      # Add row numbers to ensure unique values in case of duplicates
      mutate(
        Model_Config_Unique = paste0(Model_Config, "_", row_number()),
        # Sort by score
        Model_Config_Unique = factor(Model_Config_Unique, 
                                     levels = Model_Config_Unique[order(Score, decreasing = TRUE)])
      )
    
    # Process for each species
    species_list <- unique(model_wide$Species)
    plot_list <- list()  # Store generated plots
    
    for (species in species_list) {
      safe_species_name <- gsub("[^a-zA-Z0-9]", "_", species)
      
      # Filter for this species
      species_data <- heatmap_data %>%
        filter(Species == species) %>%
        # Ensure it's sorted by Score
        arrange(desc(Score)) %>%
        # Take top 15 models for visualisation
        head(15) %>%
        # Create a display version of Model_Config without the unique identifier
        mutate(Model_Config_Display = gsub("_[0-9]+$", "", Model_Config_Unique))
      
      # Create the heatmap
      heatmap_plot <- ggplot(species_data, aes(y = Model_Config_Display, x = "Combined Score")) +
        geom_tile(aes(fill = Score)) +
        geom_text(aes(label = round(Score, 3)), colour = "white", fontface = "bold") +
        scale_fill_viridis() +
        labs(
          title = paste("Top 15", species, "Models by Combined Score"),
          y = "",
          x = "",
          fill = "Score"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5)
        )
      
      # Save the plot
      plot_name <- paste0("heatmap_plot_", safe_species_name)
      plot_list[[plot_name]] <- heatmap_plot
      ggsave(file.path(outdir__04_top, paste0(plot_name, ".png")), heatmap_plot, width = 8, height = 10)
      
      # Create dot plot showing top models with error bars
      dot_data <- model_wide %>%
        filter(Species == species) %>%
        arrange(desc(Score)) %>%
        head(10) %>%
        mutate(
          Model_Short = paste0(BG_Method, "_", Threshold, "_", Var_Set, "_", Alg_Set),
          Model_Short = factor(Model_Short, levels = rev(Model_Short))
        )
      
      # Convert to long format for plotting
      dot_data_long <- dot_data %>%
        select(Model_Short, AUC, TSS, SD_AUC, SD_TSS, Score) %>%
        pivot_longer(
          cols = c(AUC, TSS),
          names_to = "Metric",
          values_to = "Value"
        ) %>%
        mutate(
          SD = ifelse(Metric == "AUC", SD_AUC, SD_TSS)
        )
      
      # Create the dot plot
      dot_plot <- ggplot(dot_data_long, aes(x = Value, y = Model_Short, colour = Metric)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = Value - SD, xmax = Value + SD), height = 0.2) +
        labs(
          title = paste("Top 10", species, "Models Performance with Standard Deviation"),
          x = "Performance Value",
          y = "",
          colour = "Metric"
        ) +
        scale_colour_manual(values = c("AUC" = "#440154FF", "TSS" = "#21908CFF")) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          legend.position = "bottom"
        )
      
      # Save the dot plot
      plot_name <- paste0("dot_plot_", safe_species_name)
      plot_list[[plot_name]] <- dot_plot
      ggsave(file.path(outdir__04_top, paste0(plot_name, ".png")), dot_plot, width = 10, height = 8)
      
      # Create lollipop chart for ranking visualisation
      lollipop_plot <- ggplot(dot_data, aes(x = Model_Short, y = Score)) +
        geom_segment(aes(x = Model_Short, xend = Model_Short, y = 0, yend = Score),
                     colour = "grey") +
        geom_point(size = 5, colour = "#440154FF") +
        coord_flip() +
        labs(
          title = paste("Top 10", species, "Models Ranking"),
          x = "",
          y = "Combined Score"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8)
        )
      
      # Save the lollipop plot
      plot_name <- paste0("lollipop_", safe_species_name)
      plot_list[[plot_name]] <- lollipop_plot
      ggsave(file.path(outdir__04_top, paste0(plot_name, ".png")), lollipop_plot, width = 8, height = 6)
      
      # Clean up memory after processing each species
      clear_memory()
    }
    
    # Create boxplots for performance by parameter
    # Filter to relevant metrics
    box_data <- model_wide %>%
      select(Species, BG_Method, Threshold, Var_Set, Alg_Set, AUC, TSS) %>%
      pivot_longer(
        cols = c(AUC, TSS),
        names_to = "Metric",
        values_to = "Mean"
      )
    
    # Boxplot by background method
    bg_plot <- ggplot(box_data, aes(x = BG_Method, y = Mean, fill = BG_Method)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(Metric ~ Species, scales = "free_y") +
      scale_fill_viridis_d() +
      labs(
        title = "Model Performance by Background Method",
        x = "Background Method",
        y = "Performance Metric",
        fill = "Background Method"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Boxplot by threshold
    thresh_plot <- ggplot(box_data, aes(x = Threshold, y = Mean, fill = Threshold)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(Metric ~ Species, scales = "free_y") +
      scale_fill_viridis_d() +
      labs(
        title = "Model Performance by Threshold",
        x = "Threshold",
        y = "Performance Metric",
        fill = "Threshold"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Boxplot by variable set
    var_plot <- ggplot(box_data, aes(x = Var_Set, y = Mean, fill = Var_Set)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(Metric ~ Species, scales = "free_y") +
      scale_fill_viridis_d() +
      labs(
        title = "Model Performance by Variable Set",
        x = "Variable Set",
        y = "Performance Metric",
        fill = "Variable Set"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Boxplot by algorithm set
    alg_plot <- ggplot(box_data, aes(x = Alg_Set, y = Mean, fill = Alg_Set)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(Metric ~ Species, scales = "free_y") +
      scale_fill_viridis_d() +
      labs(
        title = "Model Performance by Algorithm Set",
        x = "Algorithm Set",
        y = "Performance Metric",
        fill = "Algorithm Set"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Save the additional plots
    ggsave(file.path(outdir__04_top, "auc_tss_plot.png"), auc_tss_plot, width = 10, height = 6)
    ggsave(file.path(outdir__04_top, "bg_plot.png"), bg_plot, width = 10, height = 6)
    ggsave(file.path(outdir__04_top, "thresh_plot.png"), thresh_plot, width = 10, height = 6)
    ggsave(file.path(outdir__04_top, "var_plot.png"), var_plot, width = 10, height = 6)
    ggsave(file.path(outdir__04_top, "alg_plot.png"), alg_plot, width = 10, height = 6)
    
    # Clean up memory
    clear_memory()
    
    # Return all plots in a list
    return(c(
      list(
        auc_tss_plot = auc_tss_plot,
        bg_plot = bg_plot,
        thresh_plot = thresh_plot,
        var_plot = var_plot,
        alg_plot = alg_plot
      ),
      plot_list
    ))
  }
  
# Create species-specific explained deviance boxplots
create_explained_deviance_plots <- function(results_list, force_reprocess = FALSE) {
    # Check if output already exists
    output_file <- file.path(outdir__04_top, "explained_deviance", "combined_spp_expl_deve_boxplots_comparison.png")
    
    if (!force_reprocess && file.exists(output_file)) {
      message("Explained deviance plot already exists. Skipping plot creation.")
      return(NULL)
    }
    
    # Extract species names
    species_names <- names(results_list)
    
    # Combine results from all species
    combined_data <- do.call(rbind, lapply(species_names, function(sp) {
      data <- results_list[[sp]]$combined_results
      data$species <- sp
      return(data)
    }))
    
    # Convert algorithm names to title case (caps)
    combined_data$algorithm <- toupper(combined_data$algorithm)
    
    # Remove underscores from ensemble names
    combined_data$ensemble_name <- gsub("_", " ", combined_data$ensemble_name)
    
    # Create a faceted plot of explained deviance per species
    main_plot <- ggplot(combined_data, aes(x = algorithm, y = explained_deviance, fill = ensemble_name)) +
      geom_boxplot(lwd = 0.5) +  # Thinner black median line
      facet_wrap(~ species, ncol = 2) +  # Side-by-side panels by species
      scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9"),
                        name = "Top 3 ensembles per species") +  # colourblind-friendly palette with updated legend name
      theme_bw() +
      # Ensure consistent y-axis scale
      scale_y_continuous(limits = c(
        min(combined_data$explained_deviance, na.rm = TRUE),
        max(combined_data$explained_deviance, na.rm = TRUE)
      )) +
      labs(title = "Explained Deviance by Algorithm and Species",
           x = "Algorithm", y = "Explained Deviance") +
      theme(
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(size = 12, face = "italic", hjust = 0.5),  # Italicized species names
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.position = "bottom"
      )
    
    # Filter the data for a specific subset (e.g., first species, RF algorithm)
    inset_species <- species_names[1]
    inset_data <- subset(combined_data, species == inset_species & algorithm == "RF")
    
    # Create inset plot
    inset_plot <- ggplot(inset_data, aes(x = ensemble_name, y = explained_deviance, fill = ensemble_name)) +
      geom_boxplot(lwd = 0.5) +
      scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9")) +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 8),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        # Reduce margins to minimize dead space
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
      ) +
      labs(x = "", y = "")
    
    # Create directory if it doesn't exist
    if (!dir.exists(file.path(outdir__04_top, "explained_deviance"))) {
      dir.create(file.path(outdir__04_top, "explained_deviance"), recursive = TRUE)
    }
    
    # Save both plots together
    output_file_name <- file.path(outdir__04_top, "explained_deviance", "combined_spp_expl_deve_boxplots_comparison.png")
    
    # Save as PNG
    png(output_file_name, width = 12, height = 6, units = "in", res = 300)
    
    # Print the main plot
    print(main_plot)
    
    # Position the inset in the lower right area of the first species panel
    vp <- viewport(x = 0.39, y = 0.30, width = 0.25, height = 0.26)
    print(inset_plot, vp = vp)
    
    # Close the device
    dev.off()
    
    # Clean up memory
    clear_memory()
    
    return(list(main_plot = main_plot, inset_plot = inset_plot))
  }
    
  
#==============================================================================#
#                    ---- 6. Run full workflow ----                           
#==============================================================================#


# Function to run the entire workflow ----
run_sdm_evaluation_workflow <- function(target_alg = "alg4", force_reprocess = FALSE) {
    # Start timing
    start_time <- Sys.time()
    cat("Starting SDM evaluation workflow at", format(start_time), "\n")
    
    # Step 1: Process individual models
    cat("\n=== STEP 1: Processing Individual Models ===\n")
    individual_results <- process_all_models()
    
    # Step 2: Load all SDM models for ensemble evaluation
    cat("\n=== STEP 2: Loading SDM Models for Ensemble Evaluation ===\n")
    calo_models <- load_all_sdm_models(file.path(datafolder03, "calo"))
    vertlept_models <- load_all_sdm_models(file.path(datafolder03, "vertlept"))
    
    # Step 3: Evaluate ensemble models
    cat("\n=== STEP 3: Evaluating Ensemble Models ===\n")
    calo_results <- evaluate_sdm_ensembles(
      calo_models, 
      output_dir = file.path(outdir__04_top, "explained_deviance"),
      species_name = "Calophyllum paniculatum"
    )
    
    vertlept_results <- evaluate_sdm_ensembles(
      vertlept_models, 
      output_dir = file.path(outdir__04_top, "explained_deviance"),
      species_name = "Verticillium/Leptographium"
    )
    
    # Step 4: Rank species ensembles against each other
    cat("\n=== STEP 4: Ranking Species Ensembles ===\n")
    results_list <- list(
      "Calophyllum paniculatum" = calo_results,
      "Verticillium/Leptographium" = vertlept_results
    )
    
    species_rankings <- compare_species_ensembles(
      results_list,
      output_dir = file.path(outdir__04_top, "species_ensemble_rankings")
    )
    
    # Step 5: Process and visualise individual model metrics
    cat("\n=== STEP 5: Processing Model Metrics and Creating visualisations ===\n")
    model_wide <- process_models_metrics()
    
    if (!is.null(model_wide)) {
      performance_plots <- create_performance_plots(model_wide)
    } else {
      cat("Warning: Could not create performance plots due to missing metrics data.\n")
    }
    
    # Step 6: Create explained deviance plots
    cat("\n=== STEP 6: Creating Explained Deviance Plots ===\n")
    explained_deviance_plots <- create_explained_deviance_plots(results_list)
    
    # End timing
    end_time <- Sys.time()
    elapsed <- end_time - start_time
    cat("\nSDM evaluation workflow completed at", format(end_time), "\n")
    cat("Total time elapsed:", format(elapsed), "\n")
    
    # Return all results
    return(list(
      individual_results = individual_results,
      ensemble_results = results_list,
      species_rankings = species_rankings,
      model_metrics = model_wide,
      plots = list(
        performance_plots = if(exists("performance_plots")) performance_plots else NULL,
        explained_deviance_plots = explained_deviance_plots
      )
    ))
  }

  
#==============================================================================#
#                           ---- End of workflow ----
#==============================================================================#