# ---
# title: "04__model_perf_eval.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-02-12"
# update: "2025-12-31"
# ---

# Simplified ensemble SDM model evaluation Workflow
# This script streamlines the evaluation of species distribution models
# for two target species (calo, vertlept) with 16 models per species.

# Libraries
library(sdm)        # For SDM model handling
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping
library(ggplot2)    # For visualisation
library(viridis)    # For colourblind-friendly palettes
library(grid)       # For inset plots
library(plotly)     # For interactive plots (only for 3D plots)
library(htmlwidgets) # For saving web-based visualisations
library(tibble)     # For improved data frames


#==============================================================================#
#                           ---- 1. Workspace set up ----
#==============================================================================#

## Setwd before running source() ----
# First set working directory to "xxx > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

# Define paths
datafolder <- file.path("//xxx")

# Input data folders
datafolder03 <- file.path(datafolder, "03__model_training", "output")

# Output directories
data__04_out <- file.path(datafolder, "04__model_eval_perf", "output")
outdir__04_top <- file.path(data__04_out, "top_performing_models")

# Create output directories
for (dir in c(
  file.path(data__04_out), 
  file.path(outdir__04_top),
  file.path(outdir__04_top, "species_ensemble_rankings"),
  file.path(outdir__04_top, "explained_deviance")
)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Run the main analysis
main <- function() {
  # Define parameters
  force_reprocess <- FALSE  # Set to TRUE to force reprocessing of all models
  create_plots <- FALSE     # Set to TRUE to create PDP and variable importance plots
  
  # Run the workflow
  cat("Starting SDM evaluation workflow...\n")
  run_sdm_evaluation(force_reprocess = force_reprocess, create_plots = create_plots)
  
  cat("Workflow completed!\n")
}

# Uncomment to run the workflow
# main()

# Define the species list
species_list <- c("calo", "vertlept")
species_display_names <- c(
  "calo" = "Calophyllum paniculatum",
  "vertlept" = "Verticillium/Leptographium"
)

# Define model algorithm set - only using alg3 now
alg_set <- c("rf", "brt", "maxent")

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

# Function to clean and format variable labels
clean_labels <- function(text) {
  clean_text <- character(length(text))
  
  for(i in seq_along(text)) {
    if(text[i] %in% names(bio_lookup)) {
      clean_text[i] <- bio_lookup[text[i]]
    } else {
      clean_text[i] <- gsub("[_.]", " ", text[i])
    }
  }
  
  return(clean_text)
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
load_all_sdm_models <- function(species_directory, species_name = NULL) {
  # Get the path to the directory
  algsets_dir <- file.path(species_directory)
  
  # If species_name not provided, try to extract from directory path
  if(is.null(species_name)) {
    species_name <- basename(species_directory)
  }
  
  # List all .sdm files (non-recursive - only get top-level files)
  sdm_files <- list.files(algsets_dir, pattern = "\\.sdm$", 
                          recursive = FALSE, full.names = TRUE)
  
  if(length(sdm_files) == 0) {
    warning(paste("No SDM files found in directory:", species_directory))
    return(NULL)
  }
  
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
    
    # Remove the model variable to avoid conflicts
    rm(sdm_model)
    
    # Clean up memory after each model
    clear_memory()
  }
  
  # Return the list of models
  return(sdm_models_list)
}

# Function to calculate null deviance (for explained deviance)
calculate_null_deviance <- function(sdm_model, response = "occ", algorithm, model_id) {
  # Extract the specific model from the ensemble
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

#==============================================================================#
#                 ---- 3. Individual Model Processing ----                     
#==============================================================================#

# Process a single model to extract metrics, create plots, etc.
process_one_model <- function(model_path, create_plots = FALSE) {
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
    Mean = sapply(names(eval_metrics_desc), function(m) {
      val <- eval_metrics[[m]]
      if(!is.numeric(val)) {
        warning(paste("Non-numeric value found for metric:", m))
        return(NA)
      }
      mean(val, na.rm = TRUE)
    }),
    SD = sapply(names(eval_metrics_desc), function(m) {
      val <- eval_metrics[[m]]
      if(!is.numeric(val)) return(NA)
      sd(val, na.rm = TRUE)
    }),
    Description = unname(eval_metrics_desc),
    Species = model_components$species,
    BG_Method = model_components$bg_method,
    Threshold = model_components$threshold,
    Var_Set = model_components$var_set,
    Alg_Set = model_components$alg_set
  )
  
  # 2. Calculate explained deviance
  explained_deviance <- NA
  individual_deviances <- list()  # To store individual algorithm+model results
  
  tryCatch({
    all_devs <- c()
    
    # Get algorithm types
    algos <- names(sdm_model@models$occ)
    
    # Process each algorithm
    for (algo in algos) {
      # Initialize list for this algorithm
      individual_deviances[[algo]] <- list()
      
      # Get model IDs for this algorithm
      model_ids <- names(sdm_model@models[["occ"]][[algo]])
      
      # Process each model ID
      for (mid in model_ids) {
        # Calculate null deviance
        null_dev <- calculate_null_deviance(sdm_model, "occ", algo, mid)
        
        # Get residual deviance
        model_data <- sdm_model@models[["occ"]][[algo]][[mid]]
        residual_dev <- model_data@evaluation$training@statistics$Deviance
        
        # Calculate explained deviance
        explained_dev <- 1 - (residual_dev / null_dev)
        
        # Store individual result
        individual_deviances[[algo]][[mid]] <- explained_dev
        
        # Add to collection for mean calculation
        all_devs <- c(all_devs, explained_dev)
      }
    }
    
    # Calculate mean explained deviance across all models
    explained_deviance <- mean(all_devs, na.rm = TRUE)
    
  }, error = function(e) {
    warning(paste("Error calculating explained deviance:", e$message))
  })
  
  # Add individual values and ensemble mean to metrics_summary
  if(!is.null(metrics_summary)) {
    metrics_summary$ExplainedDeviance <- explained_deviance
    metrics_summary$IndividualDeviances <- individual_deviances
  }
  
  # Return both values (optional, if you need them outside metrics_summary)
  results <- list(
    ensemble_mean = explained_deviance,
    individual_values = individual_deviances
  )
  
  # 3. Process and save variable importance and response curves only if requested
  if(create_plots) {
    # Get variable importance
    var_imp <- getVarImp(sdm_model)
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
    
    # Get response curves
    model_ids <- getModelId(sdm_model)
    response_curves <- getResponseCurve(sdm_model, id = model_ids)
    
    # Variable Importance Plot
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
    
    # Response Curves (2D)
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
    
    # 3D Response Plots (commented out)
    # for(i in 1:(length(response_curves@variables)-1)) {
    #   for(j in (i+1):length(response_curves@variables)) {
    #     var1 <- response_curves@variables[i]
    #     var2 <- response_curves@variables[j]
    #     
    #     tryCatch({
    #       fig <- create_3d_pdp(response_curves, var1, var2)
    #       
    #       output_file <- file.path(data__04_out, model_components$species, "plots", "pdp_3d",
    #                             paste0(tools::file_path_sans_ext(model_name), "_", var1, "_vs_", var2, "_3d.html"))
    #       
    #       # Save the plotly figure
    #       saveWidget(fig, output_file)
    #     }, error = function(e) {
    #       warning(paste("Error creating 3D plot for", var1, "vs", var2, ":", e$message))
    #     })
    #   }
    # }
  }
  
  # 4. Save results
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
process_all_models <- function(force_reprocess = FALSE, create_plots = FALSE) {
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
    
    # Filter to get only 16 models per species
    if(length(model_files) > 16) {
      cat("Found", length(model_files), "models, limiting to 16 as requested\n")
      # Get only models with alg3
      model_files <- grep("alg3", model_files, value = TRUE)
      # If still more than 16, take the first 16
      if(length(model_files) > 16) {
        model_files <- model_files[1:16]
      }
    }
    
    cat("Processing", length(model_files), "models for", species, "\n")
    
    # Process each model
    species_results <- list()
    for(model_path in model_files) {
      model_name <- basename(model_path)
      cat("  Processing:", model_name, "\n")
      
      # Skip if already processed and not forcing reprocessing
      output_file <- file.path(species_dir, "metrics", 
                               paste0(tools::file_path_sans_ext(model_name), "_results.rds"))
      
      if(file.exists(output_file) && !force_reprocess) {
        cat("  Already processed, loading from file\n")
        result <- readRDS(output_file)$metrics
      } else {
        tryCatch({
          result <- process_one_model(model_path, create_plots = create_plots)
        }, error = function(e) {
          warning(paste("Error processing", model_name, ":", e$message))
          return(NULL)
        })
      }
      
      if(!is.null(result)) {
        species_results[[model_name]] <- result
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
        facet_grid(. ~ Threshold) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste(species, "Model Performance Comparison"),
             x = "Habitat Type",
             y = "AUC Score",
             colour = "Background Method")
      
      ggsave(file.path(species_dir, "plots", paste0(species, "_performance_comparison.png")),
             p_comparison, dpi = 300, width = 10, height = 6)
      
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
    
    # Create wide format of metrics for easier analysis
    # First filter to key metrics (AUC, TSS, ExplainedDeviance)
    key_metrics <- all_combined %>%
      filter(Metric %in% c("AUC", "TSS")) %>%
      select(Model, Metric, Mean, SD, Species, BG_Method, Threshold, Var_Set, Alg_Set)
    
    # Convert to wide format
    metrics_wide <- key_metrics %>%
      pivot_wider(
        id_cols = c(Model, Species, BG_Method, Threshold, Var_Set, Alg_Set),
        names_from = Metric,
        values_from = c(Mean, SD)
      ) %>%
      # Rename columns to cleaner format
      rename(
        AUC = Mean_AUC,
        TSS = Mean_TSS,
        AUC_SD = SD_AUC,
        TSS_SD = SD_TSS
      )
    
    # Add combined score (weighted AUC and TSS)
    metrics_wide$CombinedScore <- 0.6 * metrics_wide$AUC + 0.4 * metrics_wide$TSS
    
    # Add ExplainedDeviance if available
    if("ExplainedDeviance" %in% all_combined$Metric) {
      ed_metrics <- all_combined %>%
        filter(Metric == "ExplainedDeviance") %>%
        select(Model, Mean, SD)
      
      metrics_wide <- metrics_wide %>%
        left_join(ed_metrics, by = "Model") %>%
        rename(
          ExplainedDeviance = Mean,
          ExplainedDeviance_SD = SD
        )
    }
    
    # Save wide format metrics
    write.csv(metrics_wide, file.path(data__04_out, "all_species_metrics_wide.csv"), row.names = FALSE)
    
    # Create overall performance comparison
    p_overall <- ggplot(metrics_wide, 
                        aes(x = AUC, y = TSS, colour = Species, shape = BG_Method)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_colour_viridis_d() +
      theme_minimal() +
      labs(title = "Overall Model Performance Comparison",
           x = "AUC Score",
           y = "TSS Score",
           colour = "Species",
           shape = "Background Method")
    
    ggsave(file.path(data__04_out, "overall_performance_comparison.png"),
           p_overall, dpi = 300, width = 10, height = 8)
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
  
  # Convert list of results to a data frame
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
  
  # Sort by algorithm and model_id
  results_df <- results_df %>%
    arrange(algorithm, model_id)
  
  return(results_df)
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
  all_ensemble_summaries <- list()
  
  # Process each ensemble model
  for (i in seq_along(ensemble_models)) {
    # Get ensemble name
    if (!is.null(names(ensemble_models)) && names(ensemble_models)[i] != "") {
      ensemble_name <- names(ensemble_models)[i]
    } else {
      ensemble_name <- paste0("ensemble_", i)
    }
    
    # Check if this ensemble has already been processed
    results_path <- file.path(output_dir, paste0(safe_species_name, "_", ensemble_name, ".csv"))
    if (!force_reprocess && !is.null(output_dir) && file.exists(results_path)) {
      message("Ensemble ", ensemble_name, " for species ", species_name, " already processed. Loading results.")
      # Load existing results
      eval_table <- read.csv(results_path)
      all_ensemble_results[[ensemble_name]] <- eval_table
      
      # Calculate summary statistics
      summary_stats <- eval_table %>%
        group_by(algorithm) %>%
        summarize(
          n_models = n(),
          AUC_mean = mean(AUC, na.rm = TRUE),
          AUC_sd = sd(AUC, na.rm = TRUE),
          TSS_mean = mean(TSS, na.rm = TRUE),
          TSS_sd = sd(TSS, na.rm = TRUE),
          explained_deviance_mean = mean(explained_deviance, na.rm = TRUE),
          explained_deviance_sd = sd(explained_deviance, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          species = species_name,
          ensemble_name = ensemble_name
        )
      
      all_ensemble_summaries[[ensemble_name]] <- summary_stats
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
    
    # Calculate summary statistics
    summary_stats <- eval_table %>%
      group_by(algorithm) %>%
      summarize(
        n_models = n(),
        AUC_mean = mean(AUC, na.rm = TRUE),
        AUC_sd = sd(AUC, na.rm = TRUE),
        TSS_mean = mean(TSS, na.rm = TRUE),
        TSS_sd = sd(TSS, na.rm = TRUE),
        explained_deviance_mean = mean(explained_deviance, na.rm = TRUE),
        explained_deviance_sd = sd(explained_deviance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        species = species_name,
        ensemble_name = ensemble_name
      )
    
    # Store in results list
    all_ensemble_results[[ensemble_name]] <- eval_table
    all_ensemble_summaries[[ensemble_name]] <- summary_stats
    
    # Save individual results if output directory is specified
    if (!is.null(output_dir)) {
      write.csv(eval_table, results_path, row.names = FALSE)
      
      # Save summary stats
      summary_path <- file.path(output_dir, paste0(safe_species_name, "_", ensemble_name, "_summary.csv"))
      write.csv(summary_stats, summary_path, row.names = FALSE)
      
      message("Results saved to: ", results_path)
    }
  }
  
  # Combine summary statistics across all ensembles
  if (length(all_ensemble_summaries) > 0) {
    combined_summary <- do.call(rbind, all_ensemble_summaries)
    
    # Save combined summary if output directory is specified
    if (!is.null(output_dir)) {
      write.csv(combined_summary, 
                file.path(output_dir, paste0(safe_species_name, "_all_ensembles_summary.csv")), 
                row.names = FALSE)
    }
  } else {
    combined_summary <- NULL
  }
  
  return(list(
    results = all_ensemble_results,
    summary_stats = combined_summary
  ))
}

# Function to compare and rank species ensembles
rank_species_ensembles <- function(species_results, output_dir = NULL,
                                   auc_weight = 0.6, tss_weight = 0.4, ed_weight = 0.0,
                                   sd_penalty = 0.5, force_reprocess = FALSE) {
  
  # Validate inputs
  if (!is.list(species_results) || length(species_results) == 0) {
    stop("species_results must be a non-empty list where each element contains evaluation results")
  }
  
  # Check if output already exists
  if (!force_reprocess && !is.null(output_dir)) {
    ranking_file <- file.path(output_dir, "species_ensemble_ranking.csv")
    if (file.exists(ranking_file)) {
      message("Species ensemble ranking already exists. Loading from file.")
      cross_species_ranking <- read.csv(ranking_file)
      return(cross_species_ranking)
    }
  }
  
  # Normalize weights to sum to 1
  total_weight <- auc_weight + tss_weight + ed_weight
  if (total_weight != 1) {
    auc_weight <- auc_weight / total_weight
    tss_weight <- tss_weight / total_weight
    ed_weight <- ed_weight / total_weight
  }
  
  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract summary statistics from all species
  species_names <- names(species_results)
  all_summaries <- list()
  
  for (species_name in species_names) {
    # Get summary stats for this species
    if (is.null(species_results[[species_name]]$summary_stats)) {
      warning("No summary statistics found for species: ", species_name)
      next
    }
    
    # Add to list
    all_summaries[[species_name]] <- species_results[[species_name]]$summary_stats
  }
  
  # Combine summaries
  if (length(all_summaries) == 0) {
    stop("No valid summary statistics found for any species")
  }
  
  combined_summaries <- do.call(rbind, all_summaries)
  
  # Calculate ranking score
  ranking <- combined_summaries %>%
    mutate(
      # Calculate weighted score
      Score = auc_weight * AUC_mean + 
        tss_weight * TSS_mean + 
        ed_weight * explained_deviance_mean -
        sd_penalty * (auc_weight * AUC_sd + 
                        tss_weight * TSS_sd + 
                        ed_weight * explained_deviance_sd),
      
      # Create species_ensemble identifier
      species_ensemble = paste(species, ensemble_name, sep = "_")
    ) %>%
    # Rank all ensembles across species
    arrange(desc(Score)) %>%
    mutate(rank = row_number())
  
  # Save ranking if output directory provided
  if (!is.null(output_dir)) {
    write.csv(ranking, file.path(output_dir, "species_ensemble_ranking.csv"), row.names = FALSE)
    
    # Create summary by species
    species_summary <- ranking %>%
      group_by(species) %>%
      summarize(
        best_ensemble = first(ensemble_name),
        best_score = first(Score),
        mean_score = mean(Score),
        min_score = min(Score),
        max_score = max(Score),
        score_range = max_score - min_score,
        ensemble_count = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(best_score))
    
    write.csv(species_summary, file.path(output_dir, "species_performance_summary.csv"), row.names = FALSE)
    
    # Create visualisation of top models
    top_n_models <- min(10, nrow(ranking))
    
    p1 <- ranking %>%
      head(top_n_models) %>%
      ggplot(aes(x = reorder(species_ensemble, Score), y = Score, fill = species)) +
      geom_col() +
      coord_flip() +
      scale_fill_viridis_d() +
      labs(
        title = paste("Top", top_n_models, "Species Ensembles by Performance Score"),
        x = "Species_Ensemble",
        y = "Performance Score"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "top_ensembles_score.png"), p1, width = 10, height = 8)
    
    # Create AUC vs TSS scatterplot
    p2 <- combined_summaries %>%
      ggplot(aes(x = AUC_mean, y = TSS_mean, colour = species)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = ensemble_name), hjust = -0.2, vjust = 0.5, size = 3) +
      scale_colour_viridis_d() +
      labs(
        title = "Ensemble Performance: AUC vs TSS",
        x = "AUC (Mean)",
        y = "TSS (Mean)",
        colour = "Species"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "auc_vs_tss_ensembles.png"), p2, width = 10, height = 8)
  }
  
  return(ranking)
}

#==============================================================================#
#                 ---- 5. Visualisation and outputs ----                      
#==============================================================================#

# Create explained deviance boxplots for visual comparison
create_explained_deviance_plots <- function(results_list, output_dir = NULL, force_reprocess = FALSE) {
  # Check if output already exists
  if (!is.null(output_dir)) {
    output_file <- file.path(output_dir, "combined_explained_deviance_boxplots.png")
    
    if (!force_reprocess && file.exists(output_file)) {
      message("Explained deviance plot already exists. Skipping creation.")
      return(NULL)
    }
    
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Extract species names
  species_names <- names(results_list)
  
  # Extract data for plotting
  plot_data <- list()
  
  for (species in species_names) {
    # Get results for this species
    species_results <- results_list[[species]]
    
    # Extract data from each ensemble
    ensembles <- names(species_results$results)
    
    # Combine data from all ensembles for this species
    species_data <- list()
    
    for (ens in ensembles) {
      # Extract data
      ens_data <- species_results$results[[ens]] %>%
        select(algorithm, explained_deviance, ensemble_name, species)
      
      species_data[[ens]] <- ens_data
    }
    
    # Combine all data for this species
    if (length(species_data) > 0) {
      plot_data[[species]] <- do.call(rbind, species_data)
    }
  }
  
  # Combine all species data
  if (length(plot_data) == 0) {
    warning("No valid data found for creating explained deviance plots")
    return(NULL)
  }
  
  combined_data <- do.call(rbind, plot_data)
  
  # Convert algorithm names to uppercase for consistency
  combined_data$algorithm <- toupper(combined_data$algorithm)
  
  # Create boxplot
  p <- ggplot(combined_data, aes(x = algorithm, y = explained_deviance, fill = ensemble_name)) +
    geom_boxplot(lwd = 0.5) +
    facet_wrap(~ species, ncol = 2) +
    scale_fill_viridis_d(name = "Ensemble") +
    theme_bw() +
    labs(
      title = "Explained Deviance by Algorithm and Species",
      x = "Algorithm",
      y = "Explained Deviance"
    ) +
    theme(
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(size = 12, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 11),
      legend.title = element_text(size = 10),
      legend.position = "bottom"
    )
  
  # Save plot if output directory specified
  if (!is.null(output_dir)) {
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    message("Explained deviance plot saved to: ", output_file)
  }
  
  return(p)
}

# Create publication-ready tables and figures
create_publication_tables <- function(individual_metrics, ensemble_results, output_dir = NULL) {
  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Main text table: Top 3 models per species with key metrics
  if (!is.null(individual_metrics)) {
    # Ensure metrics are in wide format
    if ("Metric" %in% names(individual_metrics)) {
      # Need to convert from long to wide format
      metrics_wide <- individual_metrics %>%
        filter(Metric %in% c("AUC", "TSS")) %>%
        select(Model, Metric, Mean, SD, Species, BG_Method, Threshold, Var_Set, Alg_Set) %>%
        pivot_wider(
          id_cols = c(Model, Species, BG_Method, Threshold, Var_Set, Alg_Set),
          names_from = Metric,
          values_from = c(Mean, SD)
        ) %>%
        rename(
          AUC = Mean_AUC,
          TSS = Mean_TSS,
          AUC_SD = SD_AUC,
          TSS_SD = SD_TSS
        )
    } else {
      metrics_wide <- individual_metrics
    }
    
    # Add combined score if it doesn't exist
    if (!"CombinedScore" %in% names(metrics_wide)) {
      metrics_wide$CombinedScore <- 0.6 * metrics_wide$AUC + 0.4 * metrics_wide$TSS
    }
    
    # Get top 3 models per species
    top_models <- metrics_wide %>%
      group_by(Species) %>%
      arrange(desc(CombinedScore)) %>%
      slice_head(n = 3) %>%
      ungroup()
    
    # Format for publication - main table
    main_table <- top_models %>%
      select(Species, Model, BG_Method, Threshold, Var_Set, AUC, TSS, CombinedScore) %>%
      arrange(Species, desc(CombinedScore)) %>%
      mutate(
        # Round numeric columns
        AUC = round(AUC, 3),
        TSS = round(TSS, 3),
        CombinedScore = round(CombinedScore, 3),
        # Clean up species names for display
        Species = case_when(
          Species == "calo" ~ "Calophyllum paniculatum",
          Species == "vertlept" ~ "Verticillium/Leptographium",
          TRUE ~ Species
        ),
        # Clean up other fields
        BG_Method = case_when(
          BG_Method == "ecor" ~ "Ecoregion",
          BG_Method == "whole" ~ "Whole",
          TRUE ~ BG_Method
        ),
        Threshold = gsub("thresh_", "", Threshold),
        Var_Set = case_when(
          Var_Set == "forest_aspect" ~ "Forest + Aspect",
          Var_Set == "aspect" ~ "Aspect",
          Var_Set == "base" ~ "Base",
          TRUE ~ Var_Set
        )
      )
    
    # Save main table
    if (!is.null(output_dir)) {
      write.csv(main_table, file.path(output_dir, "top_models_main_table.csv"), row.names = FALSE)
    }
    
    # Create appendix table with all models
    appendix_table <- metrics_wide %>%
      select(Species, Model, BG_Method, Threshold, Var_Set, AUC, AUC_SD, TSS, TSS_SD, CombinedScore) %>%
      arrange(Species, desc(CombinedScore)) %>%
      mutate(
        # Round numeric columns
        AUC = round(AUC, 3),
        AUC_SD = round(AUC_SD, 3),
        TSS = round(TSS, 3),
        TSS_SD = round(TSS_SD, 3),
        CombinedScore = round(CombinedScore, 3),
        # Clean up species names for display
        Species = case_when(
          Species == "calo" ~ "Calophyllum paniculatum",
          Species == "vertlept" ~ "Verticillium/Leptographium",
          TRUE ~ Species
        ),
        # Clean up other fields
        BG_Method = case_when(
          BG_Method == "ecor" ~ "Ecoregion",
          BG_Method == "whole" ~ "Whole",
          TRUE ~ BG_Method
        ),
        Threshold = gsub("thresh_", "", Threshold),
        Var_Set = case_when(
          Var_Set == "forest_aspect" ~ "Forest + Aspect",
          Var_Set == "aspect" ~ "Aspect",
          Var_Set == "base" ~ "Base",
          TRUE ~ Var_Set
        )
      )
    
    # Save appendix table
    if (!is.null(output_dir)) {
      write.csv(appendix_table, file.path(output_dir, "all_models_appendix_table.csv"), row.names = FALSE)
    }
  }
  
  # 2. Ensemble performance table
  if (!is.null(ensemble_results)) {
    # Extract species names
    species_names <- names(ensemble_results)
    
    # Extract summary stats
    ensemble_summaries <- list()
    for (species in species_names) {
      if (!is.null(ensemble_results[[species]]$summary_stats)) {
        ensemble_summaries[[species]] <- ensemble_results[[species]]$summary_stats
      }
    }
    
    # Combine summaries if available
    if (length(ensemble_summaries) > 0) {
      combined_ensemble_summary <- do.call(rbind, ensemble_summaries)
      
      # Format for publication
      ensemble_table <- combined_ensemble_summary %>%
        select(species, ensemble_name, algorithm, AUC_mean, AUC_sd, TSS_mean, TSS_sd, 
               explained_deviance_mean, explained_deviance_sd) %>%
        arrange(species, ensemble_name, algorithm) %>%
        mutate(
          # Round numeric columns
          AUC_mean = round(AUC_mean, 3),
          AUC_sd = round(AUC_sd, 3),
          TSS_mean = round(TSS_mean, 3),
          TSS_sd = round(TSS_sd, 3),
          explained_deviance_mean = round(explained_deviance_mean, 3),
          explained_deviance_sd = round(explained_deviance_sd, 3),
          
          # Clean up species names
          species = case_when(
            species == "calo" ~ "Calophyllum paniculatum",
            species == "vertlept" ~ "Verticillium/Leptographium",
            TRUE ~ species
          ),
          
          # Uppercase algorithm
          algorithm = toupper(algorithm)
        )
      
      # Calculate average across algorithms for each ensemble
      ensemble_avg <- ensemble_table %>%
        group_by(species, ensemble_name) %>%
        summarize(
          AUC_overall = mean(AUC_mean, na.rm = TRUE),
          AUC_sd_overall = mean(AUC_sd, na.rm = TRUE),
          TSS_overall = mean(TSS_mean, na.rm = TRUE),
          TSS_sd_overall = mean(TSS_sd, na.rm = TRUE),
          ExpDeviance_overall = mean(explained_deviance_mean, na.rm = TRUE),
          ExpDeviance_sd_overall = mean(explained_deviance_sd, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(species, desc(AUC_overall))
      
      # Save ensemble tables
      if (!is.null(output_dir)) {
        write.csv(ensemble_table, file.path(output_dir, "ensemble_performance_detailed.csv"), row.names = FALSE)
        write.csv(ensemble_avg, file.path(output_dir, "ensemble_performance_summary.csv"), row.names = FALSE)
      }
    }
  }
  
  # Return the created tables
  return(list(
    main_table = if(exists("main_table")) main_table else NULL,
    appendix_table = if(exists("appendix_table")) appendix_table else NULL,
    ensemble_table = if(exists("ensemble_table")) ensemble_table else NULL,
    ensemble_avg = if(exists("ensemble_avg")) ensemble_avg else NULL
  ))
}



#==============================================================================#
#                    ---- 6. Run full workflow ----                           
#==============================================================================#

# Run the entire workflow
run_sdm_evaluation <- function(force_reprocess = FALSE, create_plots = FALSE) {
  # Start timing
  start_time <- Sys.time()
  cat("Starting SDM evaluation workflow at", format(start_time), "\n")
  
  # Step 1: Process individual models
  cat("\n=== STEP 1: Processing Individual Models ===\n")
  individual_results <- process_all_models(force_reprocess = force_reprocess, create_plots = create_plots)
  
  # Step 2: Load all SDM models for ensemble evaluation
  cat("\n=== STEP 2: Loading SDM Models for Ensemble Evaluation ===\n")
  sdm_models <- list()
  
  for (species in species_list) {
    cat("Loading models for", species, "\n")
    sdm_models[[species]] <- load_all_sdm_models(file.path(datafolder03, species))
  }
  
  # Step 3: Evaluate ensemble models
  cat("\n=== STEP 3: Evaluating Ensemble Models ===\n")
  ensemble_results <- list()
  
  for (species in species_list) {
    cat("Evaluating ensembles for", species, "\n")
    ensemble_results[[species]] <- evaluate_sdm_ensembles(
      sdm_models[[species]], 
      output_dir = file.path(outdir__04_top, "ensemble_evaluation"),
      species_name = species,
      force_reprocess = force_reprocess
    )
  }
  
  # Step 4: Rank species ensembles against each other
  cat("\n=== STEP 4: Ranking Species Ensembles ===\n")
  species_rankings <- rank_species_ensembles(
    ensemble_results,
    output_dir = file.path(outdir__04_top, "species_ensemble_rankings"),
    auc_weight = 0.6,
    tss_weight = 0.3,
    ed_weight = 0.1,
    sd_penalty = 0.5,
    force_reprocess = force_reprocess
  )
  
  # Step 5: Create explained deviance visualisation
  cat("\n=== STEP 5: Creating Explained Deviance visualisation ===\n")
  ed_plots <- create_explained_deviance_plots(
    ensemble_results,
    output_dir = file.path(outdir__04_top, "explained_deviance"),
    force_reprocess = force_reprocess
  )
  
  # Step 6: Create publication tables
  cat("\n=== STEP 6: Creating Publication Tables ===\n")
  # First load wide metrics if they exist
  metrics_wide_path <- file.path(data__04_out, "all_species_metrics_wide.csv")
  metrics_wide <- NULL
  
  if (file.exists(metrics_wide_path)) {
    metrics_wide <- read.csv(metrics_wide_path)
  }
  
  pub_tables <- create_publication_tables(
    individual_metrics = metrics_wide,
    ensemble_results = ensemble_results,
    output_dir = file.path(outdir__04_top, "publication_tables")
  )
  
  # End timing
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  cat("\nSDM evaluation workflow completed at", format(end_time), "\n")
  cat("Total time elapsed:", format(elapsed), "\n")
  
  # Return all results
  return(list(
    individual_results = individual_results,
    ensemble_results = ensemble_results,
    species_rankings = species_rankings,
    publication_tables = pub_tables,
    elapsed_time = elapsed
  ))
}


#==============================================================================#
#                        ----  Run workflow  ----
#==============================================================================#


perf_eval <- run_sdm_evaluation()


#==============================================================================#
#                        ----  End of workflow  ----
#==============================================================================#

# # # Clean up
# gc()
# rm(list = ls())