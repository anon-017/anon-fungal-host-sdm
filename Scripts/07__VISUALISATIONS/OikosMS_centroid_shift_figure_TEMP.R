# Extra temporary function to re-run centroid shift radial plots for Oikos MS submission

# ---
# title: "OikosMS_centroid_shift_figure_TEMP.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# date: "2025-12-30"
# update: "2025-12-31"
# note: action to amalgamate with existing functions
# ---

## WHAT
# - Re-extracts and re-runs the range shift analysis to produce the centroid shift radial plot figure
# 
# HOW
# - measures current and figure shift magnitude and direction and converts to compass degrees and km
# 
# WHY
# - Removal of some binary optimisation thresholds (used for initial sensitivity testing and not for final MS submission)

#==============================================================================#
#                 OIKOS MS TEMP CENTROID SHIFT PROCESSING FUNCTION ----
#==============================================================================#

# Updated extract_scenario_info function to handle the actual filename format
extract_scenario_info <- function(raster_path) {
  file_name <- basename(raster_path)
  
  # Initialize all fields
  species <- NA
  threshold <- NA
  scenario <- NA
  time_period <- NA
  is_current <- FALSE
  
  # Check if it's a current projection (contains "projection.tif")
  if (grepl("projection\\.tif", file_name)) {
    is_current <- TRUE
    scenario <- "current"
    time_period <- "current"
    
    # Extract species from current projection filename
    for (sp in c("calo", "vertlept")) {
      if (grepl(sp, file_name)) {
        species <- sp
        break
      }
    }
    
    # Extract threshold from current projection
    if (grepl("opt2", file_name)) threshold <- "opt2"
    else if (grepl("opt4", file_name)) threshold <- "opt4"
    else if (grepl("opt12", file_name)) threshold <- "opt12"
    
  } else {
    # Parse future projection filename: species_scenario_timeperiod_threshold_bin.tif
    parts <- strsplit(file_name, "_")[[1]]
    
    if (length(parts) >= 4) {
      # Extract species (first part)
      if (parts[1] %in% c("calo", "vertlept")) {
        species <- parts[1]
      }
      
      # Extract scenario (second part)
      if (parts[2] %in% c("ssp126", "ssp370", "ssp585")) {
        scenario <- parts[2]
      }
      
      # Extract time period (third part)
      if (parts[3] %in% c("2011-2040", "2041-2070", "2071-2100")) {
        time_period <- parts[3]
      }
      
      # Extract threshold (fourth part)
      if (parts[4] %in% c("opt2", "opt4", "opt12")) {
        threshold <- parts[4]
      }
    }
  }
  
  return(list(
    path = raster_path,
    species = species,
    threshold = threshold,
    scenario = scenario,
    time_period = time_period,
    is_current = is_current
  ))
}

# Updated process function that only processes opt2 for calo and opt4 for vertlept
process_all_combinations <- function() {
  # Initialize empty data frame to store all metrics
  all_metrics <- data.frame()
  
  # IMPORTANT: Define which binary optimisation threshold to use for each species
  species_threshold <- list(
    "calo" = "opt2", # maximise for TSS (sdm package opt 2 - "max(se+sp)")
    "vertlept" = "opt4" # minimum distance to ROC curve (sdm package opt 4 - "minROCdist")
  )
  
  # Get lists of current and future binary rasters
  current_files <- list.files(binary_cur_dir, pattern="\\.tif$", full.names = TRUE)
  future_files <- list.files(binary_fut_dir, pattern="\\.tif$", full.names = TRUE)
  
  cat("Found", length(current_files), "current files and", length(future_files), "future files\n")
  
  # Organize current projections by species (using only the correct threshold)
  current_rasters <- list()
  for (cur_file in current_files) {
    info <- extract_scenario_info(cur_file)
    species <- info$species
    threshold <- info$threshold
    
    # Only load if it's the correct threshold for this species
    if (!is.na(species) && !is.na(threshold) && info$is_current) {
      if (threshold == species_threshold[[species]]) {
        key <- species  # Just use species as key since we have one threshold per species
        current_rasters[[key]] <- terra::rast(cur_file)
        cat("✓ Loaded current projection for", species, "with threshold", threshold, "\n")
      } else {
        cat("  Skipping current", species, "with threshold", threshold, 
            "(using", species_threshold[[species]], "only)\n")
      }
    }
  }
  
  cat("\nProcessing future projections...\n")
  
  # Process each future prediction
  processed_count <- 0
  skipped_count <- 0
  
  for (future_file in future_files) {
    # Extract scenario information
    scenario_info <- extract_scenario_info(future_file)
    
    species <- scenario_info$species
    threshold <- scenario_info$threshold
    
    # Skip if not the correct threshold for this species
    if (!is.na(species) && !is.na(threshold)) {
      if (threshold != species_threshold[[species]]) {
        cat("  Skipping", basename(future_file), 
            "- using", species_threshold[[species]], "for", species, "\n")
        skipped_count <- skipped_count + 1
        next
      }
    }
    
    # Skip if scenario info is incomplete
    if (is.na(species) || 
        is.na(scenario_info$scenario) || 
        is.na(scenario_info$time_period) ||
        is.na(threshold) ||
        scenario_info$is_current ||
        scenario_info$scenario == "current") {
      cat("  Skipping", basename(future_file), "- incomplete info or is current projection\n")
      skipped_count <- skipped_count + 1
      next
    }
    
    # Check if we have the matching current raster
    if (!species %in% names(current_rasters)) {
      cat("  Skipping", basename(future_file), "- no matching current projection available\n")
      skipped_count <- skipped_count + 1
      next
    }
    
    # Get current raster for this species
    current_raster <- current_rasters[[species]]
    
    # Load future prediction raster
    future_raster <- terra::rast(future_file)
    
    cat("✓ Processing:", species, "-", threshold, "-", scenario_info$scenario, 
        "-", scenario_info$time_period, "\n")
    
    # Calculate metrics
    metrics <- analyze_range_shifts(
      species = species,
      current_raster = current_raster,
      future_raster = future_raster,
      threshold = threshold,
      scenario_info = scenario_info,
      dem_raster = elev_MDG
    )
    
    processed_count <- processed_count + 1
    
    # Append to all metrics
    all_metrics <- rbind(all_metrics, metrics)
  }
  
  cat("\nProcessed", processed_count, "combinations, skipped", skipped_count, "\n")
  
  # Save combined metrics
  if (nrow(all_metrics) > 0) {
    all_metrics_file <- file.path(data06_out, "output", "metrics", "extra_settings", 
                                  "summary", "calo_opt2_vertlept_opt4_metrics.csv")
    write.csv(all_metrics, all_metrics_file, row.names = FALSE)
    cat("\nAll metrics saved to:", all_metrics_file, "\n")
    
    # Create summary metrics by species and scenario
    summary_metrics <- aggregate(
      cbind(percent_change, percent_stability, percent_expansion, percent_contraction, 
            centroid_shift_distance, habitat_distance_km, habitat_exposure_km2, spatial_disruption) ~ 
        species + threshold + scenario, 
      data = all_metrics, 
      FUN = mean, 
      na.rm = TRUE
    )
    
    # Save summary metrics
    summary_file <- file.path(data06_out, "output", "metrics", "extra_settings",
                              "summary", "summary_calo_opt2_vertlept_opt4.csv")
    write.csv(summary_metrics, summary_file, row.names = FALSE)
    cat("Summary metrics saved to:", summary_file, "\n")
  } else {
    cat("No metrics were calculated. Check input files.\n")
  }
  
  return(all_metrics)
}

# Radial centroid shift plot for OikosMS

#==============================================================================#
#            RADIAL CENTROID SHIFT PLOT FROM PROCESSED CSV FILES                        #
#==============================================================================#

library(ggplot2)
library(dplyr)

datafolder <- file.path("//xxx")
# Read and filter shift data from individual CSV files - 
read_shift_data <- function(csv_directory, cutoff_date = "2025-12-29") {
  
  # Get all files in the directory
  all_files <- list.files(csv_directory, full.names = TRUE)
  
  # Filter for metric files based on naming pattern
  csv_files <- all_files[grepl("(calo_opt2|vertlept_opt4).*(metrics)", all_files)]
  
  # Filter by modification date to exclude old files
  if (length(csv_files) > 0) {
    file_info <- file.info(csv_files)
    cutoff <- as.Date(cutoff_date)
    recent_files <- csv_files[as.Date(file_info$mtime) >= cutoff]
    csv_files <- recent_files
    
    cat("Found", length(csv_files), "recent metric files (modified after", cutoff_date, ")\n")
  }
  
  if (length(csv_files) > 0) {
    cat("Files to process:\n")
    for(f in csv_files) {
      cat("  ", basename(f), "(", as.Date(file.info(f)$mtime), ")\n")
    }
  }
  
  # Read all files and combine
  all_data <- list()
  files_processed <- 0
  
  for (file in csv_files) {
    # Read the file
    df <- read.csv(file, stringsAsFactors = FALSE)
    files_processed <- files_processed + 1
    
    file_name <- basename(file)
    cat("\nProcessing:", file_name, "\n")
    
    # Check if this is the right threshold combination
    is_calo_opt2 <- grepl("calo_opt2", file_name)
    is_vertlept_opt4 <- grepl("vertlept_opt4", file_name)
    
    if (is_calo_opt2 || is_vertlept_opt4) {
      # Add/fix species and threshold based on filename
      if (is_calo_opt2) {
        df$species <- "calo"
        df$threshold <- "opt2"
      } else if (is_vertlept_opt4) {
        df$species <- "vertlept"
        df$threshold <- "opt4"
      }
      
      # Extract scenario and time period from filename if not present
      if (!"scenario" %in% names(df) || is.na(df$scenario[1])) {
        if (grepl("ssp126", file_name)) df$scenario <- "ssp126"
        else if (grepl("ssp370", file_name)) df$scenario <- "ssp370"
        else if (grepl("ssp585", file_name)) df$scenario <- "ssp585"
      }
      
      if (!"time_period" %in% names(df) || is.na(df$time_period[1])) {
        if (grepl("2011-2040", file_name)) df$time_period <- "2011-2040"
        else if (grepl("2041-2070", file_name)) df$time_period <- "2041-2070"
        else if (grepl("2071-2100", file_name)) df$time_period <- "2071-2100"
      }
      
      # Check if it's a combined file (has all data in one file)
      if (grepl("calo_opt2_vertlept_opt4_metrics", file_name)) {
        cat("  This is a combined metrics file\n")
        # For combined file, ensure species column is correct
        # Don't override species column if it already exists correctly
      }
      
      all_data[[length(all_data) + 1]] <- df
    }
  }
  
  cat("\nCombining", length(all_data), "data frames\n")
  
  # Combine all data
  if (length(all_data) > 0) {
    combined_data <- do.call(rbind, all_data)
    
    # CRITICAL FIX: Remove or correct the species_fullname column if it exists
    if ("species_fullname" %in% names(combined_data)) {
      combined_data$species_fullname <- NULL  # Remove the incorrect column
    }
    
    # Add correct species full names based on species column
    combined_data <- combined_data %>%
      mutate(
        species_fullname = case_when(
          species == "calo" ~ "C. paniculatum",
          species == "vertlept" ~ "L. calophylli",
          TRUE ~ species
        )
      )
    
    # Remove complete duplicates
    combined_data <- combined_data %>% 
      distinct(species, threshold, scenario, time_period, .keep_all = TRUE)
    
    cat("\nTotal unique rows after processing:", nrow(combined_data), "\n")
    
    # Verify we have the right combinations
    cat("\nData breakdown by species:\n")
    combined_data %>%
      group_by(species, species_fullname, threshold) %>%
      summarise(n = n(), .groups = "drop") %>%
      print()
    
    if ("scenario" %in% names(combined_data) && "time_period" %in% names(combined_data)) {
      cat("\nScenarios and time periods:\n")
      combined_data %>%
        group_by(scenario, time_period) %>%
        summarise(n = n(), .groups = "drop") %>%
        arrange(scenario, time_period) %>%
        print()
    }
    
    # Final check - we should have exactly 18 rows
    if (nrow(combined_data) != 18) {
      warning("Expected 18 rows (2 species × 3 scenarios × 3 time periods), but got ", 
              nrow(combined_data))
      
      # Show what combinations we have
      cat("\nDetailed breakdown:\n")
      combined_data %>%
        select(species, threshold, scenario, time_period) %>%
        arrange(species, scenario, time_period) %>%
        print(n = 30)
    }
    
    return(combined_data)
  } else {
    stop("No data found matching the criteria")
  }
}


# Create radial plot
create_radial_plot <- function(shift_df, output_path = NULL) {
  
  # Check that we have a data frame
  if (!is.data.frame(shift_df)) {
    stop("Input must be a data frame, got: ", class(shift_df))
  }
  
  # Check for required column names
  required_cols <- c("species", "threshold", "scenario", "time_period", 
                     "centroid_shift_distance", "centroid_shift_direction_degrees")
  available_cols <- names(shift_df)
  missing_cols <- setdiff(required_cols, available_cols)
  
  if (length(missing_cols) > 0) {
    cat("Warning: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
  }
  
  # Prepare data for plotting
  plot_data <- shift_df %>%
    dplyr::filter(!is.na(centroid_shift_distance)) %>%
    dplyr::mutate(
      # Use the actual column name: centroid_shift_direction_degrees
      direction = centroid_shift_direction_degrees,
      # Convert to compass bearing (0° = North)
      compass_angle = (90 - direction) %% 360,
      # Create clean labels WITHOUT threshold info
      species_label = dplyr::case_when(
        species == "vertlept" ~ "Wilt - L. calophylli",
        species == "calo" ~ "Host - C. paniculatum",
        TRUE ~ species
      ),
      threshold_label = toupper(threshold),  # Keep for backend only
      scenario_upper = toupper(scenario),
      time_label = gsub("-", "–", time_period),
      # CREATE species_threshold FIRST, then convert to factor
      species_threshold = species_label  # Create it from species_label
    ) %>%
    dplyr::mutate(
      # NOW convert to factor with ordered levels
      species_threshold = factor(species_threshold, 
                                 levels = c("Wilt - L. calophylli", "Host - C. paniculatum"))
    )
  
  # Check we have valid data
  if (any(is.na(plot_data$direction))) {
    warning("Some rows have missing direction data")
    plot_data <- plot_data %>% dplyr::filter(!is.na(direction))
  }
  
  cat("\nCreating plot with", nrow(plot_data), "data points\n")
  
  # Create the plot
  p <- ggplot(plot_data) +
    
    # Add centroid shift vectors - THINNER ARROWS
    geom_segment(
      aes(x = compass_angle, xend = compass_angle,
          y = 0, yend = centroid_shift_distance,
          colour = species_threshold,
          linetype = scenario),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      size = 0.7,  # Thinner arrows to make shapes more visible
      alpha = 0.6  # Slightly more opaque
    ) +
    
    # Add endpoint markers - LARGER SHAPES
    geom_point(
      aes(x = compass_angle, y = centroid_shift_distance,
          colour = species_threshold,
          shape = scenario),
      size = 3.5,  # Larger shapes for better visibility
      alpha = 0.7  # More opaque points
    ) +
    
    # Convert to polar coordinates
    coord_polar(theta = "x", start = 0) +
    
    # Compass directions
    scale_x_continuous(
      breaks = c(0, 45, 90, 135, 180, 225, 270, 315),
      labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
      limits = c(0, 360)
    ) +
    
    # Distance scale
    scale_y_continuous(
      labels = function(x) paste0(x, " km")
    ) +
    
    # Colour scheme - CORRECTED to match species_threshold values
    scale_colour_manual(
      values = c(
        "Wilt - L. calophylli" = "purple4",     # Dark purple for vertlept
        "Host - C. paniculatum" = "darkgreen"   # Dark green for calo
      ),
      name = "Species"
    ) +
    
    # Scenario shapes
    scale_shape_manual(
      values = c(
        "ssp126" = 19,  # Filled circle
        "ssp370" = 17,  # Triangle
        "ssp585" = 15   # Square
      ),
      labels = c(
        "ssp126" = "SSP1",
        "ssp370" = "SSP3",
        "ssp585" = "SSP5"
      ),
      name = "Climate scenario"
    ) +
    
    # Scenario line types
    scale_linetype_manual(
      values = c(
        "ssp126" = "solid",
        "ssp370" = "dashed",
        "ssp585" = "dotted"
      ),
      labels = c(
        "ssp126" = "SSP1",
        "ssp370" = "SSP3",
        "ssp585" = "SSP5"
      ),
      name = "Climate scenario"
    ) +
    
    # Facet by time period
    facet_wrap(~ time_label, nrow = 1) +
    
    # Theme - CLEAN WITHOUT BORDERS
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      panel.grid.major = element_line(colour = "grey90", size = 0.3),  # Lighter grid
      panel.grid.minor = element_blank(),  # Remove minor grid
      panel.border = element_blank(),  # No border
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank()  # No background on facet labels
    )
  
  # Save if output path provided
  if (!is.null(output_path)) {
    ggsave(
      filename = output_path,
      plot = p,
      width = 14,
      height = 6,
      dpi = 300,
      bg = "white"
    )
    message("Plot saved to: ", output_path)
  }
  
  return(p)
}

#==============================================================================#
#                           MAIN EXECUTION  ----                                
#==============================================================================#

# Read data with date filter (only files from Dec 29, 2025 onwards)
cat("Reading CSV files from December 29, 2025 onwards\n")
cat("Filtering for:\n")
cat("  - Calophyllum paniculatum (calo) with opt2 threshold\n")
cat("  - Leptographium calophylli (vertlept) with opt4 threshold\n")
cat("Expected: 2 species × 3 pathways × 3 periods = 18 unique combinations\n\n")

shift_data <- read_shift_data(csv_dir, cutoff_date = "2025-12-29")

# Check the final data
cat("\n=== Final Data Check ===\n")
print(head(shift_data))

# Save shift data summary (single file all in one)
write.csv(shift_data, file = file.path(csv_dir, "shifts_all_summary.csv"))

# Centroid visual output
figoutput <- file.path(datafolder, "06__rangeshift_metrics", "output", "OikosMS")

# Write centroid shift radial plot to file
create_radial_plot(shift_data, output_path = file.path(figoutput, "centroid_shifts_radial_opt2_opt4.png"))

#==============================================================================#
#                        ----  End of workflow  ----
#==============================================================================#
