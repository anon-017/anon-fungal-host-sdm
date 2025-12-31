# ---
# title: "00__species_occurrence_cleaning_bg.R"
# manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"
# corresponding_author: "xxx"
# coauthors: "xxx, K., xxx, A., xxx, M., xxx, N., xxx, R."
# update: "2025-12-31"
# notes:
# ---

## WHAT
### - This set of code downloads / loads in species occurrence data from various sources, then crops it down
### - Focal species: Calophyllum paniculatum; & wilt Genera: Verticillium; Leptographium
### - Removes duplicates from presence locations
### - Creates background points via  species-specific methods (2 degrees far, ecoregion) and whole area sampling (minus 1.5km buffer around presences)
### - Removes duplicates background points 
### - (OPTIONAL) Joins presence and background points into a single file calo_occs and vert_occs ready to be combined with bioclimatic data in next script (01__)

## HOW
### - Removes duplicate presences from within 1km cells
### - Calo: intersects presence points with Ecoregions and samples background from within the ecoregion area
### - vertlept: 2 degrees far (buffers to equivalent of 2deg from presence points), erases buffer area from study area
###   and samples background points from within area > 2 degrees away from presences.

## WHY
### - Species occurrences require cleaning and spatial/survey biases removed
### - Background points to support model building due to presence-only species data
### - Presence-background combined used to extract environmental data in order to run collinearity analysis on environmental predictors.

## WHERE
### - input = AOIs (study extents), ecoregions, elevation, GBIF, and field survey
### - intermediate = cropped data and vector versions of study area e.g. vtemplateA, ecoregsA)
### - output: out_dir = .\00_species_occurrence_cleaning > output
    #### pres_out_dir = out_dir > presence
    #### dupfree_out_dir  = pres_out_dir > dupfree (duplicates removed from presence points)
    #### bg_out_dir = out_dir > background (background data produced (currently multiple methods))
    #### plot_out_dir = out_dir > plots (plots produced that will become intro figures)
    #### occ_out_dir = out_dir > combined_occs (final combined, duplicate free, presence-background occurrence dataframes)


#==============================================================================#
#                          ---- 0. Workspace set up ----
#==============================================================================#


#### Setwd before running source() ----
# First set working directory to "Calo_SDM_Ensemble > Scripts" so source() work to load functions
setwd("~/GitHub/anon-fungal-host-sdms/Scripts")

#### Load functions.R ----
# Set up project environment and load packages and functions
source("./functions.R")

# Remove functions not relating to this script
rm(install.load.package,package_vec, process.climate.data, repair.tiff, swap_coords, thin, validate.tiff, align.forest.raster, align.raster, validate.processed.climate.files)

#### Set up base folders ----
# State main base folder for this script where everything else will become relative paths using "file.path()"
#datafolder__00 <- file.path("//xxx/00__species_occurrence_cleaning")
datafolder__00 <- file.path("//xxx/00__species_occurrence_cleaning") # vpn path

# AOIs of template extents
aois <- file.path(datafolder__00, "input", "AOI")

# Collaborator surveys with occurrence records
#surveys <- file.path(datafolder__00, "input", "surveys")

#### Create output folder ----
# Create output sub-folder
# Change date of folder if ever doing a new download to not overwrite
out_dir <- file.path(datafolder__00, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

# Create plots output sub-folder
plot_out_dir <- file.path(out_dir, "plots")
dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist


#==============================================================================#
#                    1. Load study area templates (AFR / MDG) ----
#==============================================================================#

# Template CRS (African Equal Albers / WGS 1984 Albers for Africa)
crs <- "ESRI:102022"

# Load ADMIN BOUNDARYS as a spatial object - already pre-processed in QGIS (projected to crs, resampled to 1000)

# Just African continent no islands, no MDG
templateB <- rast(file.path(aois, "AFR", "AFR_aea_template1km.tif"))

# Just MDG
templateA <- rast(file.path(aois,"MDG", "MDG_aea.tif"))

# Both extents together
templateAB <- rast(file.path(aois, "AFR", "AFR_MDG_aea.tif"))


#==============================================================================#
#                       2. GBIF - Focal Plant (CALO)  ----
#==============================================================================#

#### ***User data required for rgbif occ_download() ----
# # Sensitive info required for GBIF
# user <- "xxx"
# pwdd <- "xxx"
# email <- "xxx.xxx@xxx.xx.xx"

##### Calophyllum paniculatum ----

##### Create json request ---- 
# Create pred key for json format/talking directly with gbif (to be able to download a specific citation for writing up)
# Taxon key taken from gbif.org before running rgbif download (Calophyllum paniculatum taxonKey: 8006459)
calo_download <- occ_download(
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("taxonKey", 8006459),
  format = "SIMPLE_CSV",
  user = user,
  pwd = pwdd,
  email = email
)

##### Check status of the download ----
occ_download_wait(calo_download)

##### Get download key ----
cdownload_key <- calo_download[1]

##### Note GBIF citation ----
citation <- gbif_citation(cdownload_key)

# Print citation to quote in writing/ODMAP
print(citation)

##### Load into environment ----
calod <- occ_download_get(cdownload_key) 
calodata <- occ_download_import(calod)

##### Write/save to file (raw / unprocessed copy) ----
raw_download <- file.path(out_dir, "raw_GBIF")
dir.create(raw_download, recursive = TRUE, showWarnings = FALSE) # create if does not exist
write.csv(calodata, file = file.path(raw_download, "Calo_raw_GBIF_download_29012025.csv"))

##### *START HERE if re-running without re-downloading ----

# Prereqs - run line to load "raw_download"
raw_download <- file.path(out_dir, "raw_GBIF")

# Load back in if going back to this script at a later date:
calodata <- fread(file.path(raw_download, "Calo_raw_GBIF_download_29012025.csv"))

# Take a look at the downloaded data
calodata # 65 records found (29.01.2025 gbif updated citation download)
 
##### GBIF cleaning ----
# Create and fill a list with only the 'data' section for each species
calodata <- calodata[ , c("decimalLongitude", "decimalLatitude")]

# Clean by removing zeros/NAs from the lat/lon columns
calodata[calodata==0] <- NA

# Get rid of new NAs that were formerly 0's within the lat/lon columns
calodata <- na.omit(calodata)

# Rename lat/lon columns to make them shorter/simpler
names(calodata)[names(calodata)=="decimalLongitude"] <- "lon"
names(calodata)[names(calodata)=="decimalLatitude"] <- "lat"


#==============================================================================#
#                     3. Survey data - Focal Plant (Calo) ----
#==============================================================================#

##### Calophyllum survey data ----
# Add in other calo occurrence sources from collaborator field surveys 
CVB <- read.csv(file.path(surveys, "CVB_Calo_2018_Cpanic_only.csv"))
TEAM <- read.csv(file.path(surveys,"TEAMS_2018_Cpanic_only.csv"))
VERO <- read.csv(file.path(surveys,"VERO_MHS_2018_Cpanic_only.csv"))
AMANDA <- data.frame(lon = 49.20856, lat = -17.91439, origin = "AMANDA") # plot 5 from Betampona survey site had Calophyllum present, coordinate information taken from single location (Amanda Armstrong)
AMAN <- data.frame(lon = 49.20856, lat = -17.91439)

##### Survey cleaning ----
# Align external sources to match gbif calo records to join together
CVB_ <- CVB[, c('LONGITUD', 'LATITUDE')]
TEAM_ <- TEAM[, c('SiteLat','SiteLong')]
VERO_ <- VERO[, c('longitude', 'latitude')]
gbifcalo_ <- calodata[, c('lon','lat')]

# Update colnames
colnames(CVB_) <- c("lon","lat")
colnames(TEAM_) <- c("lon","lat")
colnames(VERO_) <- c("lon","lat")

# Add origin column to each dataset
gbifcalo_$origin <- "gbifcalo"
CVB_$origin <- "CVB"
TEAM_$origin <- "TEAM"
VERO_$origin <- "VERO"

##### Combine sources together ----
calo_join <- rbind(gbifcalo_, CVB_, TEAM_, VERO_, AMANDA)

##### Remove any NA lat/lon columns ----
calo_join <- na.omit(calo_join)

##### Valid combined Calo records ----
# Print the number of valid occurrences: 1066
cat("Number of valid, combined Calophyllum paniculatum occurrences in Madagascar:", nrow(calo_join))

##### Write Calo to file ----

# Make a spat vector out of the remaining occ records imported from gbif
calov <- vect(calo_join, geom=c("lon", "lat"), crs="EPSG:4326", keepgeom=FALSE)
calovp <- project(calov, crs, partial=FALSE) # crs = ESRI:102022

coordscalo84 <- terra::crds(calov) # crs = EPSG: 4326 (native WGS 1984)
coordscalo <- terra::crds(calovp) # crs = ESRI:102022

##### Create clean_unprocessed presence sub_out_dir for storing GBIF and survey presences ----
# sub output folder for cleaned but pre-processed data (i.e. duplicates not yet removed etc.)
pres_out_dir <- file.path(out_dir, "presence")
dir.create(pres_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

# Write to CSV
write.csv(coordscalo84, file=(file.path(out_dir, "coordscalo_WGS1984.csv")), row.names=FALSE)
write.csv(coordscalo, file=(file.path(out_dir, "coordscalo_ESRI102022.csv")), row.names=FALSE)

# Tidy up
rm(CVB, CVB_, TEAM, TEAM_, VERO, VERO_, surveys, gbifcalo_, AMANDA, AMAN, calo_join, calov, calodata, coordscalo, coordscalo84)


#==============================================================================#
#                          4. GBIF - Wilt (VERT & LEPT) ----
#==============================================================================#

##### Verticillium / Leptographium genera in Africa ---- 

##### Create json request ---- 
# Create pred key for json format/talking directly with gbif (to be able to download a specific citation for writing up)
# Taxon key taken from gbif.org before running rgbif download (Verticillium Nees, 1816 (genus) taxonKey: 9535878)
# (Leptographium Lagerb. & Melin (genus) taxonKey: 2569086)
vertlept_download <- occ_download(
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("continent", "Africa"),
  pred_or(pred("taxonKey", 9535878), pred("taxonKey", 2569086)),
  format = "SIMPLE_CSV",
  user = user,
  pwd = pwdd,
  email = email
)

##### Check status of the download ----
# Check the status of the download
occ_download_wait(vertlept_download)

##### Get download key ----
# Get the download key
vldownload_key <- vertlept_download[1]

##### Note GBIF citation ----
# Get the citation
vl_citation <- gbif_citation(vldownload_key)

# Print citation to quote in writing/ODMAP
print(vl_citation)
# [1] "GBIF Occurrence Download https://doi.org/10.15468/dl.284amw Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2025-01-29"

##### Load into environment ----
# Download to working directory
vertleptd <- occ_download_get(vldownload_key) 
vertleptdata <- occ_download_import(vertleptd)

##### Unique taxonKey ----
vertlept_taxonKeys <- unique(vertleptdata$taxonKey) # 18 unique taxonKeys

##### Write/save to file (raw / unprocessed) ----
raw_download <- file.path(out_dir, "raw_GBIF")
dir.create(raw_download, recursive = TRUE, showWarnings = FALSE) # create if does not exist
write.csv(vertleptdata, file = file.path(raw_download, "vertlept_raw_GBIF_download_29012025.csv"))


##### *START HERE if re-running without re-downloading ----
# Load back in if going back to this script at a later date:
#vertleptdata <- fread(file.path(raw_download, "vertlept_raw_GBIF_download_29012025.csv"))

##### GBIF cleaning ----
# Create and fill a dataframe with only the 'data' section
vertlept_occurrences <- vertleptdata[ , c("decimalLongitude", "decimalLatitude")]

# Update lat-lon colname
colnames(vertlept_occurrences)<- c("lon","lat")

# Clean by removing zeros/NAs from the lat/lon columns
vert_occurrences <- vertlept_occurrences %>%
  filter(lon != 0, lat != 0) %>%
  na.omit()

# Print the first few rows of the cleaned data
print(head(vertlept_occurrences))

# Check/clean vertlept data further - potential swapped columns 
# This has been commented out as the GBIF download of 29.01.2025 contains cleaning coordinates in query
# Define the valid ranges for Africa
lat_range <- c(-34.8, 37.3)
lon_range <- c(-17.5, 51.4)

# Apply the function to each row
vertlept_occurrences[c("lat", "lon")] <- t(apply(vertlept_occurrences[c("lat", "lon")], 1, function(row) swap_coords(row[1], row[2])))

# Remove points outside the valid ranges
vertlept_clean <- vertlept_occurrences[vertlept_occurrences$lat >= lat_range[1] & vertlept_occurrences$lat <= lat_range[2] & 
                              vertlept_occurrences$lon >= lon_range[1] & vertlept_occurrences$lon <= lon_range[2], ]

# Calculate number of removed points
removed_points <- nrow(vertlept_occurrences) - nrow(vertlept_clean)

# Print number of removed points
cat("Number of points removed:", removed_points, "\n") # 2 rows removed after cleaning lon lat cols

# View removed points
invalid_points <- vertlept_occurrences[!(vertlept_occurrences$lat >= lat_range[1] & vertlept_occurrences$lat <= lat_range[2] & 
                                     vertlept_occurrences$lon >= lon_range[1] & vertlept_occurrences$lon <= lon_range[2]), ]
print(invalid_points)

##### Valid vertlept records ----
cat("Number of valid Verticillium & Leptographim (genera) occurrences in Africa:", nrow(vertlept_clean)) # 85

# Make a spat vector out of the remaining occ records imported from gbif
vertleptv <- terra::vect(vertlept_clean, geom=c("lon", "lat"), crs="EPSG:4326", keepgeom=FALSE)
vertleptp <- terra::project(vertleptv, crs, partial=FALSE)

extMDG <- ext(templateA) # Creates an extent from templateA raster of just MDG area to use as the "erase" geometry below

# Need to make an extent equivalent in projection ESPG:4326
extMDG_84 <- project(extMDG, "ESRI:102022" ,"EPSG:4326")

# Crop the template raster so that it does not include Madagascar for the Verticillium/Leptographium records
vertleptnom <- terra::erase(vertleptp, extMDG)
vertleptnom84 <- terra::erase(vertleptv, extMDG_84)

# removes occurrences from island of Mada - there are still occurrences outside of AFR mainland
vertleptcrop <- terra::crop(vertleptnom, templateB, ext=FALSE) # ESRI:102022
vertleptcrop84 <- terra::crop(vertleptnom84, templateB_84, ext=FALSE) # crs = EPSG:4326

image(templateAB)
points(vertleptcrop)

##### Write vertlept to file ----
# Write intermediate file - no need to re-download GBIF again
# Extract geometry 
coordsvertlept <- terra::crds(vertleptcrop) # crs = ESRI:102022
coordsvertlept84 <- terra::crds(vertleptcrop84) # crs = EPSG:4326

# Write to CSV
write.csv(coordsvertlept, file.path(pres_out_dir, "coordsvertleptlept_ESRI102022.csv"), row.names = FALSE)
write.csv(coordsvertlept84, file.path(pres_out_dir, "coordsvertlept_WGS1984.csv"), row.names = FALSE)

# Tidy up leaving only vertleptcrop for next sections of this script
remove(vertleptcrop84, vertleptlept_taxonKeys, coordsvertlept, coordsvertlept84, vertlept_clean, vertlept_occurrences, vertleptnom, vertleptp, vertleptv, extMDG, invalid_points, removed_points, lat_range, lon_range)


#==============================================================================#
#                       5. Remove duplicate presences (1km) ----
#==============================================================================#

# # Uncomment to load back in if re-running just this section
# # Re-load "coordscalo_ESRI102022.csv" (output_dir)
# calovp <- fread(file.path(out_dir, "coordscalo_ESRI102022.csv")) # 1066
# # Convert to vector
# calovp <- vect(calovp, geom=c("x", "y"), crs="ESRI:102022", keepgeom=F)
# 
# # Re-load "coordsvertlept_ESRI102022.csv" (output_dir)
# vertleptcrop <- fread(file.path(out_dir, "coordsvertlept_ESRI102022.csv")) # 34
# # Convert to vector
# vertleptcrop <- vect(vertleptcrop, geom=c("x", "y"), crs="ESRI:102022", keepgeom=F)


##### Calophyllum ----
calocrop <- calovp # calovp comes from earlier when combined from gbif and all other survey sources

# Make sure occ vect has attributes first
coordsc <- terra::crds(calocrop)

# Extract cell numbers from calocrop
cellnumbersc <- terra::extract(templateA, calocrop, cells=T)

# Get only valid matches between calocrop locations and templateA
matched_data <- merge(
  data.frame(
    x = coordsc[,1], 
    y = coordsc[,2], 
    ID = 1:nrow(coordsc)
  ),
  cellnumbersc,
  by = "ID",
  all.x = TRUE
)

# create a df where duplication within the same cell can be removed
calo_df <- data.frame(
  x = matched_data$x,
  y = matched_data$y,
  cell = matched_data$cell,
  ID = matched_data$ID
)

# Remove any occurrences that are duplicated within the same 1km cell
calo_unique <- calo_df[!duplicated(calo_df$cell), ]

# Convert back to spatvector
calo_final <- terra::vect(calo_unique, geom=c("x", "y"), crs=terra::crs(calocrop))# returns 52 
calonodupsp <- calo_final

# Calo final records: 62

# Repeat for Verticillium & Leptographium genera presences in Africa

##### Verticillium / Leptographium ----
coordsv <- terra::crds(vertleptcrop)

# Extract cell numbers
cellnumbersv <- terra::extract(templateB, vertleptcrop, cells=TRUE) # reduces from 76 to 27

# Get only valid matches between calocrop locations and templateA
matched_data_v <- merge(
  data.frame(
    x = coordsv[,1], 
    y = coordsv[,2], 
    ID = 1:nrow(coordsv)
  ),
  cellnumbersv,
  by = "ID",
  all.x = TRUE
)

# create a df where duplication within the same cell can be removed
vertlept_df <- data.frame(
  x = matched_data_v$x,
  y = matched_data_v$y,
  cell = matched_data_v$cell,
  ID = matched_data_v$ID
)

vertlept_unique <- vertlept_df[!duplicated(vertlept_df$cell),] # Remove dups
vertleptnodups <- terra::vect(vertlept_unique, geom=c("x","y"), crs=terra::crs(vertleptcrop)) # Add col values back in
vertleptnodupsp <- project(vertleptnodups, crs) # project to AEA

# vertlept final records: 34

# Intermediate saving
# Save both sets of presence points at 1 km resolution in RData format

# Convert them back to dfs after being spatvectors
calonodupsp <- terra::as.data.frame(calonodupsp, geom="XY") # 62 final presences
vertlept_nodupsp <- terra::as.data.frame(vertleptnodupsp, geom="XY") # 34 final combined wilt presences

# Create sub-folder of presence file for duplicate free presences
# sub output folder for cleaned but pre-processed data (i.e. duplicates not yet removed etc.)
dupfree_out_dir <- file.path(pres_out_dir, "dupfree")
dir.create(dupfree_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

##### Write to file ----
save(calonodupsp, file=(file.path(dupfree_out_dir, "Calo_00_noDups.RData"))) # use if using SDM package without creating BG points first
save(vertlept_nodupsp, file=(file.path(dupfree_out_dir, "vertlept_00_noDups.RData"))) # use if using SDM package without creating BG points first


# Tidy up
rm(calo_final, cellnumbersc, cellnumbersv, calovp, calocrop ,calo_unique, calo_df, coordsv, coordsc, matched_data, matched_data_v)
rm(vertlept_df, vertlept_nodups, vertlept_unique, vertleptnodups)
rm(calo, vertlept, vertleptcrop, vertleptnodupsp)


#==============================================================================#
#                      6. Prepare data for Background pts ----
#==============================================================================#

# # Template CRS (African Equal Albers / WGS 1984 Albers for Africa)
# crs <- "ESRI:102022"
# 
# ## ***One-time only run*** (comment out after)
# # Intermediate save for reproducibility
# # Create plots output sub-folder
# int_out_dir <- file.path(datafolder__00, "intermediate")
# dir.create(int_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist
# 
# ### Vectorise each template rasters ----
# vtemplateA <- as.polygons(templateA) # v = vector
# vtemplateB <- as.polygons(templateB)
# vtemplateAB <- as.polygons(templateAB)
# 
# ## Write vector versions of study areas
# writeVector(vtemplateA, filename=file.path(int_out_dir, paste0("vtemplateA.shp")))
# writeVector(vtemplateAB, filename=file.path(int_out_dir, paste0("vtemplateAB.shp")))
# writeVector(vtemplateB, filename=file.path(int_out_dir, paste0("vtemplateB.shp")))
# 
# ## Prepare eco region data ----
# # Locate eco regions shapefile, align crs, and crop to study extent
# # Eco regions
# ecoregions <- file.path(datafolder__00, "input", "ecoregions")
# ecos <- list.files(ecoregions, pattern = "Africa_Ecoregions.shp", full.names = T)
# AFR_ecos <- terra::vect(ecos)
# 
# ### Project ecoregions to proj crs ---
# AFR_ecosp <- terra::project(AFR_ecos, templateAB)
# 
# ### Crop the ecoregions to study extents ----
# # A = MDG ; B = AFR ; AB = both
# ecoregsA <- crop(AFR_ecosp, vtemplateA)
# ecoregsB <- crop(AFR_ecosp, vtemplateB)
# ecoregsAB <- crop(AFR_ecosp, vtemplateAB)
# 
# ### Write cropped ecoregion files ----
# writeVector(ecoregsA, filename=file.path(int_out_dir, paste0("ecoregsA.shp")))
# writeVector(ecoregsB, filename=file.path(int_out_dir, paste0("ecoregsB.shp")))
# writeVector(ecoregsAB, filename=file.path(int_out_dir, paste0("ecoregsAB.shp")))
# 
# ## Prepare elevation data ----
# # elevation range 10-2087m for 8 x Verticillium records in GBIF, no other elevational info directly in GBIF DB)
# # using worldclim elevation data extract elevation values for all presences
# #
# # Tidy up
# rm(ecos, AFR_ecos, elev, elv, elevation, AFR_ecosp, ecoregsAB, templateAB, elv_AB, vtemplateAB)

# Template CRS (African Equal Albers / WGS 1984 Albers for Africa)
crs <- "ESRI:102022"

### Set bg points output ----
# subfolder location for background points and their buffer areas
bg_out_dir <- file.path(out_dir, "background")
dir.create(bg_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

# N.B. **All objects should be spat vector objects when running bg pt generation functions**
# Load back in pre-processed intermediate files e.g. cropped study areas, eco regions, elevation

### Re-load in presence points ----
# Uncomment to load back in dup-free occurrences if re-running this section
# Re-load "Calo_00_noDups.RData" (output_dir)
load(file.path(out_dir, "presence", "dupfree", "Calo_00_noDups.RData")) # calonodupsp (total: 62 duplicate free occs)

# Re-load "vertlept_00_noDups.RData" (output_dir)
load(file.path(out_dir, "presence", "dupfree", "vertlept_00_noDups.RData")) # vertleptnodupsp (total: 34 duplicate free occs)

### Vectorise presence points (spat vectors) ----
# Create vector versions of the presence points for use in this section
calop <- terra::vect(calonodupsp, geom=c("x", "y"), crs=crs) # p = presence
vertleptp <- terra::vect(vertlept_nodupsp, geom=c("x", "y"), crs=crs) # p = presence

### Load in intermediate data (ecoregion) ----
# Locate intermediate sub-folder
int_out_dir <- file.path(datafolder__00, "intermediate")
dir.create(int_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

# List files in intermediate/processed data (ecoregions, study areas)
ecoregs <- list.files(int_out_dir, pattern = "ecoregs.*\\.shp$", full.names = T)
templates <- list.files(int_out_dir, pattern = "template.*\\.shp$", full.names = T)

# Load MDG intermediate data
ecoregsA <- terra::vect(ecoregs[1]) # ecoregion
vtemplateA <- terra::vect(templates[1]) # study extent (vector)

# Load AFR intermediate data
ecoregsB <- terra::vect(ecoregs[3]) # ecoregion
vtemplateB <- terra::vect(templates[3]) # study extent (vector)

### Prepare Cryphalus data for Vert-Lept background points ----

# load in Cryphalus occurrence points from GBIF to add to wilt occs
Cryphalus <- list.files(file.path(datafolder__00, "input", "GBIF"), full.names = T)
Cry <- fread(Cryphalus)
# Make Cryphalus and Vert-Lept match to merge
names(Cry)[names(Cry) == c("decimalLatitude", "decimalLongitude")] <- c("lat", "lon")

# Swap column order to match merge
Cry <- Cry[, c(2,1)] # switch so lon first before lat in column index order

Cryv <- terra::vect(Cry, geom=c("lon", "lat"), crs="EPSG:4326")

# Remove points outside study area extent (AFR / vtemplateB)
Cryp <- terra::project(Cryv, crs)
Crypcrop <- terra::crop(Cryp, vtemplateB)

# As df for merging
Crypcropdf <- as.data.frame(Crypcrop)

# Re-load "vertlept_00_noDups.RData" (output_dir)
load(file.path(out_dir, "vertlept_00_noDups.RData")) # vertleptnodupsp (total: 34 duplicate free occs)

# Basic data wrangling to make columns and names match with Cryv
# Remove extra columns in vert-lept occurrences
vertlept <- subset(vertlept_nodupsp, select= -c(cell,ID))

# Get lat/lon coordinates of vert-lept occurrences
vertleptv <- terra::vect(vertlept, geom=c("x", "y"), crs=crs)
vertlept84 <- terra::project(vertleptv, "EPSG:4326") # convert to ESPG:4326
coordsvl_84 <- terra::crds(vertlept84) # Get lat/lon
vl84 <- as.data.frame(coordsvl_84) # new df with updated lat/lon
names(vl84)[names(vl84) == c("x", "y")] <- c("lon", "lat")

# rbind both "Cry" and "vl84" while dfs before vectorising as one merged spat vector
wilt_beetle <- rbind(vl84, Crypcropdf)

# Vectorise combined wilt and bark beetle presences
wilt_bark_beetle <- terra::vect(wilt_beetle, geom= c("lon", "lat"), crs="EPSG:4326")

# Tidy 
rm(calonodupsp, vertlept_nodupsp, ecoregs, aois, templates, coordsvl_84, Cry, Cryp, Crypcrop, Crypcropdf, Cryv,
   vertlept84, vertlept_taxonKeys, vertleptv, vertlept, vl84, wilt_beetle, Cryphalus)


#==============================================================================#
#                       6a. Calo - Create background points ----
#==============================================================================#

#   1. CALO: Background points are randomly generated within the eco region of the presences 
#       justification: Calophyllum is a specialist forest dwelling species only found within humid and sub-humid forests in Madagascar

species <- c("calo", "vertlept")

# Species / location specific:
# eco region
# Set seed for reproducibility before each set of bg point generation
## Calo BG points (620) ----
set.seed(101)

# Ecoregion buffer (ECO_ID spatial boundary of presence points)
calo_bg_ecor <- generate_ecoregion_bg_pts(calop, vtemplateA, ecoregsA, species_name = "calo", n_points = 620)
gc()

set.seed(101)
## Whole area BG points (620) ----
vertlept_rand_bg <- random_bg_whole_area(calop, vtemplateA, species = "calo", n_random=620, buffer_dist = 1500) # 1500m buffer from presences (accounts for spatial duplicates within same cell)



#==============================================================================#
#                   6b.Vert-Lept Create background points VERT/LEPT ----
#==============================================================================#


# - 2. VERT-LEPT: Geographical buffer based on dispersal limitation / previous range limits of Cryphalus bark beetle vector
#       justification: Little is known about dispersal limitation of the fungal wilt itself, literature suggests it moves max 1km
#       so accounting for the fact its likely dispersal vector is via bark beetle Cryphalus......
#       dispersal distance = Cryphalus shortest distance between two AFR occurrence points?????

## Vert-Lept 2far BG points (340) ----
set.seed(101)
# Use "process_2degfar_bg_pts" to build correct geographical buffers for each presence cluster
vertlept_2far_bg <- process_2degfar_bg_pts(wilt_bark_beetle, vtemplateB, n_clusters = 3, n_random = 340)

# Plot vertlept_2far_bg (??? re-do as a figure later)
plot(vertlept_2far_bg$remaining_area)
plot(vertlept_2far_bg$random_points, add = TRUE, col = "red")
plot(vertlept_2far_bg$buffers, add = TRUE, border = "blue", col = NA)
plot(vertlept_2far_bg$original_points, add = TRUE, col = "green")

set.seed(101)
## Vert-Lept whole area BG points (340) ----
vertlept_rand_bg <- random_bg_whole_area(vertleptp, vtemplateB, species="vertlept", n_random=340, buffer_dist = 1500) # 1500m buffer from presences (accounts for spatial duplicates within same cell)



#==============================================================================#
#                     6c. Optional plot bg points ----
#==============================================================================#

## Basic plots to check

# Load back in bg points and buffer areas for each method:
# List files: bg buffer areas
bgbuff_list <- list.files(file.path(bg_out_dir), pattern = ".shp", full.names = T)

for (i in seq_along(bgbuff_list)) {
  
  # Get the basename of the buffer 
  file_name <- tools::file_path_sans_ext(basename(bgbuff_list[i]))
  
  # Load the buffer
  current_buff <- vect(bgbuff_list[i])
  
  # Assign buffer its basename
  assign(file_name, current_buff, envir = .GlobalEnv)
  
}

# Load back in background points
bgptslist <- list.files(file.path(bg_out_dir), pattern = ".RData", full.names = TRUE)
for (i in seq_along(bgptslist)) {
  
  # Get the basename of the buffer 
  file_name <- tools::file_path_sans_ext(basename(bgptslist[i]))
  
  # Create a new environment to load the data
  temp_env <- new.env()
  
  # Load the buffer into temporary environment
  loaded_name <- load(bgptslist[i], envir = temp_env)
  
  # Get the object from temp environment and assign it with the file name
  assign(file_name, get(loaded_name, envir = temp_env), envir = .GlobalEnv)
}


##  Plots with bg pts ----
# Create calophyllum comparison plot (png to file)
# Direct comparison of ecoregions vs whole area random
# Keep constraint regions and bg points in order: (method 1, method 2, bg1, bg2..)
plot_bg_comparison_calo(calop, vtemplateA, 
                        calo_n620_ecoregion, 
                        calo_rand_bg_whole_n620, 
                        calo_n620_ecor_bg_pts, 
                        calo_rand_bg_whole_n620_bg_pts, "Calophyllum")
gc()

# Create Vert-Lept comparison plot (png to file)
# Direct comparison of 2degrees far method vs whole area random
plot_bg_comparison_vert(vertleptp, vtemplateB, 
                        vertlept_2far_n340, 
                        vertlept_rand_bg_whole_n340, 
                        vertlept_2far_n340_bg_pts, 
                        vertlept_rand_bg_whole_n340_bg_pts, "Verticillium + Leptographium")
gc()



rm(current_buff, ecoregsA, ecoregsB, temp_env)


#==============================================================================#
#                     7. Remove duplicate background points ----
#==============================================================================#


# Re-load study extent rasters if not yet in environment
# AOIs of template extents
aois <- file.path(datafolder__00, "input", "AOI")

# Just African continent no islands, no MDG
templateB <- rast(file.path(aois, "AFR", "AFR_aea_template1km.tif"))

# Just MDG
templateA <- rast(file.path(aois,"MDG", "MDG_aea.tif"))


## CALO BG PTS DUPS ----

# Rename to continue from previous section
c_back_n620_ecor <- calo_n620_ecor_bg_pts # if continuing from previous step
c_back_n620_whole <- calo_rand_bg_whole_n620_bg_pts # if continuing from previous step

### Extracting calophyllum background points from environmental cells ----
cellno_n620_ecor <- terra::extract(templateA, c_back_n620_ecor, cells=T)
cellno_n620_whole <- terra::extract(templateA, c_back_n620_whole, cells=T) 

### Remove duplicated rows by their cellnumbers ----
c_back_n620_ecor_nodup <- c_back_n620_ecor[!duplicated(cellno_n620_ecor[,'cell']),] # reduces 1 of the bg points are within 1km of one another
c_back_n620_whole_nodup <- c_back_n620_whole[!duplicated(cellno_n620_whole[,'cell']),] # reduces 1 of the bg points are within 1km of one another

### Convert duplicate free background points as dataframe ----
calo_bg_n620_ecor_nodup <- as.data.frame(c_back_n620_ecor_nodup)
calo_bg_n620_whole_nodup <- as.data.frame(c_back_n620_whole_nodup)

### Write background points at 1 km resolution in .RData format ----
save(calo_bg_n620_ecor_nodup, file=(file.path(bg_out_dir, "calo_bg_n620_ecor_nodup.RData")))
save(calo_bg_n620_whole_nodup, file=(file.path(bg_out_dir, "calo_bg_n620_whole_nodup.RData")))


## VERT & LEPT BG PTS DUPS ----

# repeat duplicate removal for Verticillium & Leptographium background points
# Rename to continue from previous section
v_back_n340_2far <- vertlept_2far_n340_bg_pts # if continuing from previous step
v_back_n340_whole <- vertlept_rand_bg_whole_n340_bg_pts # if continuing from previous step

### Extracting Vert-Lept background points from environmental cells ----
cellno_n340_2far <- terra::extract(templateB, v_back_n340_2far, cells=T)
cellno_n340_whole <- terra::extract(templateB, v_back_n340_whole, cells=T) 

### Remove duplicated rows by their cellnumbers ----
v_back_n340_2far_nodup <- v_back_n340_2far[!duplicated(cellno_n340_2far[,'cell']),] # reduces by 0 - none of the bg points are within 1km of one another
v_back_n340_whole_nodup <- v_back_n340_whole[!duplicated(cellno_n340_whole[,'cell']),] # reduces by 0 - none of the bg points are within 1km of one another

### Convert duplicate free background points as dataframe ----
vertlept_bg_n340_2far_nodup <- as.data.frame(v_back_n340_2far_nodup)
vertlept_bg_n340_whole_nodup <- as.data.frame(v_back_n340_whole_nodup)

### Write background points at 1 km resolution in .RData format ----
save(vertlept_bg_n340_2far_nodup, file=(file.path(bg_out_dir, "vertlept_bg_n340_2far_nodup.RData")))
save(vertlept_bg_n340_whole_nodup, file=(file.path(bg_out_dir, "vertlept_bg_n340_whole_nodup.RData")))


# Tidy up
rm(calo_n620_ecor_bg_pts, calo_n620_ecoregion, cellno_n340_2far, cellno_n340_whole, current_buff,
   calo_rand_bg_whole_n620, calo_rand_bg_whole_n620_bg_pts,
   ecoregsA, ecoregsB, temp_env, templateA, templateB, templateAB, v_back_n340_2far, v_back_n340_whole,
   v_back_n340_whole_nodup, v_back_n340_2far_nodup,
   vertlept_2far_n340, vertlept_2far_n340_bg_pts, vertlept_rand_bg_whole_n340, vertlept_rand_bg_whole_n340_bg_pts,
   vtemplateA, vtemplateB, wilt_bark_beetle)


#==============================================================================#
#                   8. Combining presence-background as "occs" ----
#==============================================================================#

##### Create occ_out_dir for storing combined occurrences ----
occ_out_dir <- file.path(out_dir, "combined_occs")
dir.create(occ_out_dir, recursive = TRUE, showWarnings = FALSE) # create if does not exist

## 8a. CALO OCCS ----
 # Prepare the presences data to contain a column indicating 1 for presence 
cpresences <- as.data.frame(calop, geom="XY") # 62 presences
cpresences['occ'] = 1

# reduce the columns in cpresences so that they can bind with background data
cpresences <- cpresences[-(1:2)] # this removes ID, cell etc. leaving just x,y,occ

# Prepare the background data to contain a  column indicating 0 for 'background'
calo_bg_n620_ecor_nodup$occ <- 0
calo_bg_n620_whole_nodup$occ <- 0

# Bind these two data sets together so presence and background points are together
calo_occ_ecor <- rbind(cpresences, calo_bg_n620_ecor_nodup) # 1 = bg method: ecoregions
calo_occ_whole <- rbind(cpresences, calo_bg_n620_whole_nodup) # 2 = bg method: whole area

# Save ready for the SDMs, should be 682 duplicate free occurrences (62 presences (1), 619 background (0))
save(calo_occ_ecor, file=(file.path(occ_out_dir, "Calo_occs_combined_ecor_619.RData"))) # bg method: ecoregions
save(calo_occ_whole, file=(file.path(occ_out_dir, "Calo_occs_combined_whole_619.RData"))) # bg method: whole area


## 8b.VERT & LEPT OCCS ----
# Prepare the combined presences data to contain a column indicating 1 for presence 
vpresences <- as.data.frame(vertleptp, geom="XY") # 34 presences

# create occ column
vpresences['occ'] = 1 # occ = 1 for presences
vlpresences <- vpresences[-(1:2)] # removes cell & ID columns not required for joining presences and absences

# Prepare the background data to contain a  column indicating 0 for 'background'
vertlept_bg_n340_2far_nodup$occ <- 0
vertlept_bg_n340_whole_nodup$occ <- 0

# Bind these two data sets together so presence and background points are together
vertlept_occ_2far <- rbind(vlpresences, vertlept_bg_n340_2far_nodup)
vertlept_occ_whole <- rbind(vlpresences, vertlept_bg_n340_whole_nodup)

# Save the df containing the presences (34 = 1) and background points(340 = 0))
save(vertlept_occ_2far, file=(file.path(occ_out_dir, "vertlept_occs_2far_374.RData")))
save(vertlept_occ_whole, file=(file.path(occ_out_dir, "vertlept_occs_whole_374.RData")))

# Tidy up
rm(list = ls())
gc()


#==============================================================================#
#                           ---- End of work flow ----
#==============================================================================#