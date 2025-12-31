# anon-fungal-host-sdms



---

title: "README"

manuscript: "Title: Climate change facilitates fungal pathogen expansion while driving endemic host range contractions in a tropical biodiversity hotspot"

corresponding\_author: "xxx"

coauthors: "xxx, xxx, xxx, xxx, xxx"

update: "2025-12-31"

output: Rmd readme

---



\# Abstract

Non-native plant pathogens are reshaping ecosystems globally, yet their virulence and spread potential under future climate and land-use change remains underexplored, particularly in biodiversity hotspots where endemic hosts may be highly vulnerable. In Madagascar, a vascular wilt pathogen (Leptographium calophylli, formerly Verticillium) has been increasingly observed infecting native forest trees. We modelled the future spread of this pathogen under changing climate conditions and directly compared projections to those of a susceptible endemic host tree, Calophyllum paniculatum, to assess overlap, divergence, and implications for extinction risk. Using an ensemble species distribution modelling approach (Random Forest, Boosted Regression Trees, MaxEnt), we forecasted potential distributional ranges for both pathogen and host using the Coupled Model Intercomparison Project (CMIP6) Global Climate Model (GFDL-ESM4) under three shared socioeconomic pathways (SSP1–2.6, SSP3–7.0, SSP5–8.5) and three time periods (2011–2040, 2041–2070, 2071–2100). Ensembles showed high predictive performance (AUC > 0.97), with precipitation seasonality and moisture availability in the driest month as key drivers of wilt distribution. The pathogen was predicted to retain 68.5% of the current projected distribution by 2100, with westward expansion into humid and sub-humid ecoregions, consistent with climate-facilitated invasion dynamics. In contrast, C. paniculatum was forecast to contract severely (65.9% range loss by 2100) and shift south-eastward, leading to a partial spatial decoupling. However, overlap with pathogen distributions persisted, indicating sustained mortality risk. Our results demonstrate that climate change may facilitate fungal pathogen expansion while simultaneously reducing the range of endemic hosts, an asymmetric dynamic that could intensify biodiversity loss in island ecosystems. Even without human mediated, long-distance dispersal, biological invasions can expose native species to prolonged biotic pressure. By explicitly modelling future host–pathogen interactions, this study highlights the importance of incorporating disease threats into biodiversity forecasts and conservation strategies, particularly in tropical systems facing rapid environmental change.



\# ODMAP

Following Zurell et al., (2020), we include two ODMAP csv files for openness and transparency that include model and parameter decisions and justifications:

\#### ODMAP for wilt Leptographium calophylli

\#### ODMAP for host Calophyllum paniculatum



\# R scripts outline



\#### - functions.R

This functions.R script houses all main functions for the paper and is called in via source().



\####  - 00\_\_species\_occurence\_cleaning\_bg.R

This script gathers and pre-processes all required data for the chapter, and prepares them for use within R using sdm R package. 

This includes GBIF downloads for both species of interest and the closest available wilt pathogen genera.

&nbsp;- Calophyllum paniculatum (spp)

&nbsp;- Leptographium \& Verticillium (genera) - combined into one wilt dataset representing Leptographium calophylli (ssp)

Processing of calo and wilt occurrence data, coordinate cleaning.

Removing duplicates, spatial thinning, 

Creating background points and pseudoabsences (and thinning them) (various species-specific methods: ecor = sampled within ecoregion of presence points, 2far method, and whole area)



\####  - 01\_\_extract\_species\_enviro\_data.R

Pre-processing and preparation of environmental predictor variables required for modelling.

Includes cropping, reprojecting, stacking etc. of bioclimatic, topographical and LULCC (deforestation) layers.



All 19 bioclimatic variables were downloaded from CHELSAv2 (Philipps et al., 2022) for three SSP scenarios (SSP126, SSP370 and SSP585).

Four date ranges utilised: 1981-2010, 2011-2040, 2041-2070, 2071-2100.



In addition to CHELSA(v2) variables, deforestation predictions from "forestatrisk" project from Vieilledent et al., (2018).

Five date ranges utilised: current(fcc\_123), probability of deforestation in 2020 (prob2020), and deforestations for 2040 (fcc\_2040), 2060 (fcc\_2060) and 2080 (fcc\_2080).

current (CHELSA 1981-2010 \& fcc\_123)

future 2040 (CHELSA 2011-2040 \& fcc\_2040)

future 2060 (CHELSA 2041-2070 \& fcc\_2060)

future 2080 (CHELSA 2071-2100 \& fcc\_2080)



Aspect (downloaded from worldclim world\_elev\_2.c.tif) as a proxy for microclimates



Folder structure:

input/

├── species/

│   ├── calo\_occ.RData

│   └── vert\_occ.RData

├── environmental/

│   ├── current/

│   │   ├── bioclim/

│   │   ├── forest/

│   │   └── aspect.tif

│   ├── future/

│   │   ├── ssp126/

│   │   │   ├── 2011-2040/

│   │   │   ├── 2041-2070/

│   │   │   └── 2071-2100/

│   │   ├── ssp370/

│   │   └── ssp585/

output/

├── stacked/

└── extracted/



\#### - c01\_\_extract\_species\_enviro\_data.R

Adapted script to run extraction on University cluster machine.



\####  - 02\_\_variable\_selection\_by\_species.R

Correlation matrix against all nineteen CHELSA bioclimatic variables - removes any > 0.7 correlation coefficient

Based on ratio of occurrence data, cannot use more than 6 variables for calo, or 4 variables for wilt.

In addition, selected predictors < 0.7 and < 0.5 were ranked based on VIF scoring and higher VIF scores removed until optimum predictors remained

Additional variables of forest cover and aspect (in addition to CHELSA variables) were added manually.



\####  - 03\_\_model\_building.R

Weighted mean algorithm ensembles were built

alg\_sets <- c("rf", "brt", "maxent").

Each algorithm was replicated 10 times within each ensemble using bootstrapping method.



\#### - 04\_\_integrated\_perf\_eval.R \& model\_perf\_eval.R

Model performance was assessed using 5-fold cross validation (80:20 train/test split) based on five different performance metrics:

The area under the receiver operating characteristic curve (AUC), true skill statistic (TSS), Point-biserial correlation, Residual deviance.

An additional combined/weighted ranking was calculated based on the AUC, TSS and corrected for model stability (using their standard deviations) to get top 3 models to make projections and thresholds from.

During model prediction the following performance measures were assessed:

c("AUC", "TSS", "Kappa", "COR", "Deviance")





\#### - 05\_\_model\_prediction.R

Multiple scripts used in this section of workflow:

a. 05\_full\_predsel\_env\_dat\_stacking.R -> used to stack the extract predictor variables together across whole prediction study area (MDG) - so that trained inputs match exactly

b. c05\_variable\_stack\_cropping\_MDG. R - > once environmental data was stacked correctly in a., they were cropped to MDG extent

c. 05\_\_ensemble\_projection\_feasibility.R - > semi manual script to look at initial model projection feasibility and produce thresholds for future binary predictions

d. 05\_\_model\_prediction.R - > automated for loop that for each species and forecasted year/pathway, produces ensemble prediction based on weighted mean for AUC optimised for TSS

e. 05b\_\_model\_predictions\_binary.R - > conversions of continuous model predictions from d., to binary presence absence for range shift metric calculations using threshold values produced in 05\_\_ensemble\_projection\_feasibility.R 



\#### - 06\_\_rangeshift\_metrics\_consolidated.R

Computes habitat net change, centroid shifting, direction arrows, output plots





\#### - 07\_\_visualisations.R

Visualisations and plotting script.



\#### - 08\_\_uncertainty.R

visualises uncertainty within binary projections and predictions for each spp, time, scenario, threshold



\# - End -



