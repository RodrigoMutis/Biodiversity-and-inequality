# Biodiversity-and-inequality 

# Install necessary packages 
install.packages("sf")         
install.packages("ineq")
install.packages("vegan")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tibble")
install.packages("spatialreg")
install.packages("spdep")

# Load libraries
library(sf)          # For geospatial data
library(ineq)        # For Gini Index
library(vegan)       # For biodiversity indices
library(dplyr)       # For data manipulation
library(tidyr)       # For data manipulation
library(tibble)      # For data manipulation
library(spatialreg)  # For Spatial regression
library(spdep)       # For spatial dependence

# Load bird data
birddata <- read.delim("path/to/birddata.txt")

# Clean bird data
birddata <- birddata %>%
  filter(CATEGORY %in% c("species", "issf"))

birddata <- birddata [, -c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]
birddata$OBSERVATION.COUNT <- gsub("x", "1", birddata$OBSERVATION.COUNT)
birddata$OBSERVATION.COUNT <- as.numeric(trimws(birddata$OBSERVATION.COUNT))
birddata <- birddata %>% replace_na(list(OBSERVATION.COUNT = 1))

# Convert in shape
birddata <- st_as_sf(birddata, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# Load value city data
citydata <- st_read("path/to/citydata.shp") # Polygons

# Check CRS of both shapefiles
st_crs(citydata)
st_crs(birddata)

# If CRS is different, reproject one of them
birddata <- st_transform(birddata, st_crs(citydata))

# Reproject citydata to UTM (change the EPSG code according to city zone)
citydata_utm <- st_transform(citydata, crs = 32618)

# Reproject birddata to UTM (change the EPSG code according to city zone)
birddata_utm <- st_transform(birddata, crs = 32618)

# Generates a hexagonal grid
hex_grid <- st_make_grid(citydata_utm, cellsize = 500, square = FALSE) %>% 
  st_as_sf() %>% 
  mutate(ID = row_number())

# Calculate centroids
centroids <- st_centroid(hex_grid)

centroid_coords <- st_coordinates(centroids) %>%
  as.data.frame() %>%
  rename(Longitude = X, Latitude = Y)
  
centroid_data <- hex_grid %>%
  st_drop_geometry() %>%  
  bind_cols(centroid_coords)  
  
# Intersect data  
citydata_hex <- st_intersection(citydata_utm, hex_grid) %>% 
  st_join(hex_grid, by = "ID")
  
birddata_hex <- st_intersection(birddata_utm, hex_grid) %>% 
  st_join(hex_grid, by = "ID")
  
merged_data <- st_join(citydata_hex, birddata_hex, by = "ID")

# Calculate biodiversity index 

lassify species by abundance with new naming
species_abundance <- birddata_hex %>%
  st_drop_geometry() %>%
  group_by(SCIENTIFIC.NAME) %>%
  summarise(
    total_obs = sum(OBSERVATION.COUNT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_obs)) %>%
  mutate(
    rank = row_number(),
    pct_rank = percent_rank(total_obs),
    group = case_when(
      pct_rank <= 0.10 ~ "common",        # Top 10% most abundant
      pct_rank >= 0.90 ~ "potential_error", # Bottom 10% (to exclude)
      TRUE ~ "rare"                       # Middle 80% (formerly "uncommon")
    )
  )

# Get species lists
common_species <- species_abundance %>% 
  filter(group == "common") %>% 
  pull(SCIENTIFIC.NAME)

rare_species <- species_abundance %>% 
  filter(group == "rare") %>% 
  pull(SCIENTIFIC.NAME)

# Exclude potential errors (bottom 10%)
birddata_filtered <- birddata_hex %>%
  filter(SCIENTIFIC.NAME %in% c(common_species, rare_species))

# Calculate biodiversity indices with updated names
calculate_indices <- function(data, suffix) {
  if(nrow(data) == 0) return(data.frame(ID = character()))
  
  abundance_matrix <- data %>%
    st_drop_geometry() %>%
    group_by(ID, SCIENTIFIC.NAME) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = SCIENTIFIC.NAME,
      values_from = n,
      values_fill = 0
    ) %>%
    column_to_rownames("ID")
  
  data.frame(
    ID = rownames(abundance_matrix),
    richness = vegan::specnumber(abundance_matrix),
    shannon = vegan::diversity(abundance_matrix),
    simpson = vegan::diversity(abundance_matrix, index = "simpson"),
    stringsAsFactors = FALSE
  ) %>%
    rename_with(~ paste0(., "_", suffix), -ID)
}

# Calculate indices
indices_common <- calculate_indices(
  birddata_filtered %>% filter(SCIENTIFIC.NAME %in% common_species),
  "common"
)

indices_rare <- calculate_indices(
  birddata_filtered %>% filter(SCIENTIFIC.NAME %in% rare_species),
  "rare"
)

# Merge data
hex_grid_with_indices <- hex_grid %>%
  left_join(indices_common, by = "ID") %>%
  left_join(indices_rare, by = "ID") %>%
  mutate(across(starts_with("richness"), ~ replace_na(., 0)))

Enhanced Spatial Regression Function
run_full_analysis <- function(data, prefix) {
  # Create response variables
  responses <- c("richness", "shannon", "simpson") %>%
    set_names() %>%
    map(~ paste0(., "_", prefix))
  
  # Prepare spatial weights
  coords <- st_centroid(data) %>% st_coordinates()
  nb <- knn2nb(knearneigh(coords, k = 5))
  lw <- nb2listw(nb, style = "W")
  
  # Run models for each response variable
  map_dfr(responses, function(response) {
    if(!response %in% names(data)) return(tibble())
    
    formula <- as.formula(paste(response, "~ gini + mean_valor"))
    
    # Run models with error handling
    models <- list(
      sar = tryCatch(
        lagsarlm(formula, data, lw, method = "eigen"),
        error = function(e) NULL
      ),
      sem = tryCatch(
        errorsarlm(formula, data, lw, method = "eigen"), 
        error = function(e) NULL
      ),
      sdm = tryCatch(
        lagsarlm(formula, data, lw, type = "mixed", method = "eigen"),
        error = function(e) NULL
      )
    )
    
    # Extract results
    tibble(
      model = names(models),
      aic = map_dbl(models, ~ if(!is.null(.)) AIC(.) else NA_real_),
      coefficients = map(models, ~ if(!is.null(.)) tidy(.) else NULL)
    )
  }, .id = "metric") %>%
    mutate(group = prefix)
}

# Run analyses for both groups
results <- bind_rows(
  run_full_analysis(
    final_data %>% filter(!is.na(richness_common)), 
    "common"
  ),
  run_full_analysis(
    final_data %>% filter(!is.na(richness_rare))), 
    "rare"
  )
)

# Format and display results
# Coefficients table
coefficients_table <- results %>%
  filter(!map_lgl(coefficients, is.null)) %>%
  unnest(coefficients) %>%
  select(group, metric, model, term, estimate, std.error, p.value) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    p.value = format.pval(p.value, digits = 2)
  )

# AIC comparison
aic_comparison <- results %>%
  select(group, metric, model, aic) %>%
  pivot_wider(names_from = model, values_from = aic) %>%
  arrange(group, metric)

