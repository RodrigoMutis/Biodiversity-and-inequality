## Biodiversity-and-inequality ##

# Install necessary packages 
install.packages("sf")         
install.packages("ineq")
install.packages("vegan")
install.packages("dplyr")
install.packages("tidyr")

# Load libraries
library(sf)          # For geospatial data
library(ineq)        # For Gini Index
library(vegan)       # For biodiversity indices
library(dplyr)       # For data manipulation
library(tidyr)       # For data manipulation

# Load bird data
birddata <- read.delim("path/to/birddata.txt)

# Clean bird data
birddata <- birddata %>%
  filter(CATEGORY %in% c("species", "issf"))

birddata <- birddata [, -c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]
birddata$OBSERVATION.COUNT <- gsub("x", "1", birddata$OBSERVATION.COUNT)
birddata <- birddata %>% replace_na(list(OBSERVATION.COUNT = 1))
birddata$OBSERVATION.COUNT <- as.numeric(trimws(birddata$OBSERVATION.COUNT))

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

# Create abundance matrix with explicit column naming
abundance_matrix <- birddata_hex %>%
  st_drop_geometry() %>%
  
  group_by(ID = ID.x, SCIENTIFIC.NAME) %>%
  summarise(n = n(), .groups = "drop") %>%

  pivot_wider(
    names_from = SCIENTIFIC.NAME,
    values_from = n,
    values_fill = 0
  ) %>%

  mutate(across(-ID, as.numeric)) %>%
  column_to_rownames("ID")

# Calculate biodiversity indices
indices_df <- data.frame(
  ID = rownames(abundance_matrix),
  richness = vegan::specnumber(abundance_matrix),
  shannon = vegan::diversity(abundance_matrix),
  simpson = vegan::diversity(abundance_matrix, index = "simpson"),
  stringsAsFactors = FALSE
)

# Merge with hexagon grid
hex_grid_with_indices <- hex_grid %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    indices_df %>% mutate(ID = as.character(ID)),
    by = "ID"
  ) %>%
  mutate(across(c(richness, shannon, simpson), ~ if_else(is.na(.), 0, .)))

# Calculate Gini Index

gini_por_hexagono <- citydata %>%
  st_join(hex_grid) %>%         
  st_drop_geometry() %>%        
  group_by(ID) %>%              
  summarise(
    gini = ineq::Gini(V_REF),   
    mean_valor = mean(V_REF),    
    sd_valor = sd(V_REF)        
  )
  
# Intersect Gini data  
  
hex_grid_with_gini <- hex_grid %>%
  left_join(gini_por_hexagono, by = "ID") %>%
  mutate(gini = replace_na(gini, 0)
  
