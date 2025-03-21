# Biodiversity-and-inequality

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
birddata$OBSERVATION.COUNT <- as.numeric(birddata$OBSERVATION.COUNT)
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
cellsize <- 100
hex_grid <- st_make_grid(citydata_utm, cellsize = cellsize, square = FALSE) %>% 
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
  
merged_data <- merge(citydata_hex, birddata_hex, by = "ID")

# Calculate biodiversityindex in each hexagon
