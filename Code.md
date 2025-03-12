# Biodiversity-and-inequality
# Install necessary packages 
install.packages("sf")         
install.packages("ineq")
install.packages("vegan")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggpubr")  

# Load libraries
library(sf)          # For geospatial data
library(ineq)        # For Gini Index
library(vegan)       # For biodiversity indices
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(ggpubr)      # For correlation plots

# Load the shapefiles
citydata <- st_read("path/to/citydata.shp") # Polygons
birddata <- st_read("path/to/birddata.shp") # Points

# Check CRS of both shapefiles
st_crs(citydata)
st_crs(birddata)

# If CRS is different, reproject one of them
birddata <- st_transform(birddata, st_crs(citydata))

# Perform the spatial join
points_with_values <- st_join(birddata, citydata, join = st_within)

# Extract the value column
values <- points_with_values$value

# Calculate the Gini Index for each hexagon
gini_index <- points_with_values %>%
  group_by(id) %>%  # Assuming 'id' is the unique identifier for each hexagon
  summarise(gini = Gini(value), .groups = 'drop')

# Load or create your georeferenced data
study_area <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10)))
study_area <- st_set_crs(study_area, 4326)  # Set CRS (e.g., WGS84)

# Generate hexagonal tessellation
hex_grid <- st_make_grid(study_area, cellsize = 1, square = FALSE)

# Convert the grid to an sf object for easier handling
hex_grid <- st_sf(geometry = hex_grid)

# Add a unique ID to each hexagon
hex_grid$id <- 1:nrow(hex_grid)

# Spatially join bird points with the hexagonal grid
bird_in_hex <- st_join(birddata, hex_grid, join = st_within)

# Aggregate species data within each hexagon
species_counts <- bird_in_hex %>%
  group_by(id, species) %>%  # Replace 'species' with the column name for species in your data
  summarise(count = n(), .groups = 'drop')

# Reshape data to a wide format (one row per hexagon, one column per species)
species_wide <- species_counts %>%
  tidyr::pivot_wider(names_from = species, values_from = count, values_fill = 0)

# Calculate biodiversity indices for each hexagon
biodiversity_indices <- species_wide %>%
  rowwise() %>%
  mutate(
    shannon = diversity(c_across(-id)),  # Shannon-Wiener Index
    simpson = diversity(c_across(-id)),  # Simpson's Index
    richness = sum(c_across(-id) > 0)    # Species Richness
  ) %>%
  ungroup()

# Merge biodiversity indices and Gini Index into the hexagon grid
hex_grid_with_indices <- hex_grid %>%
  left_join(biodiversity_indices, by = "id") %>%
  left_join(gini_index, by = "id")

# Calculate correlation coefficients
cor_shannon_gini <- cor(hex_grid_with_indices$shannon, hex_grid_with_indices$gini, method = "pearson")
cor_simpson_gini <- cor(hex_grid_with_indices$simpson, hex_grid_with_indices$gini, method = "pearson")
cor_richness_gini <- cor(hex_grid_with_indices$richness, hex_grid_with_indices$gini, method = "pearson")

# Print correlation coefficients
print(paste("Correlation (Shannon vs Gini):", cor_shannon_gini))
print(paste("Correlation (Simpson vs Gini):", cor_simpson_gini))
print(paste("Correlation (Richness vs Gini):", cor_richness_gini))

# Visualize correlations using scatter plots
plot_shannon_gini <- ggplot(hex_grid_with_indices, aes(x = gini, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Shannon-Wiener Index vs Gini Index",
       x = "Gini Index",
       y = "Shannon-Wiener Index") +
  theme_minimal()

plot_simpson_gini <- ggplot(hex_grid_with_indices, aes(x = gini, y = simpson)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Simpson's Index vs Gini Index",
       x = "Gini Index",
       y = "Simpson's Index") +
  theme_minimal()

plot_richness_gini <- ggplot(hex_grid_with_indices, aes(x = gini, y = richness)) +
  geom_point() +
  geom_smooth(method = "lm", col = "green") +
  labs(title = "Species Richness vs Gini Index",
       x = "Gini Index",
       y = "Species Richness") +
  theme_minimal()

# Arrange plots in a single figure
combined_plots <- ggarrange(plot_shannon_gini, plot_simpson_gini, plot_richness_gini, ncol = 2, nrow = 2)
print(combined_plots)
