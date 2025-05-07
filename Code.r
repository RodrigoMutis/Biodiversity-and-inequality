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
birddata <- read.delim("data/Bogotabirddata.txt") |> 
    #janitor::clean_names() |> 
    as_tibble()

birddata

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
citydata <- st_read("data/Bogota Value/Valor_ref_2023.shp") # Polygons

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
citydata_hex <- # st_intersection(citydata_utm, hex_grid) %>% 
  st_join(citydata_utm,hex_grid, by = "ID")
  
birddata_hex <- st_intersection(birddata_utm, hex_grid) #%>% 
  #st_join(hex_grid, birddata_utm, by = "ID")

birddata_hex |> 
    as.data.frame() |> 
    select(-geometry) |> 
    group_by(ID, SCIENTIFIC.NAME) |> 
    summarize(obs = sum(OBSERVATION.COUNT)) |> 
    right_join(hex_grid)

birdat <- aggregate(x = birddata_utm, by = hex_grid, FUN = sum)
  
merged_data <- st_join(citydata_hex, birddata_hex, by = "ID")

# Calculate biodiversity index 

# Classify species into common/rare (top 10% / middle 80%)
species_abundance <- birddata_hex %>%
  st_drop_geometry() %>%
  group_by(SCIENTIFIC.NAME) %>%
  summarise(total_obs = sum(OBSERVATION.COUNT)) %>%
  mutate(
    pct_rank = percent_rank(total_obs),
    group = case_when(
      pct_rank <= 0.10 ~ "common",
      pct_rank >= 0.90 ~ "error",
      TRUE ~ "rare"
    )
  )

common_species <- filter(species_abundance, group == "common")$SCIENTIFIC.NAME
rare_species <- filter(species_abundance, group == "rare")$SCIENTIFIC.NAME

# Calculate biodiversity indices for each group
calculate_indices <- function(data, suffix) {
  if(nrow(data) == 0) return(data.frame(ID = character()))
  
  data %>%
    st_drop_geometry() %>%
    group_by(ID, SCIENTIFIC.NAME) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = SCIENTIFIC.NAME, values_from = n, values_fill = 0) %>%
    column_to_rownames("ID") %>%
    {
      data.frame(
        ID = rownames(.),
        richness = specnumber(.),
        shannon = diversity(.),
        simpson = diversity(., index = "simpson"),
        stringsAsFactors = FALSE
      )
    } %>%
    rename_with(~ paste0(., "_", suffix), -ID)
}

indices_common <- calculate_indices(
  filter(birddata_hex, SCIENTIFIC.NAME %in% common_species),
  "common"
)

indices_rare <- calculate_indices(
  filter(birddata_hex, SCIENTIFIC.NAME %in% rare_species),
  "rare"
)
# Calculate gini index by hexagon
gini_hex <- citydata_utm %>%
  st_join(hex_grid) %>%
  st_drop_geometry() %>%
  group_by(ID) %>%
  summarise(
    gini = Gini(V_REF),
    mean_valor = mean(V_REF),
    sd_valor = sd(V_REF)
  )
# Create final data set 
final_data <- hex_grid %>%
  left_join(indices_common, by = "ID") %>%
  left_join(indices_rare, by = "ID") %>%
  left_join(gini_hex, by = "ID") %>%
  mutate(across(starts_with(c("richness", "shannon", "simpson")), ~ replace_na(., 0))) %>%
  filter(!is.na(gini)) %>%
  st_as_sf()

# Spatial regression 
run_spatial_models <- function(data, response) {
  coords <- st_centroid(data) %>% st_coordinates()
  nb <- knn2nb(knearneigh(coords, k = 5))
  lw <- nb2listw(nb, style = "W")
  
  formula <- as.formula(paste(response, "~ gini + mean_valor"))
  
  list(
    sar = tryCatch(lagsarlm(formula, data, lw, method = "eigen"), error = function(e) NULL),
    sem = tryCatch(errorsarlm(formula, data, lw, method = "eigen"), error = function(e) NULL),
    sdm = tryCatch(lagsarlm(formula, data, lw, type = "mixed", method = "eigen"), error = function(e) NULL)
  )
}

# Run models
metrics <- c("richness", "shannon", "simpson")
results <- list()

for(metric in metrics) {
  for(group in c("common", "rare")) {
    response <- paste0(metric, "_", group)
    
    if(response %in% names(final_data)) {
      filtered_data <- final_data %>% 
        filter(!is.na(get(response)))
      
      if(nrow(filtered_data) > 0) {
        results[[paste(group, metric)]] <- run_spatial_models(filtered_data, response)
      }
    }
  }
}

# Extract
extract_results <- function(models) {
  purrr::map_dfr(names(models), function(model_type) {
    model <- models[[model_type]]
 tryCatch({
      broom::tidy(model) %>%
        dplyr::mutate(
          model = model_type,
          aic = if("logLik" %in% methods(class=class(model)[1])) AIC(model) else NA_real_,
          bic = if("logLik" %in% methods(class=class(model)[1])) BIC(model) else NA_real_
        )
    }, error = function(e) {
      warning("Failed to process ", model_type, ": ", e$message)
      NULL
    })
  })
}    

all_results <- purrr::map_dfr(
  results,
  ~tryCatch({
    res <- extract_results(.)
    if(is.null(res)) {
      tibble::tibble(note = "No valid models in this group")
    } else {
      res
    }
  }, error = function(e) {
    tibble::tibble(note = paste("Processing failed:", e$message))
  }),
  .id = "analysis"
)

# Summary
coefficient_table <- all_results %>%
  select(analysis, model, term, estimate, std.error, p.value, aic, bic) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    analysis = factor(analysis, 
                      levels = c("common richness", "common shannon", "common simpson",
                                 "rare richness", "rare shannon", "rare simpson"))
  )

aic_table <- all_results %>%
  distinct(analysis, model, aic, bic) %>%
  pivot_wider(names_from = model, values_from = c(aic, bic))

# Visualization
plot_coefficients <- function(df, var) {
  df %>%
    filter(term == var) %>%
    ggplot(aes(x = analysis, y = estimate, color = model)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = estimate - 1.96*std.error,
          ymax = estimate + 1.96*std.error),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("Coefficients for", var),
         x = "Analysis Group", y = "Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

gini_plot <- plot_coefficients(coefficient_table, "gini")
valor_plot <- plot_coefficients(coefficient_table, "mean_valor")

# Results
list(
  coefficients = coefficient_table,
  aic_comparison = aic_table,
  plots = list(gini = gini_plot, valor = valor_plot)
)
