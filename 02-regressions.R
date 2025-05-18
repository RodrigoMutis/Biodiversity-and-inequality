# Load libraries
library(sf)          # For geospatial data
library(tidyverse)       # For data manipulation
library(spatialreg)  # For Spatial regression
library(spdep)       # For spatial dependence
library(tictoc)
library(naniar)


load("data/bogota.Rda")

#### Regressions ####
coords <- st_centroid(dat) %>% st_coordinates()
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

resp <- names(dat)[5:7]

df_int <- interaction(resp, c("lagsarlm", "errorsarlm", "lagsarlm"), sep = "_") |> 
    levels() |> 
    as_tibble() |> 
    separate(value, into = c("resp", "model"), sep = "_", remove = TRUE)

tic()
fit <- lagsarlm(shannon ~ gini + mean_value, dat, lw, method = "eigen")
toc() #3.04 s

tic()
fit2 <- errorsarlm(shannon ~ gini + mean_value, dat, lw, method = "eigen")
toc() #3.3s

tic()
fit3 <- lagsarlm(shannon ~ gini + mean_value, dat, lw, method = "eigen")
toc() #3s

summary(fit2)


fits <- map2(
    .x = df_int$resp,
    .y = df_int$model,
    .f = function(x,y){
        
    }
)