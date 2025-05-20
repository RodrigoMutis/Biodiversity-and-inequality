# Load libraries
library(sf)          # For geospatial data
library(tidyverse)       # For data manipulation
library(spatialreg)  # For Spatial regression
library(spdep)       # For spatial dependence
library(tictoc)
library(naniar)
library(fs)

fls <- dir_ls("data/500/") |> str_subset(pattern = "\\.Rda$")

load(fls[1])
bta <- dat
load(fls[2])
cape <- dat
load(fls[3])
ldn <- dat
load(fls[4])
la <- dat
load(fls[5])
mxc <- dat
load(fls[6])
stg <- dat
load(fls[7])
syd <- dat
rm(dat)

dat <- list(bta, cape, ldn, la, mxc, stg, syd)

dat[[3]] |> 
    ggplot() +
    geom_sf(aes(fill = gini))

dat[[7]] |> 
    #filter(shannon >0) |> 
    ggplot(aes(gini, mean_value)) + 
    geom_point() + geom_smooth(method = "lm") +
    scale_y_log10()

dat |> 
    ggplot() +
    geom_sf(aes(fill = gini))

out <- dat |> 
    map(st_drop_geometry) |> 
    map2(.y = fls,
         .f = function(x,y) {x$city <- y; return(x)} ) |> 
    bind_rows() |>
    as_tibble() |> 
    mutate(city = str_remove_all(city, "data/|\\.Rda$"))

out |> ggplot(aes(gini, shannon)) +
    geom_point(aes(color = city), alpha = 0.5) +
    geom_smooth(aes(color = city))#+ scale_x_log10()


#### Regressions ####

list_weights <- function(x){
    coords <- st_centroid(x) %>% st_coordinates()
    nb <- knn2nb(knearneigh(coords, k = 5))
    lw <- nb2listw(nb, style = "W")
    return(lw)
}

sar <- safely(function(x){
    lw <- list_weights(x)
    fit <- lagsarlm(shannon ~ gini  , x, lw, method = "eigen")
    return(fit)
})

sem <- safely(function(x){
    lw <- list_weights(x)
    fit <- errorsarlm(shannon ~ gini , x, lw, method = "eigen")
    return(fit)
})

sdm <- safely(function(x){
    lw <- list_weights(x)
    fit <- lagsarlm(shannon ~ gini, x, lw, method = "eigen", type = "mixed")
    return(fit)
})

tic()
sar_fits <- map(dat, sar)
toc() # 43s

tic()
sem_fits <- map(dat, sem)
toc() # 32

tic()
sdm_fits <- map(dat, sdm)
toc() # 40s

# sar_fits <- transpose(sar_fits)
# not_ok <- map_lgl(transpose(sem_fits)$error, is.null) ## all models fail

sar_df <- transpose(sar_fits)$result |> map(.f = broom::tidy ) |> 
    map2(.y = str_remove_all(fls, "data/|\\.Rda$"), 
         .f = function(x,y){x$city <- y; return(x)}) |> 
    map(.f = function(x) {x$model <- "sar"; return(x)})

sem_df <- transpose(sem_fits)$result |> map(.f = broom::tidy ) |> 
    map2(.y = str_remove_all(fls, "data/|\\.Rda$"), 
         .f = function(x,y){x$city <- y; return(x)})|> 
    map(.f = function(x) {x$model <- "sem"; return(x)})

sdm_df <- transpose(sdm_fits)$result |> map(.f = broom::tidy ) |> 
    map2(.y = str_remove_all(fls, "data/|\\.Rda$"), 
         .f = function(x,y){x$city <- y; return(x)})|> 
    map(.f = function(x) {x$model <- "sdm"; return(x)})

df_fit <- bind_rows(bind_rows(sar_df), bind_rows(sdm_df) , bind_rows(sem_df))

df_fit |> 
    #filter(term != "rho", term != "lambda") |> 
    mutate(p_value = case_when(
        p.value > 0.05 ~ "p > 0.05",
        p.value <= 0.05 & p.value > 0.01 ~ "p < 0.05",
        p.value <= 0.01 ~ "p < 0.01"
    )) |> 
    ggplot(aes(estimate, term)) + 
    geom_point(aes(col = p_value)) +
    geom_errorbarh(aes(xmin = estimate-std.error, xmax = estimate+std.error, color = p_value),
                   height = 0.25) +
    geom_vline(xintercept = 0, linetype = 2, color = "black") +
    facet_grid(model ~ city, scales = "free") +
    theme_light(base_size = 10)

save(sdm_fits, sar_fits, sem_fits, df_fit, file = "results_regressions.Rda")

#### Left over ####
coords <- st_centroid(mxc) %>% st_coordinates()
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")

# resp <- names(dat)[5:7]
# 
# df_int <- interaction(resp, c("lagsarlm", "errorsarlm", "lagsarlm"), sep = "_") |> 
#     levels() |> 
#     as_tibble() |> 
#     separate(value, into = c("resp", "model"), sep = "_", remove = TRUE)

tic()
fit <- lagsarlm(shannon ~ gini , mxc, lw, method = "eigen")
toc() #3.04 s

tic()
fit2 <- errorsarlm(shannon ~ gini + mean_value, bogota, lw, method = "eigen")
toc() #3.3s

tic()
fit3 <- lagsarlm(shannon ~ gini + mean_value, dat, lw, method = "eigen", type = "mixed")
toc() #3s

summary(fit); summary(fit2); summary(fit3)


fits <- map2(
    .x = df_int$resp,
    .y = df_int$model,
    .f = function(x,y){
        
    }
)