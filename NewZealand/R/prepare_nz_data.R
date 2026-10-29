# Data preparation shared by the two New Zealand model scripts
# (nz_fit_detection_only.R and nz_fit_detection_nondetection.R), so that both
# models are fitted to exactly the same data.
#
# prepare_nz_data(species_id) returns a list with
#   y           0/1 record of the species at each presence-absence site
#   x           covariate matrix (sites x 11), in the order of beta1, ..., beta11
#   covariates  names of the columns of x

library(disdat)
library(dplyr)
library(mgcv)    # GAM for latitude adjustment
library(sf)      # coordinate transformation

prepare_nz_data <- function(species_id, data_crs = 3857) {

  # Load the disdat data 
  pa  <- disPa("NZ")    # presence-absence records (one column per species)
  env <- disEnv("NZ")   # environmental covariates at the presence-absence sites
  po  <- disPo("NZ")    # presence-only records
  bg  <- disBg("NZ")    # background (pseudo-absence) sites

  stopifnot(species_id %in% names(pa), nrow(pa) == nrow(env))

  env$toxicats <- NULL
  names(env)[3:4] <- c("xcoord", "ycoord")   # columns 3 and 4 hold the site coordinates

  # Standardise covariates 
  covars_scaled <- c("deficit", "dem", "hillshade", "mas", "rain", "slope", "sseas", "tseas")
  covars  <- rbind(po[, covars_scaled], bg[, covars_scaled])

  for (v in covars_scaled)
    env[[v]] <- (env[[v]] - mean(covars[[v]])) / sd(covars[[v]])

  # ---- Site-level data for the chosen species -------------------------------------
  sites <- env |>
    mutate(
      y    = as.integer(pa[[species_id]]),   # 0/1 record
      age0 = as.integer(age == 0),           # age class dummies (class 2 is the reference)
      age1 = as.integer(age == 1)
    )

  # Latitude (degrees) from the projected coordinates, then standardised
  sites_ll <- sites |>
    st_as_sf(coords = c("xcoord", "ycoord"), crs = data_crs) |>
    st_transform(4326)
  sites$latitude        <- st_coordinates(sites_ll)[, 2]
  sites$latitude_scaled <- as.numeric(scale(sites$latitude))

  # Replace climate covariates by their residuals after removing a smooth latitudinal trend
  for (v in c("rain", "tseas", "sseas", "mas"))
    sites[[paste0(v, "_resid")]] <-
      residuals(gam(reformulate("s(latitude_scaled)", response = v), data = sites))

  # Covariates in the order of the model coefficients beta1, ..., beta11
  covariates <- c("age0", "age1", "deficit", "dem", "hillshade",
                  "mas_resid", "rain_resid", "slope", "sseas_resid", "tseas_resid",
                  "latitude_scaled")

  x <- as.matrix(sites[, covariates])
  stopifnot(!anyNA(x), all(sites$y %in% 0:1))

  list(y = sites$y, x = x, covariates = covariates)
}
