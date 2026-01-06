# =====================================================
# Script: 02_population_density.R
# Purpose: Estimation and temporal interpolation of
#          population density for Santa Marta, Colombia
# Data: Raster-based population density layers
# =====================================================

# -------------------------
# Libraries
# -------------------------
library(terra)
library(sf)
library(exactextractr)
library(purrr)
library(dplyr)

# -------------------------
# Paths (data not public)
# -------------------------
raster_folder <- "path/to/raster/"      # Folder with population density rasters
shp_path <- "path/to/MGN_URB_AREA_CENSAL.shp"

# -------------------------
# Load spatial data
# -------------------------
santa_marta <- st_read(shp_path, quiet = TRUE)

raster_files <- list.files(
  raster_folder,
  pattern = "\\.tif$",
  full.names = TRUE
)

# -------------------------
# Function to process one raster
# -------------------------
process_raster <- function(file, polygon) {
  
  r <- rast(file)
  r_masked <- mask(crop(r, polygon), polygon)
  
  pop_density <- exact_extract(r_masked, polygon, fun = "mean")
  
  # Extract year safely (assumes YYYY in filename)
  year <- as.numeric(stringr::str_extract(basename(file), "\\d{4}"))
  
  tibble(
    year = year,
    pop_density = pop_density
  )
}

# -------------------------
# Apply to all rasters
# -------------------------
df_results <- map_dfr(
  raster_files,
  process_raster,
  polygon = santa_marta
) %>%
  arrange(year)

write.csv(
  df_results,
  "population_density_santa_marta_observed.csv",
  row.names = FALSE
)

# -------------------------
# Temporal interpolation / extrapolation
# -------------------------
model <- lm(pop_density ~ poly(year, 2, raw = TRUE), data = df_results)

years_full <- seq(min(df_results$year), max(df_results$year), by = 1)

df_pred <- tibble(
  year = years_full,
  pop_density = predict(model, newdata = data.frame(year = years_full))
)

write.csv(
  df_pred,
  "population_density_santa_marta_interpolated.csv",
  row.names = FALSE
)

# -------------------------
# Visualization
# -------------------------
plot(
  df_results$year,
  df_results$pop_density,
  pch = 16,
  col = "blue",
  xlab = "Year",
  ylab = "Population density",
  main = "Observed and interpolated population density"
)

lines(
  df_pred$year,
  df_pred$pop_density,
  col = "red",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Observed", "Interpolated"),
  col = c("blue", "red"),
  pch = c(16, NA),
  lwd = c(NA, 2),
  bty = "n"
)
