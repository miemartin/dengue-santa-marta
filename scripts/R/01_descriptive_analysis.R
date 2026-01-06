# =====================================================
# Script: 01_descriptive_analysis.R
# Purpose: Descriptive analysis and visualization of dengue
#          cases and environmental variables
# Study area: Santa Marta, Colombia
# =====================================================

# -------------------------
# Libraries
# -------------------------
library(tidyverse)
library(scales)

# -------------------------
# Monthly average of relative humidity
# -------------------------
rh <- read.csv("path/to/rh.csv", sep = ";")

rh_av <- rh %>%
  group_by(year, month) %>%
  summarise(rh_av = mean(rh, na.rm = TRUE), .groups = "drop")

write.csv(rh_av, "rh_av.csv", row.names = FALSE)

# -------------------------
# Load main dataset
# -------------------------
data <- read.csv("path/to/data_lags.csv", sep = ";") %>%
  mutate(
    year_week = factor(
      paste(year, week, sep = "-"),
      levels = unique(paste(year, week, sep = "-"))
    )
  ) %>%
  na.omit()

# -------------------------
# Helper function for dual-axis plots
# -------------------------
plot_dual_axis <- function(data, x, y_bar, y_line, line_max, y_line_label) {
  
  scale_factor <- max(data[[y_bar]], na.rm = TRUE) / line_max
  
  ggplot(data, aes(x = .data[[x]])) +
    geom_bar(aes(y = .data[[y_bar]]),
             stat = "identity", fill = "#959494") +
    geom_line(aes(y = .data[[y_line]] * scale_factor, group = 1),
              color = "red", size = 1) +
    scale_y_continuous(
      name = "Dengue Cases",
      sec.axis = sec_axis(~ . / scale_factor, name = y_line_label)
    ) +
    scale_x_discrete(breaks = levels(data[[x]])[seq(1, nrow(data), by = 40)]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Year–Epidemiological Week")
}

# -------------------------
# Population density
# -------------------------
popden_plot <- plot_dual_axis(
  data = data,
  x = "year_week",
  y_bar = "dengue",
  y_line = "popden",
  line_max = 350,
  y_line_label = "Population density (per km²)"
)

ggsave("pop_density_plot.tiff", popden_plot)

# -------------------------
# NDVI mode
# -------------------------
ndvi_plot <- plot_dual_axis(
  data = data,
  x = "year_week",
  y_bar = "dengue",
  y_line = "NDVI_mode",
  line_max = 0.7,
  y_line_label = "NDVI mode"
)

ggsave("NDVImode.tiff", ndvi_plot)

# -------------------------
# Maximum LST day
# -------------------------
lstd_plot <- plot_dual_axis(
  data = data,
  x = "year_week",
  y_bar = "dengue",
  y_line = "LSTD_max",
  line_max = 50,
  y_line_label = "Maximum LST day (°C)"
)

# -------------------------
# NDWI mode (custom scaling)
# -------------------------
ndwi_min <- -0.75
ndwi_max <- 0

scale_factor_ndwi <- max(data$dengue, na.rm = TRUE) / abs(ndwi_min)

ndwi_plot <- ggplot(data, aes(x = year_week)) +
  geom_bar(aes(y = dengue), stat = "identity", fill = "#959494") +
  geom_line(
    aes(y = (NDWI_mode - ndwi_min) * scale_factor_ndwi, group = 1),
    color = "red", size = 1
  ) +
  scale_y_continuous(
    name = "Dengue Cases",
    sec.axis = sec_axis(
      ~ . / scale_factor_ndwi + ndwi_min,
      name = "NDWI mode",
      breaks = seq(ndwi_min, ndwi_max, by = 0.1)
    )
  ) +
  scale_x_discrete(breaks = levels(data$year_week)[seq(1, nrow(data), by = 40)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year–Epidemiological Week")

# -------------------------
# Precipitation variability (SD)
# -------------------------
pp_plot <- plot_dual_axis(
  data = data,
  x = "year_week",
  y_bar = "dengue",
  y_line = "pp_sd",
  line_max = 95,
  y_line_label = "Precipitation SD (mm)"
)

# -------------------------
# Descriptive statistics
# -------------------------
dl <- data %>%
  mutate(
    year = factor(year),
    week = factor(week)
  )

summary(dl)

desc_stats <- psych::describe(dl)
write.csv(desc_stats, "descriptive_statistics.csv")
