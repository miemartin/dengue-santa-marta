# =====================================================
# Script: 03_lag_calculation.R
# Purpose: Creation of lagged environmental and
#          demographic predictors for dengue analysis
# =====================================================

# -------------------------
# Libraries
# -------------------------
library(dplyr)

# -------------------------
# Load dataset
# -------------------------
data <- read.csv("path/to/data.csv", sep = ";")

# -------------------------
# Identify variables to lag
# (exclude time and outcome variables)
# -------------------------
vars_to_lag <- setdiff(
  names(data),
  c("year", "week", "date", "dengue")
)

# -------------------------
# Create lagged variables (t-1 to t-4)
# -------------------------
max_lag <- 4

data_lags <- data %>%
  arrange(year, week) %>%
  mutate(
    across(
      all_of(vars_to_lag),
      list(
        lag1 = ~ lag(.x, 1),
        lag2 = ~ lag(.x, 2),
        lag3 = ~ lag(.x, 3),
        lag4 = ~ lag(.x, 4)
      ),
      .names = "{.col}_{.fn}"
    )
  )

# -------------------------
# Save dataset with lags
# -------------------------
write.csv(
  data_lags,
  "data_lags.csv",
  row.names = FALSE
)

# -------------------------
# Quick check
# -------------------------
head(data_lags)
