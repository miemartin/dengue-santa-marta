# =====================================================
# Script: 04_statistical_modelling.R
# Purpose: Statistical modelling of dengue cases using
#          climatic, environmental, and demographic predictors
# Study area: Santa Marta, Colombia
# =====================================================

# -------------------------
# Libraries
# -------------------------
library(tidyverse)
library(psych)
library(lme4)
library(MASS)
library(car)
library(DHARMa)
library(performance)
library(ggeffects)
library(gridExtra)

# -------------------------
# Data loading
# -------------------------
data_l <- read.csv("path/to/data_lags.csv", sep = ";")

# -------------------------
# Data preparation
# -------------------------
data_l <- data_l %>%
  mutate(
    year = factor(year),
    week = factor(week)
  )

# Remove rows with missing values introduced by lagging
dl <- na.omit(data_l)

# -------------------------
# Exploratory analysis
# -------------------------
summary(dl)
describe(dl)

# Temporal patterns
ggplot(dl, aes(x = week, y = dengue)) +
  geom_bar(stat = "identity", fill = "grey70") +
  theme_classic() +
  labs(x = "Week", y = "Dengue cases")

ggplot(dl, aes(x = year, y = dengue)) +
  geom_bar(stat = "identity", fill = "grey70") +
  theme_classic() +
  labs(x = "Year", y = "Dengue cases")

# -------------------------
# Correlation analysis (Spearman)
# -------------------------
cor_matrix <- cor(
  dl %>% select(-dengue, -year, -week),
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(cor_matrix, "correlation_spearman.csv")

# -------------------------
# Standardization of predictors
# (only variables considered for modelling)
# -------------------------
dl <- dl %>%
  mutate(across(
    c(popden, NDVI_min, NDVI, EVI_max,
      LSTN_max, LSTN_mode,
      NDVI_min_1, NDVI_mode_2,
      EVI_max_3, EVI_sum_1,
      LSTD_min_1, LSTD_max_1, LSTD_sd_3,
      LSTN_min_4, LSTN_mode_3, LSTN_sum_4,
      NDWI_mode_1, pp_min_3, pp_4),
    scale
  ))

# -------------------------
# Model comparison: Poisson vs Negative Binomial
# -------------------------
glm_pois <- glm(
  dengue ~ popden + NDVI + LSTN_max + LSTD_max_1 +
    LSTD_sd_3 + LSTN_min_4 + NDWI_mode_1 + pp_4,
  family = poisson(),
  data = dl
)

glmm_pois <- glmer(
  dengue ~ popden + NDVI + LSTN_max + LSTD_max_1 +
    LSTD_sd_3 + LSTN_min_4 + NDWI_mode_1 + pp_4 +
    (1 | year),
  family = poisson,
  data = dl,
  control = glmerControl(optimizer = "bobyqa")
)

glmm_nb <- glmer.nb(
  dengue ~ popden + NDVI + LSTN_max + LSTD_max_1 +
    LSTD_sd_3 + LSTN_min_4 + NDWI_mode_1 + pp_4 +
    (1 | year),
  data = dl
)

AIC(glm_pois, glmm_pois, glmm_nb)

# -------------------------
# Final main-effects model
# -------------------------
final_model <- glmer.nb(
  dengue ~ popden + NDVI + LSTD_max_1 +
    LSTN_min_4 + NDWI_mode_1 + pp_4 +
    (1 | year),
  data = dl
)

summary(final_model)

# -------------------------
# Multicollinearity check
# -------------------------
vif(
  lm(
    dengue ~ popden + NDVI + LSTD_max_1 +
      LSTN_min_4 + NDWI_mode_1 + pp_4,
    data = dl
  )
)

# -------------------------
# Interaction models (biologically motivated)
# -------------------------
interaction_combined <- glmer.nb(
  dengue ~ 
    popden * NDVI +
    NDVI * LSTD_max_1 +
    LSTN_min_4 +
    NDWI_mode_1 +
    pp_4 +
    (1 | year),
  data = dl,
  control = glmerControl(optimizer = "bobyqa")
)

null_model <- glmer.nb(
  dengue ~ 1 + (1 | year),
  data = dl,
  na.action = na.fail
)

AIC(final_model, interaction_combined, null_model)
summary(interaction_combined)

# -------------------------
# Model diagnostics
# -------------------------
sim_res <- simulateResiduals(interaction_combined, n = 1000)

plot(sim_res)
acf(sim_res$scaledResiduals)

testDispersion(sim_res)
testZeroInflation(sim_res)

# -------------------------
# Model performance
# -------------------------
r2_nakagawa(interaction_combined)
r2_nakagawa(null_model)

# -------------------------
# Observed vs predicted (annual aggregation)
# -------------------------
dl$predicted <- predict(final_model, type = "response", re.form = NA)

dl_year <- dl %>%
  group_by(year) %>%
  summarise(
    observed = sum(dengue),
    predicted = sum(predicted),
    .groups = "drop"
  )

ggplot(dl_year, aes(x = year)) +
  geom_line(aes(y = observed, color = "Observed", group = 1)) +
  geom_line(aes(y = predicted, color = "Predicted", group = 1)) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_classic() +
  labs(x = "Year", y = "Dengue cases", color = NULL)

# -------------------------
# Marginal effects and interactions
# -------------------------
p_NDVI <- ggpredict(interaction_combined, terms = "NDVI")
p_LSTDmax <- ggpredict(interaction_combined, terms = "LSTD_max_1")
p_LSTNmin <- ggpredict(interaction_combined, terms = "LSTN_min_4")
p_PP <- ggpredict(interaction_combined, terms = "pp_4")

p_NDVI_pop <- ggpredict(
  interaction_combined,
  terms = c("NDVI", "popden [-1, 1]")
)

p_NDVI_LSTD <- ggpredict(
  interaction_combined,
  terms = c("NDVI", "LSTD_max_1 [-1, 1]")
)

plot_effect <- function(pred, xlab) {
  ggplot(pred, aes(x = x, y = predicted)) +
    geom_line(size = 1.1) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.2, fill = "grey70"
    ) +
    labs(x = xlab, y = "Predicted dengue cases") +
    theme_classic()
}

gA <- plot_effect(p_NDVI, "NDVI")
gB <- plot_effect(p_LSTDmax, "LSTD max (lag 1)")
gC <- plot_effect(p_LSTNmin, "LSTN min (lag 4)")
gD <- plot_effect(p_PP, "Precipitation (lag 4)")

gE <- ggplot(p_NDVI_pop,
             aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2, color = NA
  ) +
  theme_classic() +
  labs(x = "NDVI", y = "Predicted dengue cases",
       color = "Population density",
       fill = "Population density")

gF <- ggplot(p_NDVI_LSTD,
             aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2, color = NA
  ) +
  theme_classic() +
  labs(x = "NDVI", y = "Predicted dengue cases",
       color = "LSTD max (lag 1)",
       fill  = "LSTD max (lag 1)")

figure3 <- grid.arrange(gA, gB, gC, gD, gE, gF, ncol = 3)

ggsave(
  "Figure3.tif",
  figure3,
  width = 12,
  height = 8,
  dpi = 600
)
