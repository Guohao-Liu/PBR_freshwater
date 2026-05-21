
library(rstudioapi)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(readr)
library(readxl)
library(tidyr)
library(purrr)
library(rlang)

script_path <- getActiveDocumentContext()$path
script_dir  <- dirname(script_path)
setwd(script_dir)


z_within_group <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    rep(0, length(x))
  } else {
    as.numeric(scale(x))
  }
}

dat <- read_xlsx(file.path("..", "000_input_file", "All_data.xlsx")) %>%
  mutate(
    Richness   = as.numeric(Richness),
    FRic       = as.numeric(FRic),
    TP         = as.numeric(TP),
    logTP      = log(TP + 1),
    bio1_1500  = as.numeric(bio1_1500),
    bio12_1500 = as.numeric(bio12_1500),
    AL         = readr::parse_number(as.character(AL))
  ) %>%
  filter(
    !is.na(Richness),
    !is.na(FRic),
    !is.na(TP),
    !is.na(logTP),
    !is.na(Trophic_level),
    !is.na(Ecosystem),
    !is.na(Dataset),
    !is.na(bio1_1500),
    !is.na(bio12_1500),
    !is.na(bio4_1500),
    !is.na(bio15_1500),
    !is.na(Elevation),
    !is.na(AL)
  ) %>%
  group_by(Dataset) %>%
  mutate(
    Richness_z = z_within_group(Richness),
    FRic_z     = z_within_group(FRic)
  ) %>%
  ungroup() %>%
  mutate(
    Dataset       = as.factor(Dataset),
    Trophic_level = as.factor(Trophic_level),
    Ecosystem     = as.factor(Ecosystem)
  )


# Step 1
m0 <- lmer(Richness_z ~ logTP + I(logTP^2) + (1 + logTP || Dataset), REML = FALSE, data = dat)

m_land <- lmer(Richness_z ~ logTP + I(logTP^2) + Landuse_class + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_bio1 <- lmer(Richness_z ~ logTP + I(logTP^2) + bio1_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_bio12 <- lmer(Richness_z ~ logTP + I(logTP^2) + bio12_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_bio4 <- lmer(Richness_z ~ logTP + I(logTP^2) + bio4_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_bio15 <- lmer(Richness_z ~ logTP + I(logTP^2) + bio15_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_elev <- lmer(Richness_z ~ logTP + I(logTP^2) + Elevation + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_trophic <- lmer(Richness_z ~ logTP + I(logTP^2) + Trophic_level + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_ecosystem <- lmer(Richness_z ~ logTP + I(logTP^2) + Ecosystem + (1 + logTP || Dataset), REML = FALSE, data = dat)
m_lat <- lmer(Richness_z ~ logTP + I(logTP^2) + AL + (1 + logTP || Dataset), REML = FALSE, data = dat)

AIC(m0, m_land, m_bio1, m_bio12, m_bio4, m_bio15, m_elev, m_trophic, m_ecosystem, m_lat)



mf0 <- lmer(FRic_z ~ logTP + I(logTP^2) + (1 + logTP || Dataset), REML = FALSE, data = dat)

mf_land <- lmer(FRic_z ~ logTP + I(logTP^2) + Landuse_class + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_bio1 <- lmer(FRic_z ~ logTP + I(logTP^2) + bio1_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_bio12 <- lmer(FRic_z ~ logTP + I(logTP^2) + bio12_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_bio4 <- lmer(FRic_z ~ logTP + I(logTP^2) + bio4_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_bio15 <- lmer(FRic_z ~ logTP + I(logTP^2) + bio15_1500 + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_elev <- lmer(FRic_z ~ logTP + I(logTP^2) + Elevation + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_trophic <- lmer(FRic_z ~ logTP + I(logTP^2) + Trophic_level + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_ecosystem <- lmer(FRic_z ~ logTP + I(logTP^2) + Ecosystem + (1 + logTP || Dataset), REML = FALSE, data = dat)
mf_lat <- lmer(FRic_z ~ logTP + I(logTP^2) + AL + (1 + logTP || Dataset), REML = FALSE, data = dat)

AIC(mf0, mf_land, mf_bio1, mf_bio12, mf_bio4, mf_bio15, mf_elev, mf_trophic, mf_ecosystem, mf_lat)


# Step 2
library(lme4)
library(MuMIn)
library(performance)

options(na.action = "na.fail")

global_rich <- lmer(
  Richness_z ~ logTP + I(logTP^2) +
    Landuse_class +  Ecosystem + 
    bio1_1500 + bio12_1500 +  bio15_1500 + Elevation +
    (1 + logTP || Dataset),
  REML = FALSE,
  data = dat
)


rich_dredge <- dredge(global_rich, rank = "AIC")

rich_dredge

model_performance(global_rich)

### Best final model TD
Richness_best <- lmer(
  Richness_z ~ logTP + I(logTP^2) +
    Landuse_class + 
    bio1_1500 + bio12_1500 + bio15_1500 + Elevation +
    (1 + logTP || Dataset),
  REML = FALSE,
  data = dat
)
check_collinearity(Richness_best)
summary(Richness_best)
model_performance(Richness_best)


options(na.action = "na.fail")

global_Frich <- lmer(
  FRic_z ~ logTP + I(logTP^2) +
    Landuse_class +  Ecosystem + 
    bio12_1500 + bio4_1500 +  bio15_1500 + Elevation +
    (1 + logTP || Dataset),
  REML = FALSE,
  data = dat
)

model_performance(global_Frich)
MuMIn::r.squaredGLMM(global_Frich)

Frich_dredge <- dredge(global_Frich, rank = "AIC")

Frich_dredge


rich_dredge_df  <- as.data.frame(rich_dredge)
Frich_dredge_df <- as.data.frame(Frich_dredge)


### Best final model FD
FRic_best <- lmer(
  FRic_z ~ logTP + I(logTP^2) +
    Landuse_class +
    bio12_1500 + bio4_1500 + bio15_1500 + Elevation +
    (1 + logTP || Dataset) ,
  REML = FALSE,
  data = dat
)
check_collinearity(FRic_best)
model_performance(FRic_best)
MuMIn::r.squaredGLMM(FRic_best)



library(broom.mixed)
library(openxlsx)

rich_fix  <- tidy(Richness_best, effects = "fixed")
rich_ran  <- tidy(Richness_best, effects = "ran_pars")
rich_fit  <- glance(Richness_best)

fric_fix  <- tidy(FRic_best, effects = "fixed")
fric_ran  <- tidy(FRic_best, effects = "ran_pars")
fric_fit  <- glance(FRic_best)

write.xlsx(
  list(
    Richness_fixed = rich_fix,
    Richness_random = rich_ran,
    Richness_fit = rich_fit,
    FRic_fixed = fric_fix,
    FRic_random = fric_ran,
    FRic_fit = fric_fit
  ),
  file = "Supplementary Data2.xlsx",
  overwrite = TRUE
)
