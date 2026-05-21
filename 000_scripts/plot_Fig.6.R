library(dplyr)
library(purrr)
library(tidyr)
library(lmerTest)
library(MuMIn)
library(openxlsx)
library(ggplot2)
library(ggeffects)
library(patchwork)


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

dat_env <- dat %>%
  mutate(
    bio1_1500_z  = as.numeric(scale(bio1_1500)),
    bio12_1500_z = as.numeric(scale(bio12_1500)),
    bio4_1500_z  = as.numeric(scale(bio4_1500)),
    bio15_1500_z = as.numeric(scale(bio15_1500)),
    Elevation_z  = as.numeric(scale(Elevation)),
    Dataset = as.factor(Dataset),
    Landuse_class = as.factor(Landuse_class)
  )

fit_tp_env_table <- function(data, predictor) {
  
  form_lin <- as.formula(
    paste0("logTP ~ ", predictor, " + (1 | Dataset)")
  )
  
  form_quad <- as.formula(
    paste0("logTP ~ ", predictor, " + I(", predictor, "^2) + (1 | Dataset)")
  )
  
  m_lin <- lmer(form_lin, data = data, REML = FALSE)
  m_quad <- lmer(form_quad, data = data, REML = FALSE)
  
  aic_tab <- AIC(m_lin, m_quad)
  aic_lin <- aic_tab$AIC[1]
  aic_quad <- aic_tab$AIC[2]
  
  best_model <- if ((aic_lin - aic_quad) > 2) "Quadratic" else "Linear"
  chosen_model <- if (best_model == "Quadratic") m_quad else m_lin
  
  sm <- summary(chosen_model)
  coef_tab <- as.data.frame(sm$coefficients)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  r2_tab <- MuMIn::r.squaredGLMM(chosen_model)
  
  coef_tab %>%
    transmute(
      response = "logTP",
      predictor = predictor,
      model_type = best_model,
      term = term,
      coefficient = round(Estimate, 4),
      std_error = round(`Std. Error`, 4),
      t_value = round(`t value`, 4),
      denominator_df = round(df, 4),
      p_value = round(`Pr(>|t|)`, 4),
      pseudo_R2_marginal = round(r2_tab[1], 4),
      pseudo_R2_conditional = round(r2_tab[2], 4),
      AIC_linear = round(aic_lin, 4),
      AIC_quadratic = round(aic_quad, 4)
    )
}

fit_tp_factor_table <- function(data, predictor) {
  
  form_fac <- as.formula(
    paste0("logTP ~ ", predictor, " + (1 | Dataset)")
  )
  
  m_fac <- lmer(form_fac, data = data, REML = FALSE)
  
  sm <- summary(m_fac)
  coef_tab <- as.data.frame(sm$coefficients)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  r2_tab <- MuMIn::r.squaredGLMM(m_fac)
  
  coef_tab %>%
    transmute(
      response = "logTP",
      predictor = predictor,
      model_type = "Factor",
      term = term,
      coefficient = round(Estimate, 4),
      std_error = round(`Std. Error`, 4),
      t_value = round(`t value`, 4),
      denominator_df = round(df, 4),
      p_value = round(`Pr(>|t|)`, 4),
      pseudo_R2_marginal = round(r2_tab[1], 4),
      pseudo_R2_conditional = round(r2_tab[2], 4),
      AIC_linear = NA_real_,
      AIC_quadratic = NA_real_
    )
}

tp_cont_vars <- c(
  "bio1_1500_z",
  "bio12_1500_z",
  "bio4_1500_z",
  "bio15_1500_z",
  "Elevation_z"
)

tp_factor_vars <- c("Landuse_class")

tp_table_cont <- map_dfr(
  tp_cont_vars,
  ~ fit_tp_env_table(dat_env, .x)
)

tp_table_fac <- map_dfr(
  tp_factor_vars,
  ~ fit_tp_factor_table(dat_env, .x)
)

tp_table2 <- bind_rows(tp_table_cont, tp_table_fac)

print(tp_table2)

tp_table2_noint <- tp_table2 %>%
  filter(term != "(Intercept)")

print(tp_table2_noint)

tp_sig_check <- tp_table2 %>%
  mutate(significant = ifelse(p_value < 0.05, "Yes", "No"))

print(tp_sig_check)

dat_env <- dat %>%
  mutate(
    bio12_1500_z = as.numeric(scale(bio12_1500)),
    bio15_1500_z = as.numeric(scale(bio15_1500)),
    Elevation_z  = as.numeric(scale(Elevation)),
    Dataset = as.factor(Dataset),
    Landuse_class = as.factor(Landuse_class)
  )

m_tp_bio12 <- lmer(
  logTP ~ bio12_1500_z + I(bio12_1500_z^2) + (1 | Dataset),
  data = dat_env,
  REML = FALSE
)

m_tp_bio15 <- lmer(
  logTP ~ bio15_1500_z + I(bio15_1500_z^2) + (1 | Dataset),
  data = dat_env,
  REML = FALSE
)

m_tp_elev <- lmer(
  logTP ~ Elevation_z + I(Elevation_z^2) + (1 | Dataset),
  data = dat_env,
  REML = FALSE
)

m_tp_landuse <- lmer(
  logTP ~ Landuse_class + (1 | Dataset),
  data = dat_env,
  REML = FALSE
)

r2_bio12   <- MuMIn::r.squaredGLMM(m_tp_bio12)[2]
r2_bio15   <- MuMIn::r.squaredGLMM(m_tp_bio15)[2]
r2_elev    <- MuMIn::r.squaredGLMM(m_tp_elev)[2]
r2_landuse <- MuMIn::r.squaredGLMM(m_tp_landuse)[2]

pred_bio12 <- as.data.frame(ggpredict(m_tp_bio12, terms = "bio12_1500_z [all]"))
pred_bio15 <- as.data.frame(ggpredict(m_tp_bio15, terms = "bio15_1500_z [all]"))
pred_elev  <- as.data.frame(ggpredict(m_tp_elev,  terms = "Elevation_z [all]"))

n_groups <- length(unique(dat_env$Dataset))
dataset_cols <- colorRampPalette(
  c("#2E5090", "#7FB3D5", "#ECECEC", "#E67E22", "#D73027")
)(n_groups)
names(dataset_cols) <- sort(unique(as.character(dat_env$Dataset)))

landuse_levels <- levels(dat_env$Landuse_class)

landuse_cols <- c(
  Forest = "#C6DBEF",
  Human   = "#E67E22",
  Mixed   = "#7FB3D5"
)

plot_tp_env <- function(data, xvar, pred_df, xlab, r2_value) {
  
  x_max <- max(data[[xvar]], na.rm = TRUE)
  y_max <- max(data$logTP, na.rm = TRUE)
  
  ggplot(data, aes(x = .data[[xvar]], y = logTP, color = Dataset)) +
    geom_point(alpha = 0.03, size = 0.4) +
    geom_ribbon(
      data = pred_df,
      aes(x = x, ymin = conf.low, ymax = conf.high),
      inherit.aes = FALSE,
      fill = "grey50", alpha = 0.2
    ) +
    geom_line(
      data = pred_df,
      aes(x = x, y = predicted),
      inherit.aes = FALSE,
      color = "black", linewidth = 1.3
    ) +
    annotate(
      "text",
      x = x_max,
      y = y_max,
      label = paste0("Pseudo R\u00B2 = ", sprintf("%.3f", r2_value)),
      hjust = 1.02, vjust = 1.2, size = 4.2
    ) +
    scale_color_manual(values = dataset_cols) +
    labs(x = xlab, y = "log(TP)") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
}

p_tp_landuse <- ggplot(dat_env, aes(x = Landuse_class, y = logTP, fill = Landuse_class)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "black",
    alpha = 0.3
  ) +
  geom_jitter(
    aes(color = Landuse_class),
    width = 0.15,
    alpha = 0.03,
    size = 0.5
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Pseudo R\u00B2 = ", sprintf("%.3f", r2_landuse)),
    hjust = 1.1, vjust = 1.4, size = 4.2
  ) +
  scale_fill_manual(values = landuse_cols) +
  scale_color_manual(values = landuse_cols) +
  labs(x = "Land-use", y = "log(TP)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

p_tp_bio12 <- plot_tp_env(
  data = dat_env,
  xvar = "bio12_1500_z",
  pred_df = pred_bio12,
  xlab = "Z-Mean Annual Precipitation",
  r2_value = r2_bio12
)

p_tp_bio15 <- plot_tp_env(
  data = dat_env,
  xvar = "bio15_1500_z",
  pred_df = pred_bio15,
  xlab = "Z-Precipitation Seasonality",
  r2_value = r2_bio15
)

p_tp_elev <- plot_tp_env(
  data = dat_env,
  xvar = "Elevation_z",
  pred_df = pred_elev,
  xlab = "Z-Elevation",
  r2_value = r2_elev
)

tp_4panel <- (p_tp_bio12 | p_tp_bio15) /
  (p_tp_elev  | p_tp_landuse)

tp_4panel <- tp_4panel + plot_annotation(tag_levels = "a")

print(tp_4panel)
