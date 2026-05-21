
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


# LM or Quadratic 
# ---------- Richness_z ----------
mod_rich_lin <- lmer(
  Richness_z ~ logTP +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

mod_rich_quad <- lmer(
  Richness_z ~ logTP + I(logTP^2) +
    (1 + logTP  || Dataset),
  data = dat,
  REML = FALSE
)


anova(mod_rich_lin, mod_rich_quad)
AIC(mod_rich_lin, mod_rich_quad)
BIC(mod_rich_lin, mod_rich_quad)


# ---------- FRic ----------
mod_fric_lin <- lmer(
  FRic_z ~ logTP +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

mod_fric_quad <- lmer(
  FRic_z ~ logTP + I(logTP^2) +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

anova(mod_fric_lin, mod_fric_quad)
AIC(mod_fric_lin, mod_fric_quad)
BIC(mod_fric_lin, mod_fric_quad)


### Fig 3
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggeffects)
library(patchwork)

get_math_labels_list <- function(model) {
  coefs <- summary(model)$coefficients
  b1 <- coefs["logTP", "Estimate"]; s1 <- coefs["logTP", "Std. Error"]
  b2 <- coefs["I(logTP^2)", "Estimate"]; s2 <- coefs["I(logTP^2)", "Std. Error"]
  
  f_fmt <- function(val, se) {
    b_sci <- formatC(val, format = "e", digits = 1)
    b_p <- strsplit(b_sci, "e")[[1]]
    s_sci <- formatC(se, format = "e", digits = 1)
    s_p <- strsplit(s_sci, "e")[[1]]
    sprintf("%.1f %%*%% 10^%d ~ (SE == %.1f %%*%% 10^%d)", 
            as.numeric(b_p[1]), as.integer(b_p[2]),
            as.numeric(s_p[1]), as.integer(s_p[2]))
  }
  
  c(
    line1 = paste0("italic(beta)[1]==", f_fmt(b1, s1)),
    line2 = paste0("italic(beta)[2]==", f_fmt(b2, s2))
  )
}

blue_pal <- function(n) colorRampPalette(c("#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))(n)
orange_pal <- function(n) colorRampPalette(c("#FDBE85", "#FD8D3C", "#E6550D", "#A63603"))(n)

make_gradient_plot <- function(data, response, y_lab, model, 
                               line_palette_fun, 
                               text_col = "black") {
  
  data <- data %>% mutate(Dataset = as.factor(Dataset))
  tp_range <- data %>%
    group_by(Dataset) %>%
    summarise(min_t = min(logTP, na.rm = TRUE),
              max_t = max(logTP, na.rm = TRUE), .groups = "drop")
  
  pred_ind <- map_dfr(seq_len(nrow(tp_range)), function(i) {
    df_sub <- data.frame(
      Dataset = tp_range$Dataset[i],
      logTP = seq(tp_range$min_t[i], tp_range$max_t[i], length.out = 50)
    )
    df_sub$pred <- predict(model, newdata = df_sub, re.form = NULL)
    df_sub
  })
  
  pred_overall <- as.data.frame(ggpredict(model, terms = "logTP [all]"))
  dataset_levels <- levels(data$Dataset)
  color_values <- line_palette_fun(length(dataset_levels))
  names(color_values) <- dataset_levels
  x_rng <- range(data$logTP, na.rm = TRUE)
  label_list <- get_math_labels_list(model)
  x_pos <- x_rng[2] - 0.05 * diff(x_rng) 
  
  ggplot() +
    geom_point(data = data, aes(x = logTP, y = .data[[response]], color = Dataset),
               alpha = 0.02, size = 1) +
    geom_line(data = pred_ind, aes(x = logTP, y = pred, group = Dataset, color = Dataset),
              linewidth = 0.7, alpha = 0.5) +
    geom_ribbon(data = pred_overall, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "grey20", alpha = 0.12) +
    geom_line(data = pred_overall, aes(x = x, y = predicted),
              color = "black", linewidth = 1.4) +
    annotate("text", x = x_pos, y = 4.2, 
             label = label_list["line1"], parse = TRUE, 
             size = 5, hjust = 1, color = text_col) +
    annotate("text", x = x_pos, y = 3.3, 
             label = label_list["line2"], parse = TRUE, 
             size = 5, hjust = 1, color = text_col) +
    
    scale_color_manual(values = color_values) +
    coord_cartesian(ylim = c(-4, 4)) +
    scale_y_continuous(breaks = seq(-4, 4, by = 2)) +
    labs(x = "log(TP)", y = y_lab) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1, 1, 1, 1), "lines"))
}

p_rich <- make_gradient_plot(
  data = dat,
  response = "Richness_z", 
  y_lab = "Z-Richness",
  model = mod_rich_quad,
  line_palette_fun = blue_pal
)

p_fric <- make_gradient_plot(
  data = dat,
  response = "FRic_z", 
  y_lab = "Z-FRic",
  model = mod_fric_quad,
  line_palette_fun = orange_pal
)

final_plot <- p_rich | p_fric
final_plot <- final_plot + plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 24))

print(final_plot)




### Fig 4
get_math_labels_list <- function(model) {
  coefs <- summary(model)$coefficients
  b1 <- coefs["logTP", "Estimate"]; s1 <- coefs["logTP", "Std. Error"]
  b2 <- coefs["I(logTP^2)", "Estimate"]; s2 <- coefs["I(logTP^2)", "Std. Error"]
  
  f_fmt <- function(val, se) {
    b_sci <- formatC(val, format = "e", digits = 1)
    b_p <- strsplit(b_sci, "e")[[1]]
    s_sci <- formatC(se, format = "e", digits = 1)
    s_p <- strsplit(s_sci, "e")[[1]]
    
    sprintf("%.1f %%*%% 10^%d ~ (SE == %.1f %%*%% 10^%d)", 
            as.numeric(b_p[1]), as.integer(b_p[2]),
            as.numeric(s_p[1]), as.integer(s_p[2]))
  }
  
  c(
    line1 = paste0("italic(beta)[1]==", f_fmt(b1, s1)),
    line2 = paste0("italic(beta)[2]==", f_fmt(b2, s2))
  )
}

fit_quad_mixed <- function(data, response) {
  form1 <- as.formula(
    paste0(response, " ~ logTP + I(logTP^2) + (1 + logTP || Dataset)")
  )
  
  lmer(
    form1,
    data = data,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
}

make_combined_plot <- function(data, response, group_var, y_lab) {
  
  group_sym <- rlang::sym(group_var)
  groups <- sort(unique(as.character(data[[group_var]])))
  pal_list <- list(
    c("#084594", "#4292C6", "#C6DBEF"), 
    c("#8C2D04", "#EC7014", "#FDD0A2") 
  )
  
  dataset_color_vec <- c()
  group_main_color_vec <- c()
  
  for (i in seq_along(groups)) {
    g_name <- groups[i]
    g_datasets <- sort(unique(as.character(data$Dataset[data[[group_var]] == g_name])))
    g_cols <- colorRampPalette(pal_list[[i]])(length(g_datasets))
    names(g_cols) <- g_datasets
    dataset_color_vec <- c(dataset_color_vec, g_cols)
    group_main_color_vec[g_name] <- pal_list[[i]][1]
  }
  
  nested_res <- data %>%
    group_by(!!group_sym) %>%
    nest() %>%
    mutate(
      model = map(data, ~ fit_quad_mixed(.x, response)),
      
      pred_ind = map2(data, model, ~ {
        tp_range <- .x %>%
          group_by(Dataset) %>%
          summarise(
            min_t = min(logTP, na.rm = TRUE),
            max_t = max(logTP, na.rm = TRUE),
            .groups = "drop"
          )
        
        ind_list <- lapply(1:nrow(tp_range), function(j) {
          df_sub <- data.frame(
            Dataset = tp_range$Dataset[j],
            logTP = seq(tp_range$min_t[j], tp_range$max_t[j], length.out = 50)
          )
          df_sub$pred <- predict(.y, newdata = df_sub, re.form = NULL)
          df_sub
        })
        do.call(rbind, ind_list)
      }),
      
      pred_overall = map(model, ~ as.data.frame(ggpredict(.x, terms = "logTP [all]"))),
      label_list = map(model, get_math_labels_list)
    )
  
  df_points <- nested_res %>%
    select(!!group_sym, data) %>%
    unnest(data) %>%
    mutate(col = dataset_color_vec[as.character(Dataset)])
  
  df_ind <- nested_res %>%
    select(!!group_sym, pred_ind) %>%
    unnest(pred_ind) %>%
    mutate(col = dataset_color_vec[as.character(Dataset)])
  
  df_overall <- nested_res %>%
    select(!!group_sym, pred_overall) %>%
    unnest(pred_overall) %>%
    mutate(col = group_main_color_vec[as.character(!!group_sym)])
  df_labels <- nested_res %>%
    select(!!group_sym, label_list) %>%
    unnest_longer(label_list, indices_to = "line_id") %>%
    ungroup() %>%
    mutate(
      x = max(df_points$logTP, na.rm = TRUE),
      group_idx = as.numeric(factor(!!group_sym)),
      y = 3.8 - (group_idx - 1) * 1.4 - if_else(line_id == "line2", 0.5, 0),
      col = group_main_color_vec[as.character(!!group_sym)],
      final_label = if_else(
        line_id == "line1",
        paste0("bold('", !!group_sym, "'):~", label_list),
        label_list
      )
    )
  
  ggplot() +
    geom_point(
      data = df_points,
      aes(x = logTP, y = .data[[response]], color = col),
      alpha = 0.05, size = 0.3
    ) +
    geom_line(
      data = df_ind,
      aes(x = logTP, y = pred, color = col, group = Dataset),
      linewidth = 0.5, alpha = 0.4
    ) +
    geom_ribbon(
      data = df_overall,
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = col),
      alpha = 0.15
    ) +
    geom_line(
      data = df_overall,
      aes(x = x, y = predicted, color = col),
      linewidth = 1.2
    ) +
    geom_text(
      data = df_labels,
      aes(x = x, y = y, label = final_label, color = col),
      parse = TRUE, hjust = 1, size = 3.2
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    coord_cartesian(ylim = c(-4, 4)) +
    scale_y_continuous(breaks = seq(-4, 4, by = 2)) +
    labs(x = "log(TP)", y = y_lab) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none"
    )
}

p_trophic_rich <- make_combined_plot(dat, "Richness_z", "Trophic_level", "Z-Richness")
p_trophic_fric <- make_combined_plot(dat, "FRic_z", "Trophic_level", "Z-FRic")
p_ecosystem_rich <- make_combined_plot(dat, "Richness_z", "Ecosystem", "Z-Richness")
p_ecosystem_fric <- make_combined_plot(dat, "FRic_z", "Ecosystem", "Z-FRic")

final_plot <- (p_trophic_rich | p_trophic_fric) / 
  (p_ecosystem_rich | p_ecosystem_fric) + 
  plot_annotation(tag_levels = "a")

final_plot




### Fig S1
organism_order <- c(
  "Aquatic Plants",
  "Benthic Diatoms",
  "Phytoplankton",
  "Fish",
  "Benthic Macroinvertebrates",
  "Zooplankton"
)

organism_labels <- c(
  "Aquatic Plants" = "Aquatic Plants",
  "Benthic Diatoms" = "Diatom",
  "Phytoplankton" = "Phytoplankton",
  "Fish" = "Fish",
  "Benthic Macroinvertebrates" = "Macroinvertebrates",
  "Zooplankton" = "Zooplankton"
)

dat <- dat %>%
  mutate(
    Dataset = as.character(Dataset),
    Organism = factor(Organism, levels = organism_order)
  ) %>%
  filter(!is.na(Organism))

fit_quad_mixed <- function(data, response) {
  form1 <- as.formula(
    paste0(response, " ~ logTP + I(logTP^2) + (1 + logTP || Dataset)")
  )
  
  lmer(
    form1,
    data = data,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa", 
      optCtrl = list(maxfun = 2e5),
      check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
      check.conv.grad = .makeCC(action = "ignore", tol = 1e-3)
    )
  )
}

get_math_labels_two_lines <- function(model) {
  coefs <- summary(model)$coefficients
  b1 <- coefs["logTP", "Estimate"]; s1 <- coefs["logTP", "Std. Error"]
  b2 <- coefs["I(logTP^2)", "Estimate"]; s2 <- coefs["I(logTP^2)", "Std. Error"]
  
  f_fmt <- function(val, se, idx) {
    b_sci <- formatC(val, format = "e", digits = 1); b_p <- strsplit(b_sci, "e")[[1]]
    s_sci <- formatC(se, format = "e", digits = 1); s_p <- strsplit(s_sci, "e")[[1]]
    sprintf("italic(beta)[%s]==%.1f %%*%% 10^%d ~ (SE == %.1f %%*%% 10^%d)",
            idx, as.numeric(b_p[1]), as.integer(b_p[2]), as.numeric(s_p[1]), as.integer(s_p[2]))
  }
  paste0("atop(", f_fmt(b1, s1, "1"), ",", f_fmt(b2, s2, "2"), ")")
}

make_base_palette <- function(response) {
  if (response == "Richness_z") {
    c("#08306B", "#4292C6", "#C6DBEF")
  } else if (response == "FRic_z") {
    c("#8C2D04", "#EC7014", "#FDD0A2")
  } else {
    c("#444444", "#888888", "#DDDDDD")
  }
}

make_organism_plot <- function(data, response, y_lab, ncol = 3) {
  
  group_sym <- sym("Organism")
  base_pal <- make_base_palette(response)
  main_col <- base_pal[1]
  groups <- levels(droplevels(data$Organism))
  dataset_color_vec <- c()
  for (g in groups) {
    g_dat <- data %>% filter(Organism == g)
    g_datasets <- sort(unique(as.character(g_dat$Dataset)))
    n_ds <- length(g_datasets)
    pal <- colorRampPalette(base_pal)(max(n_ds, 3))[1:n_ds]
    names(pal) <- g_datasets
    dataset_color_vec <- c(dataset_color_vec, pal)
  }
  
  nested_res <- data %>%
    group_by(!!group_sym) %>%
    nest() %>%
    mutate(
      model = map(data, ~ fit_quad_mixed(.x, response)),
      
      pred_ind = map2(data, model, ~{
        tp_range <- .x %>%
          group_by(Dataset) %>%
          summarise(min_t = min(logTP, na.rm = TRUE),
                    max_t = max(logTP, na.rm = TRUE), .groups = "drop")
        
        ind_list <- lapply(seq_len(nrow(tp_range)), function(j) {
          df_sub <- data.frame(
            Dataset = tp_range$Dataset[j],
            logTP = seq(tp_range$min_t[j], tp_range$max_t[j], length.out = 80)
          )
          tryCatch({
            df_sub$pred <- predict(.y, newdata = df_sub, re.form = NULL, allow.new.levels = TRUE)
            df_sub
          }, error = function(e) return(NULL))
        })
        do.call(rbind, Filter(Negate(is.null), ind_list))
      }),
      
      pred_overall = map(model, ~ as.data.frame(ggpredict(.x, terms = "logTP [all]"))),
      label_str = map_chr(model, get_math_labels_two_lines)
    )
  
  df_points <- nested_res %>% select(!!group_sym, data) %>% unnest(data) %>%
    mutate(col = dataset_color_vec[as.character(Dataset)])
  
  df_ind <- nested_res %>% select(!!group_sym, pred_ind) %>% unnest(pred_ind) %>%
    mutate(col = dataset_color_vec[as.character(Dataset)])
  
  df_overall <- nested_res %>% select(!!group_sym, pred_overall) %>% unnest(pred_overall) %>%
    mutate(col = main_col)
  
  df_labels <- nested_res %>% ungroup() %>%
    mutate(x = map_dbl(data, ~max(.x$logTP, na.rm = TRUE)),
           y = 3.6, col = main_col)
  
  ggplot() +
    geom_point(data = df_points, aes(x = logTP, y = .data[[response]], color = col),
               alpha = 0.2, size = 0.5) +
    geom_line(data = df_ind, aes(x = logTP, y = pred, color = col, group = Dataset),
              linewidth = 0.6, alpha = 0.6) +
    geom_ribbon(data = df_overall, aes(x = x, ymin = conf.low, ymax = conf.high, fill = col),
                alpha = 0.15) +
    geom_line(data = df_overall, aes(x = x, y = predicted, color = col),
              linewidth = 1.3, alpha = 0.8) +
    geom_text(data = df_labels, aes(x = x, y = y, label = label_str, color = col),
              parse = TRUE, hjust = 1.05, size = 2.8) +
    scale_color_identity() +
    scale_fill_identity() +
    facet_wrap(~ Organism, ncol = ncol, labeller = as_labeller(organism_labels)) +
    coord_cartesian(ylim = c(-4, 4.5)) +
    scale_y_continuous(breaks = seq(-4, 4, by = 2)) +
    labs(x = "log(TP)", y = y_lab) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

p_organism_rich <- make_organism_plot(dat, "Richness_z", "Z-Richness")
p_organism_fric <- make_organism_plot(dat, "FRic_z", "Z-FRic")

p_organism_combined <- p_organism_rich / p_organism_fric + plot_annotation(tag_levels = "a")

print(p_organism_combined)

