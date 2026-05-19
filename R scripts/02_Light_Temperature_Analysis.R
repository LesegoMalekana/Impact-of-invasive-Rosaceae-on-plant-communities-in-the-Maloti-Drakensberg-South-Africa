#This script is for temperature and light analyses.


# =========================
# Microclimate analyses
# Temperature and light, 7-day-thinned daily means
# =========================

# Step 0: Packages
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(broom)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(stringr)
library(openxlsx)
library(patchwork)

# Step 1: Styling
pal_trt <- c(Control = "#0072B2", Invaded = "#D55E00")

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 0.6),
    axis.line.y = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", size = 11, angle = 90),
    panel.spacing.y = unit(1.2, "lines")
  )

# Step 2: Load data
file_path <- "Your location of Cleaned_Light_Temperature_Standardised_Publish.csv"

lt_raw <- readr::read_csv(file_path, show_col_types = FALSE)

# Step 3: Clean labels and create date
lt <- lt_raw %>%
  dplyr::mutate(
    treatment = dplyr::if_else(treatment == "Uninvaded", "Control", treatment),
    treatment = factor(treatment, levels = c("Control", "Invaded")),
    site = factor(site, levels = c("Low", "Middle", "Top")),
    date = as.Date(datetime)
  )

# Step 4: Daily means
daily <- lt %>%
  dplyr::group_by(site, treatment, date) %>%
  dplyr::summarise(
    temp_mean = mean(temperature_c, na.rm = TRUE),
    light_mean = mean(light_lux, na.rm = TRUE),
    light_log1p_mean = mean(log1p(light_lux), na.rm = TRUE),
    .groups = "drop"
  )

cat("\n===== Daily counts before thinning =====\n")
daily %>%
  dplyr::count(site, treatment, name = "n_days") %>%
  dplyr::arrange(site, treatment) %>%
  print(n = Inf)

# Step 5: Thin to every 7th day
daily_7d <- daily %>%
  dplyr::arrange(site, treatment, date) %>%
  dplyr::group_by(site, treatment) %>%
  dplyr::mutate(day_index = dplyr::row_number()) %>%
  dplyr::filter(day_index %% 7 == 1) %>%
  dplyr::ungroup()

cat("\n===== Daily counts after 7-day thinning =====\n")
daily_7d %>%
  dplyr::count(site, treatment, name = "n_days_retained") %>%
  dplyr::arrange(site, treatment) %>%
  print(n = Inf)

# Step 6: Overall models
mod_temp <- stats::lm(temp_mean ~ treatment * site, data = daily_7d)
mod_light <- stats::lm(light_log1p_mean ~ treatment * site, data = daily_7d)

anova_temp <- car::Anova(mod_temp, type = 2)
anova_light <- car::Anova(mod_light, type = 2)

print(anova_temp)
print(anova_light)

# Step 7: Helper functions
p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

make_sig_label <- function(anova_obj) {
  tab <- as.data.frame(anova_obj)
  p_col <- names(tab)[grepl("Pr\\(>.*\\)", names(tab))]

  if (length(p_col) == 0) {
    stop("Could not find p-value column.")
  }

  tab$term <- rownames(tab)

  p_trt  <- tab[[p_col]][match("treatment", tab$term)]
  p_site <- tab[[p_col]][match("site", tab$term)]
  p_int  <- tab[[p_col]][match("treatment:site", tab$term)]

  paste0(
    "Site ", p_to_stars(p_site),
    ", Treatment ", p_to_stars(p_trt),
    ", Site × Treatment ", p_to_stars(p_int)
  )
}

get_cld_df <- function(model, metric_name) {
  em <- emmeans::emmeans(model, ~ site * treatment)

  multcomp::cld(
    em,
    Letters = letters,
    adjust = "bonferroni"
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      .group = stringr::str_trim(.group),
      Metric = metric_name,
      site = factor(site, levels = c("Low", "Middle", "Top")),
      treatment = factor(treatment, levels = c("Control", "Invaded"))
    )
}

anova_to_df <- function(anova_obj, metric_name) {
  tab <- as.data.frame(anova_obj)
  tab$Term <- rownames(tab)
  tab$Metric <- metric_name
  rownames(tab) <- NULL

  tab %>%
    dplyr::relocate(Metric, Term)
}

get_site_contrasts <- function(model, metric_name) {
  em <- emmeans::emmeans(model, ~ treatment | site)

  emmeans::contrast(em, method = "pairwise", adjust = "bonferroni") %>%
    as.data.frame() %>%
    dplyr::mutate(Metric = metric_name, .before = 1)
}

# Step 8: Compact letters
cld_temp <- get_cld_df(mod_temp, "Daily mean temperature (°C)")
cld_light <- get_cld_df(mod_light, "Daily mean light [log(1 + lux)]")

cld_all <- dplyr::bind_rows(cld_temp, cld_light)

# Step 9: Plotting data
plot_df <- daily_7d %>%
  dplyr::select(site, treatment, temp_mean, light_log1p_mean) %>%
  tidyr::pivot_longer(
    cols = c(temp_mean, light_log1p_mean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      temp_mean = "Daily mean temperature (°C)",
      light_log1p_mean = "Daily mean light [log(1 + lux)]"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Daily mean temperature (°C)",
        "Daily mean light [log(1 + lux)]"
      )
    )
  )

# Step 10: Label positions
y_pos <- plot_df %>%
  dplyr::group_by(Metric, site, treatment) %>%
  dplyr::summarise(y = max(Value, na.rm = TRUE), .groups = "drop")

metric_pad <- plot_df %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    rng = max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pad = pmax(0.10 * rng, 0.18)
  )

letters_df <- cld_all %>%
  dplyr::left_join(y_pos, by = c("Metric", "site", "treatment")) %>%
  dplyr::left_join(metric_pad, by = "Metric") %>%
  dplyr::mutate(y = y + pad)

ann_df <- tibble::tibble(
  Metric = c(
    "Daily mean temperature (°C)",
    "Daily mean light [log(1 + lux)]"
  ),
  label = c(
    make_sig_label(anova_temp),
    make_sig_label(anova_light)
  )
)

# Step 11: Main microclimate figure
p_micro <- ggplot(plot_df, aes(x = site, y = Value, fill = treatment)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.9,
    outlier.shape = 16,
    outlier.size = 1.8,
    colour = "black"
  ) +
  geom_text(
    data = letters_df,
    aes(x = site, y = y, group = treatment, label = .group),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  geom_text(
    data = ann_df,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.02,
    vjust = 1.2,
    size = 3.0
  ) +
  facet_grid(
    Metric ~ .,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(values = pal_trt) +
  labs(
    x = "Site",
    y = NULL,
    fill = "Treatment"
  ) +
  coord_cartesian(clip = "off") +
  theme_pub +
  theme(
    plot.margin = margin(12, 60, 12, 12)
  )

print(p_micro)

# Step 12: Welch t-test tables
make_inv_first <- function(x) {
  if (is.factor(x)) {
    stats::relevel(x, ref = "Invaded")
  } else {
    factor(x, levels = c("Invaded", "Control"))
  }
}

t_temp_7d <- daily_7d %>%
  dplyr::mutate(treatment = make_inv_first(treatment)) %>%
  dplyr::group_by(site) %>%
  dplyr::do({
    tt <- t.test(temp_mean ~ treatment, data = ., var.equal = FALSE)
    broom::tidy(tt)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    t_value = statistic,
    df = parameter,
    p_value = p.value,
    ci_low = conf.low,
    ci_high = conf.high,
    diff_C = estimate
  ) %>%
  dplyr::transmute(
    site,
    mean_invaded_C = estimate1,
    mean_control_C = estimate2,
    diff_C,
    ci_low_C = ci_low,
    ci_high_C = ci_high,
    t_value,
    df,
    p_value,
    pct_change = 100 * diff_C / mean_control_C
  )

t_light_7d <- daily_7d %>%
  dplyr::mutate(treatment = make_inv_first(treatment)) %>%
  dplyr::group_by(site) %>%
  dplyr::do({
    tt <- t.test(light_log1p_mean ~ treatment, data = ., var.equal = FALSE)
    broom::tidy(tt)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    t_value = statistic,
    df = parameter,
    p_value = p.value,
    ci_low = conf.low,
    ci_high = conf.high,
    diff_log = estimate
  ) %>%
  dplyr::mutate(
    ratio = exp(diff_log),
    ratio_lwr = exp(ci_low),
    ratio_upr = exp(ci_high),
    pct_change = 100 * (ratio - 1),
    pct_lwr = 100 * (ratio_lwr - 1),
    pct_upr = 100 * (ratio_upr - 1)
  ) %>%
  dplyr::transmute(
    site,
    diff_log,
    ci_low_log = ci_low,
    ci_high_log = ci_high,
    t_value,
    df,
    p_value,
    pct_change,
    pct_lwr,
    pct_upr
  )

print(t_temp_7d)
print(t_light_7d)

# Step 13: Export tables
anova_export <- dplyr::bind_rows(
  anova_to_df(anova_temp, "Daily mean temperature (°C)"),
  anova_to_df(anova_light, "Daily mean light [log(1 + lux)]")
)

contrasts_export <- dplyr::bind_rows(
  get_site_contrasts(mod_temp, "Daily mean temperature (°C)"),
  get_site_contrasts(mod_light, "Daily mean light [log(1 + lux)]")
)

cld_export <- cld_all %>%
  dplyr::select(Metric, site, treatment, .group)

out_dir <- dirname(file_path)
out_xlsx <- file.path(out_dir, "Microclimate_model_results_and_Welch_tests.xlsx")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Model_tests")
openxlsx::addWorksheet(wb, "Within_site_contrasts")
openxlsx::addWorksheet(wb, "Compact_letters")
openxlsx::addWorksheet(wb, "Welch_temperature")
openxlsx::addWorksheet(wb, "Welch_light")

openxlsx::writeData(wb, "Model_tests", anova_export)
openxlsx::writeData(wb, "Within_site_contrasts", contrasts_export)
openxlsx::writeData(wb, "Compact_letters", cld_export)
openxlsx::writeData(wb, "Welch_temperature", t_temp_7d)
openxlsx::writeData(wb, "Welch_light", t_light_7d)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Excel file saved to:\n", out_xlsx, "\n")

# Step 14: Export figure
ggsave(
  filename = file.path(out_dir, "Figure_Microclimate_Temperature_Light.png"),
  plot = p_micro,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Figure_Microclimate_Temperature_Light.pdf"),
  plot = p_micro,
  width = 9,
  height = 8,
  units = "in",
  bg = "white"
)

