# Step 0: Load required packages
library(readr)
library(dplyr)

# Step 1: Define file path
file_path <- "Your file location"

# Step 2: Read the CSV file (no changes yet)
lt_data <- read_csv(file_path)

# Step 3: Explore dataset structure and summary
# --- Basic structure ---
cat("===== STRUCTURE =====\n")
str(lt_data)

# --- Column names ---
cat("\n===== COLUMN NAMES =====\n")
names(lt_data)

# --- First few rows ---
cat("\n===== HEAD (First 10 rows) =====\n")
print(head(lt_data, 10))

# --- Summary statistics ---
cat("\n===== SUMMARY =====\n")
summary(lt_data)

# --- Check data types and missing values ---
cat("\n===== MISSING VALUES PER COLUMN =====\n")
sapply(lt_data, function(x) sum(is.na(x)))

# --- Quick glimpse (tibble view) ---
cat("\n===== GLIMPSE =====\n")
glimpse(lt_data)



# Step 4: Deeper exploration (no mutations, just prints)

library(dplyr)
library(lubridate)

cat("\n===== DISTINCT LEVELS =====\n")
cat("Sites: "); print(sort(unique(lt_data$site)))
cat("Treatments: "); print(sort(unique(lt_data$treatment)))

cat("\n===== ROW COUNTS BY SITE & TREATMENT =====\n")
print(
  lt_data %>%
    count(site, treatment, name = "n") %>%
    arrange(site, treatment)
)

cat("\n===== DATE RANGE (OVERALL) =====\n")
overall_range <- range(lt_data$datetime)
print(overall_range)

cat("\n===== DATE RANGE BY SITE & TREATMENT =====\n")
print(
  lt_data %>%
    group_by(site, treatment) %>%
    summarise(start = min(datetime), end = max(datetime), .groups = "drop") %>%
    arrange(site, treatment)
)

cat("\n===== SAMPLING INTERVAL CHECK (mins) =====\n")
# Compute diffs (mins) per site/treatment to inspect typical logging interval
print(
  lt_data %>%
    arrange(site, treatment, datetime) %>%
    group_by(site, treatment) %>%
    mutate(dt_min = as.numeric(difftime(datetime, lag(datetime), units = "mins"))) %>%
    summarise(
      n_gaps = sum(!is.na(dt_min)),
      min_gap = min(dt_min, na.rm = TRUE),
      q25_gap = quantile(dt_min, 0.25, na.rm = TRUE),
      median_gap = median(dt_min, na.rm = TRUE),
      q75_gap = quantile(dt_min, 0.75, na.rm = TRUE),
      max_gap = max(dt_min, na.rm = TRUE),
      .groups = "drop"
    )
)

cat("\n===== BASIC RANGE BY SITE & TREATMENT =====\n")
print(
  lt_data %>%
    group_by(site, treatment) %>%
    summarise(
      n = dplyr::n(),
      temp_min = min(temperature_c, na.rm = TRUE),
      temp_q25 = quantile(temperature_c, 0.25, na.rm = TRUE),
      temp_med = median(temperature_c, na.rm = TRUE),
      temp_q75 = quantile(temperature_c, 0.75, na.rm = TRUE),
      temp_max = max(temperature_c, na.rm = TRUE),
      lux_min = min(light_lux, na.rm = TRUE),
      lux_q25 = quantile(light_lux, 0.25, na.rm = TRUE),
      lux_med = median(light_lux, na.rm = TRUE),
      lux_q75 = quantile(light_lux, 0.75, na.rm = TRUE),
      lux_max = max(light_lux, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(site, treatment)
)

cat("\n===== DUPLICATE TIMESTAMPS WITHIN SITE & TREATMENT? =====\n")
print(
  lt_data %>%
    count(site, treatment, datetime) %>%
    filter(n > 1) %>%
    arrange(site, treatment, datetime) %>%
    head(20)








  # Step 0: Libraries & theme ----------------------------------------------------
  library(dplyr)
  library(lubridate)
  library(readr)
  library(ggplot2)
  library(broom)
  library(tidyr)

  # Styling per your preferences
  pal_trt <- c(Control = "#7CAE00", Invaded = "#F8766D")
  theme_pub <- theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(),
      axis.ticks = element_line(),
      legend.position = "right"
    )

  # Step 1: Load -----------------------------------------------------------------
  file_path <- "Your file location"
  lt_raw <- read_csv(file_path, show_col_types = FALSE)

  # Step 2: Standardise treatment label (Uninvaded -> Control), add date --------
  lt <- lt_raw %>%
    mutate(
      treatment = recode(treatment, "Uninvaded" = "Control"),
      treatment = factor(treatment, levels = c("Control", "Invaded")),
      date = as.Date(datetime)
    )

  # Step 3: Daily means per site × treatment × date ------------------------------
  # (Light kept in original and log1p scales; we’ll use log1p for inference/plots)
  daily <- lt %>%
    group_by(site, treatment, date) %>%
    summarise(
      temp_mean = mean(temperature_c, na.rm = TRUE),
      light_mean = mean(light_lux, na.rm = TRUE),
      light_log1p_mean = mean(log1p(light_lux), na.rm = TRUE),
      .groups = "drop"
    )

  # Quick counts (how many days per cell)
  days_count <- daily %>% count(site, treatment, name = "n_days")
  print(days_count)

  # Step 4: Plots — daily mean Temperature (°C) ----------------------------------
  p_temp_box <- ggplot(daily, aes(x = treatment, y = temp_mean, fill = treatment)) +
    geom_boxplot(outlier.shape = 16, alpha = 0.8, width = 0.7) +
    facet_wrap(~ site, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = pal_trt) +
    labs(x = NULL, y = "Daily mean temperature (°C)", fill = "Treatment",
         title = "Daily mean temperature by treatment (boxplots)") +
    theme_pub
  print(p_temp_box)

  # Alternative: mean ± 95% CI across days
  p_temp_ci <- ggplot(daily, aes(x = treatment, y = temp_mean, colour = treatment)) +
    stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(width = 0.4)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2,
                 position = position_dodge(width = 0.4)) +
    facet_wrap(~ site, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = pal_trt) +
    labs(x = NULL, y = "Daily mean temperature (°C)", colour = "Treatment",
         title = "Daily mean temperature (mean ± 95% CI across days)") +
    theme_pub
  print(p_temp_ci)

  # Step 5: Plots — daily mean Light on log1p scale ------------------------------
  # Boxplots on log1p scale
  p_light_box <- ggplot(daily, aes(x = treatment, y = light_log1p_mean, fill = treatment)) +
    geom_boxplot(outlier.shape = 16, alpha = 0.8, width = 0.7) +
    facet_wrap(~ site, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = pal_trt) +
    labs(x = NULL, y = "Daily mean log(1 + light [lux])", fill = "Treatment",
         title = "Daily mean light (log1p) by treatment (boxplots)") +
    theme_pub
  print(p_light_box)

  # Mean ± 95% CI on log1p scale
  p_light_ci <- ggplot(daily, aes(x = treatment, y = light_log1p_mean, colour = treatment)) +
    stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(width = 0.4)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2,
                 position = position_dodge(width = 0.4)) +
    facet_wrap(~ site, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = pal_trt) +
    labs(x = NULL, y = "Daily mean log(1 + light [lux])", colour = "Treatment",
         title = "Daily mean light (log1p) (mean ± 95% CI across days)") +
    theme_pub
  print(p_light_ci)

  # Step 6: Inference — Welch t-tests by site ------------------------------------
  # Temperature (°C) on daily means
  t_temp <- daily %>%
    group_by(site) %>%
    do({
      tt <- t.test(temp_mean ~ treatment, data = ., var.equal = FALSE)
      broom::tidy(tt)
    }) %>% ungroup() %>%
    mutate(comparison = "Invaded - Control (°C)")
  print(t_temp)

  # Light on log1p scale (differences on log scale)
  t_light_log <- daily %>%
    group_by(site) %>%
    do({
      tt <- t.test(light_log1p_mean ~ treatment, data = ., var.equal = FALSE)
      broom::tidy(tt)
    }) %>% ungroup() %>%
    mutate(comparison = "Invaded - Control on log(1+lux)")
  print(t_light_log)

  # Step 7: Effect sizes ----------------------------------------------------------
  # 7a) Temperature: % change relative to Control on the original scale
  eff_temp <- daily %>%
    group_by(site, treatment) %>%
    summarise(mean_temp = mean(temp_mean), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = treatment, values_from = mean_temp) %>%
    mutate(temp_pct_change = 100 * (Invaded - Control) / Control) %>%
    select(site, Control_temp = Control, Invaded_temp = Invaded, temp_pct_change)

  print(eff_temp)

  # 7b) Light: % change as a **ratio** from log1p differences, back-transformed
  # From t-tests we have estimate and CI on the log1p scale. Convert to ratios.
  eff_light <- t_light_log %>%
    transmute(
      site,
      log_diff = estimate,                 # (Invaded - Control) on log(1+lux)
      log_diff_lwr = conf.low,
      log_diff_upr = conf.high,
      ratio = exp(log_diff) - 1,          # fractional change relative to Control
      ratio_lwr = exp(log_diff_lwr) - 1,
      ratio_upr = exp(log_diff_upr) - 1,
      light_pct_change = 100 * ratio,
      light_pct_change_lwr = 100 * ratio_lwr,
      light_pct_change_upr = 100 * ratio_upr
    )

  print(eff_light)

  # Optional: also report raw-scale % change using daily means on original scale
  eff_light_raw <- daily %>%
    group_by(site, treatment) %>%
    summarise(mean_light = mean(light_mean), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = mean_light) %>%
    mutate(light_pct_change_raw = 100 * (Invaded - Control) / Control) %>%
    select(site, Control_light = Control, Invaded_light = Invaded, light_pct_change_raw)

  print(eff_light_raw)











  # ---- Step 0: Libraries & styling --------------------------------------------
  library(dplyr)
  library(lubridate)
  library(readr)
  library(ggplot2)
  library(broom)
  library(tidyr)
  library(stringr)

  # Colours (your scheme)
  pal_trt <- c(Control = "#7CAE00", Invaded = "#F8766D")

  # Theme: ONLY left (y) and bottom (x) axis lines visible
  theme_pub <- theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.6),
      axis.line.y = element_line(colour = "black", linewidth = 0.6),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  # ---- Step 1: Load ------------------------------------------------------------
  file_path <- "C:/Users/lmalekana/OneDrive/UFS - Masters/Light and temperature/Cleaned_Light_Temperature_Standardised.csv"
  lt_raw <- read_csv(file_path, show_col_types = FALSE)

  # ---- Step 2: Relabel & add date ---------------------------------------------
  lt <- lt_raw %>%
    mutate(
      treatment = if_else(treatment == "Uninvaded", "Control", treatment),
      treatment = factor(treatment, levels = c("Control", "Invaded")),
      date = as.Date(datetime)
    )

  # ---- Step 3: Daily means (per site × treatment × date) ----------------------
  daily <- lt %>%
    group_by(site, treatment, date) %>%
    summarise(
      temp_mean = mean(temperature_c, na.rm = TRUE),
      light_mean = mean(light_lux, na.rm = TRUE),
      light_log1p_mean = mean(log1p(light_lux), na.rm = TRUE),
      .groups = "drop"
    )

  # (Optional) sanity check
  daily %>% count(site, treatment, name = "n_days") %>% arrange(site, treatment) %>% print()

  # ---- Step 4: Welch t-tests (daily means; by site) ---------------------------
  # Temperature (°C)
  t_temp <- daily %>%
    group_by(site) %>%
    do(broom::tidy(t.test(temp_mean ~ treatment, data = ., var.equal = FALSE))) %>%
    ungroup() %>%
    mutate(measure = "Temperature (°C)", comparison = "Invaded - Control")
  print(t_temp)

  # Light on log1p scale
  t_light_log <- daily %>%
    group_by(site) %>%
    do(broom::tidy(t.test(light_log1p_mean ~ treatment, data = ., var.equal = FALSE))) %>%
    ungroup() %>%
    mutate(measure = "log(1+light)", comparison = "Invaded - Control")
  print(t_light_log)

  # ---- Step 5: Effect sizes ----------------------------------------------------
  # 5a) Temperature: % change from Control
  eff_temp <- daily %>%
    group_by(site, treatment) %>%
    summarise(mean_temp = mean(temp_mean), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = mean_temp) %>%
    mutate(temp_pct_change = 100 * (Invaded - Control) / Control) %>%
    select(site, Control_temp = Control, Invaded_temp = Invaded, temp_pct_change)
  print(eff_temp)

  # 5b) Light: % change from Control using back-transformed log1p difference
  eff_light <- daily %>%
    group_by(site, treatment) %>%
    summarise(mean_log = mean(light_log1p_mean), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = mean_log) %>%
    transmute(
      site,
      log_diff = Invaded - Control,
      ratio = exp(log_diff) - 1,
      light_pct_change = 100 * ratio
    )
  print(eff_light)

  # (Optionally also show raw-scale % change for intuition; not used for inference)
  eff_light_raw <- daily %>%
    group_by(site, treatment) %>%
    summarise(mean_light = mean(light_mean), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = mean_light) %>%
    mutate(light_pct_change_raw = 100 * (Invaded - Control) / Control)
  print(eff_light_raw)

  # ---- Step 6: Plot summaries with mean ± 95% CI ------------------------------
  # Build summary table with CIs for plotting
  summ <- daily %>%
    group_by(site, treatment) %>%
    summarise(
      n = dplyr::n(),
      temp_mean = mean(temp_mean),
      temp_sd   = sd(temp_mean),
      temp_se   = temp_sd / sqrt(n),
      temp_ci   = qt(0.975, df = n - 1) * temp_se,
      light_log_mean = mean(light_log1p_mean),
      light_log_sd   = sd(light_log1p_mean),
      light_log_se   = light_log_sd / sqrt(n),
      light_log_ci   = qt(0.975, df = n - 1) * light_log_se,
      # Back-transform light for plotting on original lux scale
      light_bt_mean  = exp(light_log_mean) - 1,
      light_bt_lwr   = exp(light_log_mean - light_log_ci) - 1,
      light_bt_upr   = exp(light_log_mean + light_log_ci) - 1,
      .groups = "drop"
    )

  # A) Temperature: mean ± 95% CI
  p_temp_ci <- ggplot(
    summ,
    aes(x = site, y = temp_mean, colour = treatment)
  ) +
    geom_point(position = position_dodge(width = 0.45), size = 2.8) +
    geom_errorbar(
      aes(ymin = temp_mean - temp_ci, ymax = temp_mean + temp_ci),
      width = 0.15,
      position = position_dodge(width = 0.45)
    ) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(
      title = "A. Daily Mean Temperature (mean ± 95% CI)",
      x = "Site",
      y = "Daily Mean Temperature (°C)"
    ) +
    theme_pub
  print(p_temp_ci)

  # B) Light (back-transformed): mean ± 95% CI
  p_light_ci <- ggplot(
    summ,
    aes(x = site, y = light_bt_mean, colour = treatment)
  ) +
    geom_point(position = position_dodge(width = 0.45), size = 2.8) +
    geom_errorbar(
      aes(ymin = light_bt_lwr, ymax = light_bt_upr),
      width = 0.15,
      position = position_dodge(width = 0.45)
    ) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(
      title = "B. Daily Mean Light (back-transformed from log(1+lux))",
      x = "Site",
      y = "Daily Mean Light (Lux, daylight)"
    ) +
    theme_pub
  print(p_light_ci)

  # ---- (Optional) Save figures --------------------------------------------------
  # ggsave("daily_temp_ci.png", p_temp_ci, width = 8, height = 3.2, dpi = 300)
  # ggsave("daily_light_ci.png", p_light_ci, width = 8, height = 3.2, dpi = 300)



  # ---- Recompute summary with 95% CI (same logic) ------------------------------
  summ <- daily %>%
    group_by(site, treatment) %>%
    summarise(
      n = dplyr::n(),
      temp_mean = mean(temp_mean),
      temp_sd   = sd(temp_mean),
      temp_se   = temp_sd / sqrt(n),
      temp_ci   = qt(0.975, df = n - 1) * temp_se,   # 95% CI
      light_log_mean = mean(light_log1p_mean),
      light_log_sd   = sd(light_log1p_mean),
      light_log_se   = light_log_sd / sqrt(n),
      light_log_ci   = qt(0.975, df = n - 1) * light_log_se,
      light_bt_mean  = exp(light_log_mean) - 1,
      light_bt_lwr   = exp(light_log_mean - light_log_ci) - 1,
      light_bt_upr   = exp(light_log_mean + light_log_ci) - 1,
      .groups = "drop"
    )

  # Common dodge + little y padding so caps aren’t cut off
  pd <- position_dodge(width = 0.5)

  # ---- A) Temperature: mean ± 95% CI (beefier error bars) ----------------------
  p_temp_ci <- ggplot(summ, aes(x = site, y = temp_mean, colour = treatment)) +
    geom_errorbar(
      aes(ymin = temp_mean - temp_ci, ymax = temp_mean + temp_ci),
      width = 0.18,                      # cap width
      linewidth = 0.9,                   # thicker so it shows
      position = pd
    ) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(
      title = "A. Daily Mean Temperature (mean ± 95% CI)",
      x = "Site",
      y = "Daily Mean Temperature (°C)"
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub

  print(p_temp_ci)

  # ---- B) Light (back-transformed): mean ± 95% CI (same visual tweaks) ---------
  p_light_ci <- ggplot(summ, aes(x = site, y = light_bt_mean, colour = treatment)) +
    geom_errorbar(
      aes(ymin = light_bt_lwr, ymax = light_bt_upr),
      width = 0.18,
      linewidth = 0.9,
      position = pd
    ) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(
      title = "B. Daily Mean Light (back-transformed from log(1+lux))",
      x = "Site",
      y = "Daily Mean Light (Lux, daylight)"
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub

  print(p_light_ci)


  # === Inspect summary values used in plots =====================================
  summ_print <- summ %>%
    select(
      site, treatment,
      n,
      temp_mean,
      temp_ci,
      light_bt_mean,
      light_bt_lwr,
      light_bt_upr
    ) %>%
    arrange(site, treatment) %>%
    mutate(
      temp_lower = temp_mean - temp_ci,
      temp_upper = temp_mean + temp_ci
    ) %>%
    select(
      site, treatment, n,
      temp_mean, temp_lower, temp_upper,
      light_bt_mean, light_bt_lwr, light_bt_upr
    )

  print(summ_print, n = Inf)







  # ---- FIXED summary with 95% CI ----------------------------------------------
  summ <- daily %>%
    group_by(site, treatment) %>%
    summarise(
      n = dplyr::n(),
      temp_mean_val = mean(temp_mean),
      temp_sd_val   = sd(temp_mean),
      temp_se       = temp_sd_val / sqrt(n),
      temp_ci       = qt(0.975, df = n - 1) * temp_se,

      light_log_mean = mean(light_log1p_mean),
      light_log_sd   = sd(light_log1p_mean),
      light_log_se   = light_log_sd / sqrt(n),
      light_log_ci   = qt(0.975, df = n - 1) * light_log_se,

      # Back-transform light for plotting on original lux scale
      light_bt_mean  = exp(light_log_mean) - 1,
      light_bt_lwr   = exp(light_log_mean - light_log_ci) - 1,
      light_bt_upr   = exp(light_log_mean + light_log_ci) - 1,
      .groups = "drop"
    )

  # ---- Print the values (now with non-NA temp CIs) -----------------------------
  summ_print <- summ %>%
    transmute(
      site, treatment, n,
      temp_mean = temp_mean_val,
      temp_lower = temp_mean_val - temp_ci,
      temp_upper = temp_mean_val + temp_ci,
      light_bt_mean, light_bt_lwr, light_bt_upr
    ) %>%
    arrange(site, treatment)

  print(summ_print, n = Inf)

  # ---- Plots (mean ± 95% CI), using the corrected names ------------------------
  pd <- position_dodge(width = 0.5)

  p_temp_ci <- ggplot(summ, aes(x = site, y = temp_mean_val, colour = treatment)) +
    geom_errorbar(
      aes(ymin = temp_mean_val - temp_ci, ymax = temp_mean_val + temp_ci),
      width = 0.18, linewidth = 0.9, position = pd
    ) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(title = "A. Daily Mean Temperature (mean ± 95% CI)",
         x = "Site", y = "Daily Mean Temperature (°C)") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub
  print(p_temp_ci)

  p_light_ci <- ggplot(summ, aes(x = site, y = light_bt_mean, colour = treatment)) +
    geom_errorbar(
      aes(ymin = light_bt_lwr, ymax = light_bt_upr),
      width = 0.18, linewidth = 0.9, position = pd
    ) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(title = "B. Daily Mean Light (back-transformed from log(1+lux))",
         x = "Site", y = "Daily Mean Light (Lux, daylight)") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub
  print(p_light_ci)















  # Uses your existing 'daily' object
  summ <- daily %>%
    group_by(site, treatment) %>%
    summarise(
      n = dplyr::n(),
      # Temperature
      temp_mean_val = mean(temp_mean),
      temp_sd_val   = sd(temp_mean),
      temp_se       = temp_sd_val / sqrt(n),
      temp_ci       = qt(0.975, df = n - 1) * temp_se,
      # Light on log(1+lux)
      light_log_mean = mean(light_log1p_mean),
      light_log_sd   = sd(light_log1p_mean),
      light_log_se   = light_log_sd / sqrt(n),
      light_log_ci   = qt(0.975, df = n - 1) * light_log_se,
      .groups = "drop"
    )
  summ_print <- summ %>%
    transmute(
      site, treatment, n,
      temp_mean = temp_mean_val,
      temp_lower = temp_mean_val - temp_ci,
      temp_upper = temp_mean_val + temp_ci,
      light_log_mean,
      light_log_lower = light_log_mean - light_log_ci,
      light_log_upper = light_log_mean + light_log_ci
    ) %>%
    arrange(site, treatment)

  print(summ_print, n = Inf)
  pd <- position_dodge(width = 0.5)

  # A) Temperature (mean ± 95% CI)
  p_temp_ci <- ggplot(summ, aes(x = site, y = temp_mean_val, colour = treatment)) +
    geom_errorbar(aes(ymin = temp_mean_val - temp_ci, ymax = temp_mean_val + temp_ci),
                  width = 0.18, linewidth = 0.9, position = pd) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(title = "A. Daily Mean Temperature (mean ± 95% CI)",
         x = "Site", y = "Daily Mean Temperature (°C)") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub
  print(p_temp_ci)

  # B) Light on log(1+lux) (mean ± 95% CI)
  p_light_ci <- ggplot(summ, aes(x = site, y = light_log_mean, colour = treatment)) +
    geom_errorbar(aes(ymin = light_log_mean - light_log_ci, ymax = light_log_mean + light_log_ci),
                  width = 0.18, linewidth = 0.9, position = pd) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(title = "B. Daily Mean Light (log(1 + lux); mean ± 95% CI)",
         x = "Site", y = "Daily Mean Light (log(1 + lux))") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub
  print(p_light_ci)



  # ---- Combined stacked plot: Temperature (A) + Light (B) ----------------------
  library(patchwork)

  pd <- position_dodge(width = 0.5)

  # A) Temperature
  p_temp <- ggplot(summ, aes(x = site, y = temp_mean_val, colour = treatment)) +
    geom_errorbar(aes(ymin = temp_mean_val - temp_ci, ymax = temp_mean_val + temp_ci),
                  width = 0.18, linewidth = 0.9, position = pd) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(x = NULL, y = "Daily Mean Temperature (°C)") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub +
    theme(legend.position = "none") +         # remove individual legend
    annotate("text", x = -Inf, y = Inf, label = "A", fontface = "bold",
             hjust = -0.4, vjust = 1.6, size = 5)

  # B) Light (log scale)
  p_light <- ggplot(summ, aes(x = site, y = light_log_mean, colour = treatment)) +
    geom_errorbar(aes(ymin = light_log_mean - light_log_ci, ymax = light_log_mean + light_log_ci),
                  width = 0.18, linewidth = 0.9, position = pd) +
    geom_point(position = pd, size = 3.2) +
    scale_colour_manual(values = pal_trt, name = "Treatment") +
    labs(x = "Site", y = "Daily Mean Light (log(1 + lux))") +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    theme_pub +
    annotate("text", x = -Inf, y = Inf, label = "B", fontface = "bold",
             hjust = -0.4, vjust = 1.6, size = 5)

  # Combine with shared legend
  combined_plot <- p_temp + p_light +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 10))

  print(combined_plot)

  # Optional save
  # ggsave("Temp_Light_Stacked.png", combined_plot, width = 6, height = 7, dpi = 300)







  # ===== T-tests per site with Invaded − Control direction ======================
  library(dplyr)
  library(broom)
  library(tidyr)

  # Helper: ensure Invaded is the first level so t.test returns (Invaded - Control)
  make_inv_first <- function(x) {
    if (is.factor(x)) {
      stats::relevel(x, ref = "Invaded")
    } else {
      factor(x, levels = c("Invaded", "Control"))
    }
  }

  ## 1) Temperature (°C) ----------------------------------------------------------
  t_temp_tbl <- daily %>%
    mutate(treatment = make_inv_first(treatment)) %>%
    group_by(site) %>%
    do({
      tt <- t.test(temp_mean ~ treatment, data = ., var.equal = FALSE)
      broom::tidy(tt)
    }) %>%
    ungroup() %>%
    # Rename & add % change vs Control
    rename(
      t_value   = statistic,
      df        = parameter,
      p_value   = p.value,
      ci_low    = conf.low,
      ci_high   = conf.high,
      diff_C    = estimate   # this is Invaded - Control (°C)
    ) %>%
    # broom::tidy also gives estimate1 = mean in level 1 (Invaded), estimate2 = Control
    transmute(
      site,
      mean_invaded_C = estimate1,
      mean_control_C = estimate2,
      diff_C,
      ci_low_C = ci_low,
      ci_high_C = ci_high,
      t_value, df, p_value,
      pct_change = 100 * diff_C / mean_control_C  # (Invaded - Control)/Control * 100
    )

  cat("\n===== Welch t-tests: Temperature (°C) — Invaded minus Control =====\n")
  print(t_temp_tbl, n = Inf, digits = 4)

  ## 2) Light on log(1+lux) ------------------------------------------------------
  t_light_tbl <- daily %>%
    mutate(treatment = make_inv_first(treatment)) %>%
    group_by(site) %>%
    do({
      tt <- t.test(light_log1p_mean ~ treatment, data = ., var.equal = FALSE)
      broom::tidy(tt)
    }) %>%
    ungroup() %>%
    rename(
      t_value   = statistic,
      df        = parameter,
      p_value   = p.value,
      ci_low    = conf.low,
      ci_high   = conf.high,
      diff_log  = estimate   # this is Invaded - Control on log(1+lux)
    ) %>%
    # Back-transform to ratios and % change (with CI)
    mutate(
      ratio       = exp(diff_log),            # Invaded / Control on (1+lux) scale
      ratio_lwr   = exp(ci_low),
      ratio_upr   = exp(ci_high),
      pct_change  = 100 * (ratio - 1),
      pct_lwr     = 100 * (ratio_lwr - 1),
      pct_upr     = 100 * (ratio_upr - 1)
    ) %>%
    transmute(
      site,
      diff_log, ci_low_log = ci_low, ci_high_log = ci_high,
      t_value, df, p_value,
      pct_change, pct_lwr, pct_upr
    )

  cat("\n===== Welch t-tests: Light (log(1+lux)) — Invaded minus Control =====\n")
  print(t_light_tbl, n = Inf, digits = 4)



  # ---- Effect size plots (Welch t-tests, Invaded − Control) --------------------
  library(ggplot2)
  library(patchwork)

  # Colours for sites
  pal_site <- c(Low = "#1b9e77", Middle = "#d95f02", Top = "#7570b3")

  # A) Temperature (% change), horizontal, coloured by site (no legend)
  p_ttemp <- ggplot(t_temp_tbl, aes(x = pct_change, y = site, colour = site)) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = 100 * ci_low_C / mean_control_C,
                       xmax = 100 * ci_high_C / mean_control_C),
                   height = 0.18, linewidth = 0.9) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = pal_site) +
    labs(x = "Effect size (% Temperature change, Invaded – Control)", y = NULL) +
    theme_pub +
    theme(legend.position = "none") +
    annotate("text", x = -Inf, y = Inf, label = "A",
             fontface = "bold", hjust = -0.4, vjust = 1.6, size = 5)

  # B) Light (% change), horizontal, coloured by site (no legend)
  p_tlight <- ggplot(t_light_tbl, aes(x = pct_change, y = site, colour = site)) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = pct_lwr, xmax = pct_upr),
                   height = 0.18, linewidth = 0.9) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = pal_site) +
    labs(x = "Effect size (% Light change, Invaded – Control)", y = NULL) +
    theme_pub +
    theme(legend.position = "none") +
    annotate("text", x = -Inf, y = Inf, label = "B",
             fontface = "bold", hjust = -0.4, vjust = 1.6, size = 5)

  # Stack vertically (no plot_annotation)
  effect_t_tests <- p_ttemp / p_tlight

  print(effect_t_tests)

  # Optional save
  # ggsave("EffectSizes_Temp_Light_tTests.png", effect_t_tests, width = 7, height = 6, dpi_

