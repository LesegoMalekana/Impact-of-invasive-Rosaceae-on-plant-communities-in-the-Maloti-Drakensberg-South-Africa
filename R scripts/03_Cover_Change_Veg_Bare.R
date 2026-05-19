# ==========================
# Appendix figure: vegetation and bare ground cover
# 2022 only; Invaded balanced to 6 plots per site
# ==========================

# Step 0: Packages
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
})

# Step 1: Colours and theme
pal_trt <- c(Control = "#0072B2", Invaded = "#D55E00")

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 0.6),
    axis.line.y = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

# Step 2: Load data
in_file <- "Your location of Cover_Change_Perc"
sheet_use <- readxl::excel_sheets(in_file)[1]
raw <- readxl::read_excel(in_file, sheet = sheet_use)

# Step 3: Prepare 2022 cover data
dat_2022 <- raw %>%
  dplyr::mutate(
    Date = suppressWarnings(as.integer(Date)),
    Site = factor(Site, levels = c("Low", "Middle", "Top")),
    Treatment = stringr::str_to_title(trimws(Treatment)),
    Treatment = dplyr::case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Invaded" ~ "Invaded",
      Treatment == "Cleared" ~ "Invaded",
      TRUE ~ NA_character_
    ),
    Treatment = factor(Treatment, levels = c("Control", "Invaded")),
    Veg_sum = Native + Alien + Grasses
  ) %>%
  dplyr::filter(Date == 2022, !is.na(Treatment), !is.na(Site))

# Step 4: Balance to 6 plots per Site × Treatment
set.seed(42)
target_n <- 6

choose_plots <- function(df, target) {
  plots <- unique(df$Plot)
  sample(plots, size = min(length(plots), target), replace = FALSE)
}

chosen_tbl <- dat_2022 %>%
  dplyr::group_by(Site, Treatment) %>%
  dplyr::summarise(
    chosen_plots = list(choose_plots(dplyr::cur_data_all(), target_n)),
    .groups = "drop"
  ) %>%
  tidyr::unnest_longer(chosen_plots, values_to = "Plot") %>%
  dplyr::arrange(Site, Treatment, Plot)

dat_bal <- dat_2022 %>%
  dplyr::inner_join(chosen_tbl, by = c("Site", "Treatment", "Plot"))

print(
  dat_bal %>%
    dplyr::count(Site, Treatment, Plot) %>%
    dplyr::count(Site, Treatment) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = n, values_fill = 0)
)

# Step 5: Summary table
sum_cover <- dat_bal %>%
  dplyr::group_by(Site, Treatment) %>%
  dplyr::summarise(
    N = dplyr::n_distinct(Plot),
    vegetation_mean = mean(Veg_sum, na.rm = TRUE),
    vegetation_se = sd(Veg_sum, na.rm = TRUE) / sqrt(N),
    bare_ground_mean = mean(Bare_Ground, na.rm = TRUE),
    bare_ground_se = sd(Bare_Ground, na.rm = TRUE) / sqrt(N),
    .groups = "drop"
  )

print(sum_cover)

# Step 6: Vegetation cover plot
p_veg <- ggplot(dat_bal, aes(x = Site, y = Veg_sum, fill = Treatment)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = 16,
    outlier.alpha = 0.45,
    colour = "black"
  ) +
  geom_point(
    aes(colour = Treatment),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    alpha = 0.60,
    size = 1.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal_trt) +
  scale_colour_manual(values = pal_trt) +
  labs(
    x = NULL,
    y = "Vegetation cover (%)",
    fill = "Treatment"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Step 7: Bare ground plot
p_bare <- ggplot(dat_bal, aes(x = Site, y = Bare_Ground, fill = Treatment)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = 16,
    outlier.alpha = 0.45,
    colour = "black"
  ) +
  geom_point(
    aes(colour = Treatment),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    alpha = 0.60,
    size = 1.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal_trt) +
  scale_colour_manual(values = pal_trt) +
  labs(
    x = "Site",
    y = "Bare ground (%)",
    fill = "Treatment"
  ) +
  theme_pub

# Step 8: Combine plots
cover_fig <- (p_veg / p_bare) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(1, 1)
  )

cover_fig <- cover_fig &
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

print(cover_fig)

# Step 9: Export summary and chosen plots
out_dir <- dirname(in_file)
out_xlsx <- file.path(out_dir, "Cover_Summaries_For_Paper.xlsx")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "2022_Cover_Summary")
openxlsx::addWorksheet(wb, "2022_Chosen_Plots")

openxlsx::writeData(wb, "2022_Cover_Summary", sum_cover)
openxlsx::writeData(wb, "2022_Chosen_Plots", chosen_tbl)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# Step 10: Export figure
ggsave(
  filename = file.path(out_dir, "Appendix_Cover_Vegetation_BareGround.png"),
  plot = cover_fig,
  width = 8.5,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Appendix_Cover_Vegetation_BareGround.pdf"),
  plot = cover_fig,
  width = 8.5,
  height = 7.5,
  units = "in",
  bg = "white"
)

