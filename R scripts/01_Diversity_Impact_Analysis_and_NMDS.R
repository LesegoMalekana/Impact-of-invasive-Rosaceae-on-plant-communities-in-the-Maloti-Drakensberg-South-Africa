# =========================
# Impact figure: diversity metrics and NMDS
# October 2022 only
# =========================

# Step 0: Packages
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(emmeans)
library(MASS)
library(car)
library(stringr)
library(openxlsx)
library(multcomp)
library(multcompView)

# Step 1: Colours and theme
pal_trt <- c(Control = "#0072B2", Invaded = "#D55E00")

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
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

# Step 2: Load and prepare data
file_path <- "Your file location of Diversity_Impact_Analysis_and_NMDS"

data <- readr::read_csv(file_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    Site = factor(Site, levels = c("Low", "Middle", "Top")),
    Treatment = dplyr::if_else(Treatment == "Uncleared", "Invaded", Treatment),
    Treatment = factor(Treatment, levels = c("Control", "Invaded")),
    Date = factor(Date),
    Plot_ID = factor(Plot_ID)
  ) %>%
  dplyr::filter(Date == "October_2022", !is.na(Treatment))

# Step 3: Build species matrix and metadata
data_long <- data %>%
  tidyr::pivot_longer(
    cols = c(Natives_sp, Alien_sp),
    names_to = "Species_Type",
    values_to = "Species"
  ) %>%
  tidyr::pivot_longer(
    cols = c(Natives_nr, Alien_nr),
    names_to = "Abundance_Type",
    values_to = "Abundance"
  ) %>%
  dplyr::filter(!is.na(Species), !is.na(Abundance)) %>%
  dplyr::mutate(Abundance = as.numeric(Abundance))

species_matrix_2022 <- data_long %>%
  tidyr::complete(Plot_ID, Species, fill = list(Abundance = 0)) %>%
  tidyr::pivot_wider(
    id_cols = Plot_ID,
    names_from = Species,
    values_from = Abundance,
    values_fn = sum
  ) %>%
  tibble::column_to_rownames("Plot_ID") %>%
  dplyr::select(dplyr::where(is.numeric))

metadata_2022 <- data %>%
  dplyr::select(Plot_ID, Site, Treatment) %>%
  dplyr::distinct() %>%
  dplyr::filter(Plot_ID %in% rownames(species_matrix_2022)) %>%
  dplyr::mutate(Plot_ID = as.character(Plot_ID))

species_matrix_2022 <- species_matrix_2022[metadata_2022$Plot_ID, , drop = FALSE]
stopifnot(identical(rownames(species_matrix_2022), metadata_2022$Plot_ID))

# Step 4: NMDS ordination
set.seed(42)

nmds_3d <- vegan::metaMDS(
  species_matrix_2022,
  distance = "bray",
  k = 3,
  trymax = 100
)

nmds_scores <- as.data.frame(vegan::scores(nmds_3d, display = "sites")) %>%
  tibble::rownames_to_column("Plot_ID") %>%
  dplyr::left_join(metadata_2022, by = "Plot_ID")

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_point(size = 2.5, alpha = 0.9) +
  stat_ellipse(type = "t", linewidth = 0.7) +
  facet_wrap(~ Site, nrow = 1) +
  scale_colour_manual(values = pal_trt) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Treatment",
    title = paste0(
      "NMDS ordination, Bray-Curtis dissimilarity; stress = ",
      round(nmds_3d$stress, 3)
    )
  ) +
  theme_pub +
  theme(
    strip.text = element_text(face = "bold")
  )

print(p_nmds)
cat("NMDS stress:", round(nmds_3d$stress, 3), "\n")

# Step 5: PERMANOVA and dispersion checks
adonis_main <- vegan::adonis2(
  species_matrix_2022 ~ Treatment + Site,
  data = metadata_2022,
  method = "bray",
  permutations = 999,
  by = "margin"
)

adonis_int <- vegan::adonis2(
  species_matrix_2022 ~ Treatment * Site,
  data = metadata_2022,
  method = "bray",
  permutations = 999,
  by = "margin"
)

dist_matrix <- vegan::vegdist(species_matrix_2022, method = "bray")
bd_treatment <- vegan::betadisper(dist_matrix, metadata_2022$Treatment)
bd_site <- vegan::betadisper(dist_matrix, metadata_2022$Site)

disp_treatment <- vegan::permutest(bd_treatment)
disp_site <- vegan::permutest(bd_site)

print(adonis_main)
print(adonis_int)
print(disp_treatment)
print(disp_site)

# Step 6: Diversity metrics
H <- vegan::diversity(species_matrix_2022, index = "shannon")
S <- vegan::specnumber(species_matrix_2022)
J <- ifelse(S > 0, H / log(S), NA_real_)
N <- rowSums(species_matrix_2022)

div_df <- metadata_2022 %>%
  dplyr::mutate(
    Richness = S[Plot_ID],
    Shannon = H[Plot_ID],
    Evenness = J[Plot_ID],
    Abundance = N[Plot_ID]
  )

# Step 7: Beta diversity
calc_beta_pairwise_long <- function(pa_mat, meta) {
  d <- vegan::vegdist(pa_mat, method = "bray", binary = TRUE)
  D <- as.matrix(d)

  meta2 <- meta %>%
    dplyr::mutate(
      Plot_ID = as.character(Plot_ID),
      Site = as.character(Site),
      Treatment = as.character(Treatment)
    )

  out <- list()

  for (g in unique(meta2$Site)) {
    for (tr in c("Control", "Invaded")) {
      plots <- meta2$Plot_ID[meta2$Site == g & meta2$Treatment == tr]
      plots <- intersect(plots, rownames(D))

      if (length(plots) >= 2) {
        sub <- D[plots, plots, drop = FALSE]
        vals <- sub[lower.tri(sub)]

        out[[length(out) + 1]] <- data.frame(
          Site = g,
          Treatment = tr,
          d = vals
        )
      }
    }
  }

  dplyr::bind_rows(out)
}

pa_mat <- (species_matrix_2022 > 0) * 1

beta_long <- calc_beta_pairwise_long(pa_mat, metadata_2022) %>%
  dplyr::mutate(
    Site = factor(Site, levels = c("Low", "Middle", "Top")),
    Treatment = factor(Treatment, levels = c("Control", "Invaded"))
  ) %>%
  dplyr::filter(is.finite(d))

# Step 8: Fit models
mod_rich <- MASS::glm.nb(Richness ~ Treatment * Site, data = div_df)
mod_abun <- MASS::glm.nb(Abundance ~ Treatment * Site, data = div_df)

mod_shan <- stats::lm(Shannon ~ Treatment * Site, data = div_df)
mod_even <- stats::lm(Evenness ~ Treatment * Site, data = div_df)

mod_beta <- stats::aov(d ~ Treatment * Site, data = beta_long)

# Step 9: Overall model tests
anova_rich <- car::Anova(mod_rich, type = 2)
anova_abun <- car::Anova(mod_abun, type = 2)
anova_shan <- car::Anova(mod_shan, type = 2)
anova_even <- car::Anova(mod_even, type = 2)
anova_beta <- car::Anova(mod_beta, type = 2)

print(anova_rich)
print(anova_abun)
print(anova_shan)
print(anova_even)
print(anova_beta)

# Step 10: Helper functions
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
    stop("Could not find p-value column in ANOVA table.")
  }

  tab$term <- rownames(tab)

  p_trt  <- tab[[p_col]][match("Treatment", tab$term)]
  p_site <- tab[[p_col]][match("Site", tab$term)]
  p_int  <- tab[[p_col]][match("Treatment:Site", tab$term)]

  paste0(
    "Site ", p_to_stars(p_site),
    ", Treatment ", p_to_stars(p_trt),
    ", Site × Treatment ", p_to_stars(p_int)
  )
}

get_cld_df <- function(model, metric_name, type_response = FALSE) {
  if (type_response) {
    em <- emmeans::emmeans(model, ~ Site * Treatment, type = "response")
  } else {
    em <- emmeans::emmeans(model, ~ Site * Treatment)
  }

  cld_df <- multcomp::cld(
    em,
    Letters = letters,
    adjust = "bonferroni"
  ) %>%
    as.data.frame()

  cld_df %>%
    dplyr::mutate(
      .group = stringr::str_trim(.group),
      Metric = metric_name
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

get_site_contrasts <- function(model, metric_name, type_response = FALSE) {
  if (type_response) {
    em <- emmeans::emmeans(model, ~ Treatment | Site, type = "response")
  } else {
    em <- emmeans::emmeans(model, ~ Treatment | Site)
  }

  emmeans::contrast(em, method = "pairwise", adjust = "bonferroni") %>%
    as.data.frame() %>%
    dplyr::mutate(Metric = metric_name, .before = 1)
}

# Step 11: Compact letters
cld_rich <- get_cld_df(mod_rich, "Species richness", type_response = TRUE)
cld_abun <- get_cld_df(mod_abun, "Species abundance", type_response = TRUE)
cld_shan <- get_cld_df(mod_shan, "Shannon diversity", type_response = FALSE)
cld_even <- get_cld_df(mod_even, "Species evenness", type_response = FALSE)
cld_beta <- get_cld_df(mod_beta, "Beta diversity", type_response = FALSE)

cld_all <- dplyr::bind_rows(cld_rich, cld_abun, cld_shan, cld_even, cld_beta) %>%
  dplyr::mutate(
    Site = factor(Site, levels = c("Low", "Middle", "Top")),
    Treatment = factor(Treatment, levels = c("Control", "Invaded"))
  )

# Step 12: Plot data
alpha_long <- div_df %>%
  dplyr::select(Plot_ID, Site, Treatment, Richness, Abundance, Shannon, Evenness) %>%
  tidyr::pivot_longer(
    cols = c(Richness, Abundance, Shannon, Evenness),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      Richness  = "Species richness",
      Abundance = "Species abundance",
      Shannon   = "Shannon diversity",
      Evenness  = "Species evenness"
    )
  )

beta_plot <- beta_long %>%
  dplyr::transmute(
    Site,
    Treatment,
    Metric = "Beta diversity",
    Value = d
  )

plot_df <- dplyr::bind_rows(alpha_long, beta_plot) %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c(
        "Species richness",
        "Species abundance",
        "Shannon diversity",
        "Species evenness",
        "Beta diversity"
      )
    ),
    Site = factor(Site, levels = c("Low", "Middle", "Top")),
    Treatment = factor(Treatment, levels = c("Control", "Invaded"))
  )

# Step 13: Label positions
y_pos <- plot_df %>%
  dplyr::group_by(Metric, Site, Treatment) %>%
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
  dplyr::left_join(y_pos, by = c("Metric", "Site", "Treatment")) %>%
  dplyr::left_join(metric_pad, by = "Metric") %>%
  dplyr::mutate(y = y + pad)

# Step 14: Panel annotation text
ann_df <- tibble::tibble(
  Metric = c(
    "Species richness",
    "Species abundance",
    "Shannon diversity",
    "Species evenness",
    "Beta diversity"
  ),
  label = c(
    make_sig_label(anova_rich),
    make_sig_label(anova_abun),
    make_sig_label(anova_shan),
    make_sig_label(anova_even),
    make_sig_label(anova_beta)
  )
)

# Step 15: Diversity figure
p_main <- ggplot(plot_df, aes(x = Site, y = Value, fill = Treatment)) +
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
    aes(x = Site, y = y, group = Treatment, label = .group),
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
    legend.position = "right",
    plot.margin = margin(12, 60, 12, 12)
  )

print(p_main)

# Step 16: Export tables
anova_export <- dplyr::bind_rows(
  anova_to_df(anova_rich, "Species richness"),
  anova_to_df(anova_abun, "Species abundance"),
  anova_to_df(anova_shan, "Shannon diversity"),
  anova_to_df(anova_even, "Species evenness"),
  anova_to_df(anova_beta, "Beta diversity")
)

contrasts_export <- dplyr::bind_rows(
  get_site_contrasts(mod_rich, "Species richness", type_response = TRUE),
  get_site_contrasts(mod_abun, "Species abundance", type_response = TRUE),
  get_site_contrasts(mod_shan, "Shannon diversity", type_response = FALSE),
  get_site_contrasts(mod_even, "Species evenness", type_response = FALSE),
  get_site_contrasts(mod_beta, "Beta diversity", type_response = FALSE)
)

cld_export <- cld_all %>%
  dplyr::select(Metric, Site, Treatment, .group)

out_xlsx <- " If you want to save the outputs - Impact_model_results_and_letters.xlsx"

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Model_tests")
openxlsx::addWorksheet(wb, "Within_site_contrasts")
openxlsx::addWorksheet(wb, "Compact_letters")
openxlsx::addWorksheet(wb, "PERMANOVA")
openxlsx::addWorksheet(wb, "Dispersion_tests")

openxlsx::writeData(wb, "Model_tests", anova_export)
openxlsx::writeData(wb, "Within_site_contrasts", contrasts_export)
openxlsx::writeData(wb, "Compact_letters", cld_export)
openxlsx::writeData(wb, "PERMANOVA", as.data.frame(adonis_int))
openxlsx::writeData(wb, "Dispersion_tests", list(
  Treatment = capture.output(print(disp_treatment)),
  Site = capture.output(print(disp_site))
))

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Excel file saved to:\n", out_xlsx, "\n")

# Step 17: Export figures
ggsave(
  filename = " Save the picture - Figure_Impact_Diversity_Metrics.png",
  plot = p_main,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Save the picture - Figure_Impact_Diversity_Metrics.pdf",
  plot = p_main,
  width = 10,
  height = 14,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "Save the picture - Figure_NMDS_2022.png",
  plot = p_nmds,
  width = 8.5,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Save the picture - Figure_NMDS_2022.pdf",
  plot = p_nmds,
  width = 8.5,
  height = 4.8,
  units = "in",
  bg = "white"
)

