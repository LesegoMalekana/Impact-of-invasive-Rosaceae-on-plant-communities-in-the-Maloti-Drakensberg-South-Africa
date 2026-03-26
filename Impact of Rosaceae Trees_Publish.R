# Step 0: Libraries and theme
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ggrepel)
library(emmeans)

# Colours and plot theme
pal_trt <- c(Control = "#7CAE00", Invaded = "#F8766D")
pal_metric <- c("Species Richness" = "#1b9e77",   # teal
                "Species Abundance" = "#d95f02",  # orange
                "Shannon Diversity" = "#7570b3",  # purple
                "Species Evenness"  = "#e7298a",  # pink
                "Beta Diversity"    = "#1f78b4")  # blue

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    legend.position = "right"
  )

# Step 1: Load and prepare data
file_path <- "Your local file location"
data <- read_csv(file_path, show_col_types = FALSE) %>%
  mutate(
    Site      = factor(Site),                 # expected: Low / Middle / Top
    Treatment = recode(Treatment, Uncleared = "Invaded"),
    Treatment = factor(Treatment, levels = c("Control","Invaded")),
    Date      = factor(Date),
    Plot_ID   = factor(Plot_ID)
  ) %>%
  filter(Date == "October_2022", !is.na(Treatment))

# Step 2: Build species matrix (plots × species; abundance) and metadata
data_long <- data %>%
  pivot_longer(cols = c(Natives_sp, Alien_sp),
               names_to = "Species_Type", values_to = "Species") %>%
  pivot_longer(cols = c(Natives_nr, Alien_nr),
               names_to = "Abundance_Type", values_to = "Abundance") %>%
  filter(!is.na(Species), !is.na(Abundance)) %>%
  mutate(Abundance = as.numeric(Abundance))

species_matrix_2022 <- data_long %>%
  complete(Plot_ID, Species, fill = list(Abundance = 0)) %>%
  pivot_wider(id_cols = Plot_ID,
              names_from = Species,
              values_from = Abundance,
              values_fn = sum) %>%
  column_to_rownames("Plot_ID") %>%
  select(where(is.numeric))

metadata_2022 <- data %>%
  select(Plot_ID, Site, Treatment) %>%
  distinct() %>%
  filter(Plot_ID %in% rownames(species_matrix_2022)) %>%
  mutate(
    Plot_ID = as.character(Plot_ID),
    Site = factor(Site, levels = c("Low","Middle","Top"))
  )

species_matrix_2022 <- species_matrix_2022[metadata_2022$Plot_ID, , drop = FALSE]
stopifnot(identical(rownames(species_matrix_2022), metadata_2022$Plot_ID))

# Step 3: NMDS ordinations (k=3) and 2D projections by Site
set.seed(42)
nmds_3d <- metaMDS(species_matrix_2022, distance = "bray", k = 3, trymax = 100)

plot_nmds_by_site <- function(nmds_result, dim_x = 1, dim_y = 2, title_text = "", meta = metadata_2022) {
  scr <- as.data.frame(scores(nmds_result, display = "sites"))
  scr$Plot_ID <- rownames(scr)
  df <- left_join(scr, meta, by = "Plot_ID")

  ggplot(df, aes_string(x = paste0("NMDS", dim_x),
                        y = paste0("NMDS", dim_y),
                        colour = "Treatment")) +
    geom_point(size = 2.5) +
    stat_ellipse(type = "t", linewidth = 0.6) +
    facet_wrap(~ Site) +
    scale_colour_manual(values = pal_trt) +
    labs(
      title = paste0(title_text, "  (Stress = ", round(nmds_result$stress, 3), ")"),
      x = paste0("NMDS", dim_x),
      y = paste0("NMDS", dim_y)
    ) +
    theme_pub
}

p_nmds12 <- plot_nmds_by_site(nmds_3d, 1, 2, "NMDS1 vs NMDS2")
p_nmds13 <- plot_nmds_by_site(nmds_3d, 1, 3, "NMDS1 vs NMDS3")
p_nmds23 <- plot_nmds_by_site(nmds_3d, 2, 3, "NMDS2 vs NMDS3")
print(p_nmds12); print(p_nmds13); print(p_nmds23)
cat("NMDS (3D) stress:", round(nmds_3d$stress, 3), "\n")

cat("NMDS (3D) stress:", formatC(nmds_3d$stress, format = "f", digits = 3), "\n")

# Step 4: PERMANOVA (+ interaction) and homogeneity of dispersions
adonis_main <- adonis2(
  species_matrix_2022 ~ Treatment + Site,
  data = metadata_2022,
  method = "bray", permutations = 999, by = "margin"
)
adonis_int <- adonis2(
  species_matrix_2022 ~ Treatment * Site,
  data = metadata_2022,
  method = "bray", permutations = 999, by = "margin"
)
print(adonis_main)
print(adonis_int)

dist_matrix <- vegdist(species_matrix_2022, method = "bray")
bd_Treat <- betadisper(dist_matrix, metadata_2022$Treatment)
bd_Site  <- betadisper(dist_matrix, metadata_2022$Site)
cat("\nHomogeneity of dispersions (Treatment):\n"); print(permutest(bd_Treat))
cat("\nHomogeneity of dispersions (Site):\n");      print(permutest(bd_Site))

# Step 5: Diversity metrics per plot (alpha diversity and abundance)
H  <- diversity(species_matrix_2022, index = "shannon")     # ln-based
S  <- specnumber(species_matrix_2022)
J  <- ifelse(S > 0, H / log(S), NA_real_)                   # Pielou’s evenness
N  <- rowSums(species_matrix_2022)

div_df <- metadata_2022 %>%
  mutate(
    Richness = S[Plot_ID],
    Shannon  = H[Plot_ID],
    Evenness = J[Plot_ID],
    Abundance = N[Plot_ID]
  )
div_df2 <- div_df %>%
  mutate(ElevGroup = as.character(Site),
         Treatment = as.character(Treatment))

# Step 6: ANOVAs with interaction (Richness, Abundance, Shannon, Evenness, Beta)

do_anova_interaction <- function(metric, df) {
  f <- as.formula(paste0(metric, " ~ Treatment * Site"))
  fit <- aov(f, data = df)
  cat("\n===== ANOVA for", metric, "=====\n")
  print(summary(fit))
  em <- emmeans(fit, ~ Treatment | Site)
  print(contrast(em, method = "pairwise", adjust = "bonferroni"))
  invisible(fit)
}
fit_R <- do_anova_interaction("Richness",  div_df)
fit_A <- do_anova_interaction("Abundance", div_df)
fit_H <- do_anova_interaction("Shannon",   div_df)
fit_J <- do_anova_interaction("Evenness",  div_df)

# Beta diversity ANOVA using within-treatment pairwise dissimilarities as observations
calc_beta_pairwise_long <- function(pa_mat, meta) {
  d <- vegdist(pa_mat, method = "bray", binary = TRUE)   # Sørensen dissimilarity
  D <- as.matrix(d)
  meta <- meta %>%
    mutate(Plot_ID = as.character(Plot_ID),
           ElevGroup = as.character(Site),
           Treatment = as.character(Treatment))
  out <- list()
  for (g in unique(meta$ElevGroup)) {
    for (tr in c("Control","Invaded")) {
      plots <- meta$Plot_ID[meta$ElevGroup == g & meta$Treatment == tr]
      plots <- intersect(plots, rownames(D))
      if (length(plots) >= 2) {
        sub <- D[plots, plots, drop = FALSE]
        vals <- sub[lower.tri(sub)]
        out[[length(out)+1]] <- data.frame(
          ElevGroup = g, Site = g, Treatment = tr, d = vals
        )
      }
    }
  }
  bind_rows(out)
}
pa_mat <- (species_matrix_2022 > 0) * 1
beta_long <- calc_beta_pairwise_long(pa_mat, metadata_2022) %>% filter(is.finite(d))

fit_B <- aov(d ~ Treatment * Site, data = beta_long)
cat("\n===== ANOVA for Beta Diversity (pairwise dissimilarities) =====\n")
print(summary(fit_B))
em_B <- emmeans(fit_B, ~ Treatment | Site)
print(contrast(em_B, method = "pairwise", adjust = "bonferroni"))

# Step 7: Effect plots using LRR (Richness, Abundance, Shannon, Evenness, Beta)

# Step 7: Effect sizes (log-response ratios for alpha & beta diversity)

# Helper: compute LRR + CI for alpha metrics
calc_LRR_alpha <- function(df, metric) {
  df %>%
    group_by(ElevGroup, Treatment) %>%
    summarise(mean_val = mean(.data[[metric]], na.rm = TRUE),
              sd_val   = sd(.data[[metric]], na.rm = TRUE),
              n        = sum(!is.na(.data[[metric]])), .groups = "drop") %>%
    pivot_wider(names_from = Treatment, values_from = c(mean_val, sd_val, n)) %>%
    mutate(
      LRR   = log(mean_val_Invaded / mean_val_Control),
      SE    = sqrt( (sd_val_Control^2 / (pmax(n_Control,1) * mean_val_Control^2)) +
                      (sd_val_Invaded^2 / (pmax(n_Invaded,1) * mean_val_Invaded^2)) ),
      CI_lo = LRR - 1.96 * SE,
      CI_hi = LRR + 1.96 * SE
    ) %>%
    transmute(ElevGroup,
              Metric = metric,
              LRR, CI_lo, CI_hi)
}

# Alpha metrics
lrr_rich    <- calc_LRR_alpha(div_df2, "Richness")   %>% mutate(Metric = "Species Richness")
lrr_abun    <- calc_LRR_alpha(div_df2, "Abundance")  %>% mutate(Metric = "Species Abundance")
lrr_shannon <- calc_LRR_alpha(div_df2, "Shannon")    %>% mutate(Metric = "Shannon Diversity")
lrr_even    <- calc_LRR_alpha(div_df2, "Evenness")   %>% mutate(Metric = "Species Evenness")

# Beta diversity (Sørensen dissimilarity means)
beta_means <- calc_beta_means(pa_mat, metadata_2022)
lrr_beta <- beta_means %>%
  mutate(
    LRR   = log(mean_d_Invaded / mean_d_Control),
    SE    = sqrt((sd_d_Control^2 / (pmax(n_pairs_Control,1) * mean_d_Control^2)) +
                   (sd_d_Invaded^2 / (pmax(n_pairs_Invaded,1) * mean_d_Invaded^2))),
    CI_lo = LRR - 1.96 * SE,
    CI_hi = LRR + 1.96 * SE,
    Metric = "Beta Diversity"
  ) %>%
  select(ElevGroup, Metric, LRR, CI_lo, CI_hi)

# Combine
lrr_all <- bind_rows(lrr_rich, lrr_abun, lrr_shannon, lrr_even, lrr_beta) %>%
  mutate(
    ElevGroup = factor(ElevGroup, levels = c("Low","Middle","Top")),
    Metric    = factor(Metric, levels = c("Species Richness","Species Abundance",
                                          "Shannon Diversity","Species Evenness","Beta Diversity"))
  )

# Distinct palette (avoid confusion between Richness & Beta)
pal_metric <- c("Species Richness" = "#1b9e77",   # teal
                "Species Abundance" = "#d95f02",  # orange
                "Shannon Diversity" = "#7570b3",  # purple
                "Species Evenness"  = "#e7298a",  # pink
                "Beta Diversity"    = "#1f78b4")  # blue

# Plot

p_lrr <- ggplot(lrr_all, aes(y = ElevGroup, x = LRR, colour = Metric)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.2, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "solid", colour = "black") +
  geom_hline(yintercept = c(1.5, 2.5), linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = pal_metric) +
  labs(
    title = "Effect of invasion (log-response ratios)",
    subtitle = "LRR = ln(Invaded / Control); LRR < 0 indicates loss under invasion",
    x = "Effect of Invasion",
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank()
  )

print(p_lrr)



# Step 8: Scatterplot of paired richness loss vs Sørensen similarity

# Sørensen dissimilarity (binary Bray–Curtis)
sor_dist <- vegdist(pa_mat, method = "bray", binary = TRUE)
sor_mat  <- as.matrix(sor_dist)
rn <- rownames(sor_mat)

# Build all Control–Invaded pairs within each Site
pair_df <- metadata_2022 %>%
  transmute(Plot_ID, Site = as.character(Site),
            ElevGroup = as.character(Site),
            Treatment = as.character(Treatment)) %>%
  group_by(Site, ElevGroup) %>%
  summarise(pairs = list(tidyr::crossing(
    Control = Plot_ID[Treatment == "Control"],
    Invaded = Plot_ID[Treatment == "Invaded"]
  )), .groups = "drop") %>%
  tidyr::unnest(pairs) %>%
  mutate(
    Diss = mapply(function(i, j) {
      if (i %in% rn && j %in% rn) sor_mat[i, j] else NA_real_
    }, Control, Invaded),
    Sorensen_pct = (1 - Diss) * 100
  ) %>%
  filter(is.finite(Sorensen_pct))

# Add richness values
plot_rich <- div_df2 %>%
  transmute(Plot_ID = as.character(Plot_ID),
            ElevGroup = as.character(ElevGroup),
            Treatment = as.character(Treatment),
            Richness  = as.numeric(Richness))

pair_compare <- pair_df %>%
  left_join(plot_rich %>% select(Plot_ID, Rich_C = Richness, ElevGroup),
            by = c("Control" = "Plot_ID", "ElevGroup")) %>%
  left_join(plot_rich %>% select(Plot_ID, Rich_I = Richness),
            by = c("Invaded" = "Plot_ID")) %>%
  mutate(
    Rich_loss_pct = (Rich_C - Rich_I) / pmax(Rich_C, 1e-9) * 100
  ) %>%
  filter(is.finite(Rich_loss_pct))

# Scatterplot
p_pair <- ggplot(pair_compare,
                 aes(x = Rich_loss_pct, y = Sorensen_pct, colour = ElevGroup)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  scale_colour_brewer(palette = "Set2") +
  labs(
    title = "Plot-pair: Richness loss vs Sørensen similarity",
    x = "Richness loss per Control–Invaded pair (%)",
    y = "Sørensen similarity per pair (%)",
    colour = "Elevation group"
  ) + theme_pub

print(p_pair)

library(broom)
library(purrr)
library(dplyr)
library(tidyr)

# Fit separate linear models per ElevGroup
lm_results <- pair_compare %>%
  group_by(ElevGroup) %>%
  nest() %>%
  mutate(
    model  = map(data, ~ lm(Sorensen_pct ~ Rich_loss_pct, data = .x)),
    tidy   = map(model, broom::tidy),
    glance = map(model, broom::glance)
  )

# Extract slopes
slopes <- lm_results %>%
  unnest(tidy) %>%
  filter(term == "Rich_loss_pct") %>%
  select(ElevGroup, slope = estimate, std.error, p.value)

# Extract R² values
r2_vals <- lm_results %>%
  unnest(glance) %>%
  select(ElevGroup, r.squared, adj.r.squared)

# Combine slope + R²
lm_summary <- left_join(slopes, r2_vals, by = "ElevGroup")

print(lm_summary)



# ---- Step X1: Indicator species per elevation (Control vs Invaded) ----
suppressPackageStartupMessages({
  library(indicspecies)
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(stringr)
})

# Helper to run IndVal within one elevation level
run_indval_one_site <- function(site_level, comm_mat, meta, nperm = 999, presence_absence = FALSE) {
  rows <- which(meta$Site == site_level)
  if (length(rows) < 4) return(NULL)  # need a few plots

  comm <- comm_mat[rows, , drop = FALSE]
  grp  <- droplevels(meta$Treatment[rows])

  # Drop species with zero total in this stratum
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  if (ncol(comm) == 0) return(NULL)

  # Optional: presence/absence (toggle)
  if (presence_absence) comm[] <- (comm > 0) * 1

  set.seed(123)
  res <- multipatt(comm, grp, func = "IndVal.g", duleg = TRUE,
                   control = how(nperm = nperm))

  sumres <- summary(res, alpha = 0.05, indvalcomp = TRUE)
  sig <- sumres$sign
  if (is.null(sig) || nrow(sig) == 0) {
    return(tibble(Site = site_level, Species = character(0),
                  Associated_with = character(0), IndVal = numeric(0),
                  Sensitivity_s = numeric(0), Specificity_r = numeric(0),
                  p_value = numeric(0)))
  }

  # Tidy
  sig_df <- sig %>% as.data.frame() %>% rownames_to_column("Species")
  metric_cols <- c("Species","stat","p.value","s","r","index")
  comb_cols <- setdiff(names(sig_df), metric_cols)

  # Map numeric columns ("1","2") to group labels (Control/Invaded)
  idx_num <- suppressWarnings(as.integer(comb_cols))
  if (all(!is.na(idx_num))) {
    lvl <- levels(grp)
    names(sig_df)[match(comb_cols, names(sig_df))] <- lvl[idx_num]
    comb_cols <- lvl[idx_num]
  }

  out <- sig_df %>%
    filter(p.value <= 0.05) %>%
    rowwise() %>%
    mutate(Associated_with = paste(comb_cols[as.logical(c_across(all_of(comb_cols)))], collapse = " + ")) %>%
    ungroup() %>%
    transmute(
      Site          = site_level,
      Species,
      Associated_with,
      IndVal        = stat,
      Sensitivity_s = s,
      Specificity_r = r,
      p_value       = p.value
    ) %>%
    arrange(p_value, desc(IndVal))

  out
}

# Run for all Site levels in your metadata_2022
indval_results <- levels(metadata_2022$Site) %>%
  map(~ run_indval_one_site(.x, species_matrix_2022, metadata_2022,
                            nperm = 999, presence_absence = FALSE)) %>%
  list_rbind()

# Quick peek (top 10 per elevation)
indval_results %>%
  group_by(Site) %>%
  slice_min(p_value, n = 10, with_ties = FALSE) %>%
  ungroup()

# Optional export
# write.csv(indval_results, "IndVal_Oct2022_byElevation.csv", row.names = FALSE)





library(dplyr); library(tidyr); library(purrr); library(tibble)

# Presence-absence matrix by Site × Treatment
presence_pa <- as.data.frame((species_matrix_2022 > 0) * 1)
presence_pa$Plot_ID <- rownames(species_matrix_2022)
pres_long <- presence_pa %>%
  pivot_longer(-Plot_ID, names_to = "Species", values_to = "Present") %>%
  left_join(metadata_2022, by = "Plot_ID")

# Helper: for one elevation, which spp are in Control but not in Invaded?
lost_in_invaded_one <- function(site_level) {
  d <- pres_long %>% filter(Site == site_level)
  ctl <- d %>% filter(Treatment == "Control") %>%
    group_by(Species) %>% summarise(n_ctl = sum(Present), .groups = "drop")
  inv <- d %>% filter(Treatment == "Invaded") %>%
    group_by(Species) %>% summarise(n_inv = sum(Present), .groups = "drop")
  full <- ctl %>% full_join(inv, by = "Species") %>%
    mutate(n_ctl = coalesce(n_ctl, 0L), n_inv = coalesce(n_inv, 0L),
           Lost_in_invaded = (n_ctl > 0 & n_inv == 0))
  full %>%
    filter(Lost_in_invaded) %>%
    arrange(desc(n_ctl)) %>%
    mutate(Site = site_level, .before = 1)
}

lost_lists <- levels(metadata_2022$Site) %>%
  map_df(lost_in_invaded_one)

# Top lost species per elevation (rank by how many Control plots they occurred in)
lost_lists %>%
  group_by(Site) %>%
  slice_max(n_ctl, n = 10, with_ties = FALSE) %>%
  ungroup()


library(betapart)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

beta_by_site <- map_df(levels(metadata_2022$Site), function(site_level) {
  rows <- which(metadata_2022$Site == site_level)
  if (length(rows) < 2) {
    return(tibble(Site = site_level, beta_jac = NA_real_, beta_jtu = NA_real_, beta_jne = NA_real_,
                  n_pairs_between = 0L, note = "Too few plots"))
  }
  # Subset community + treatments
  comm <- (species_matrix_2022[rows, , drop = FALSE] > 0) * 1
  trt  <- droplevels(metadata_2022$Treatment[rows])

  # Drop all-zero species
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  if (ncol(comm) == 0) {
    return(tibble(Site = site_level, beta_jac = NA_real_, beta_jtu = NA_real_, beta_jne = NA_real_,
                  n_pairs_between = 0L, note = "All-zero species after subsetting"))
  }
  # Drop all-zero plots (keep trt aligned)
  keep_rows <- rowSums(comm) > 0
  comm <- comm[keep_rows, , drop = FALSE]
  trt  <- trt[keep_rows]

  # Quick counts after filtering
  cnt <- table(trt)
  n_ctl <- unname(if ("Control" %in% names(cnt)) cnt[["Control"]] else 0L)
  n_inv <- unname(if ("Invaded" %in% names(cnt)) cnt[["Invaded"]] else 0L)

  # Need at least one Control and one Invaded
  if (n_ctl == 0L || n_inv == 0L) {
    return(tibble(Site = site_level, beta_jac = NA_real_, beta_jtu = NA_real_, beta_jne = NA_real_,
                  n_pairs_between = 0L,
                  note = sprintf("Only one treatment present after filtering (Control=%d, Invaded=%d)", n_ctl, n_inv)))
  }

  # Pairwise Jaccard partition (presence/absence)
  part <- beta.pair(comm, index.family = "jaccard")  # dist objects

  # Convert to matrices for masking
  mjac <- as.matrix(part$beta.jac)
  mjtu <- as.matrix(part$beta.jtu)
  mjne <- as.matrix(part$beta.jne)

  # Mask: only i<j AND between-treatment
  upper <- row(mjac) < col(mjac)
  between <- outer(as.character(trt), as.character(trt), `!=`)
  mask <- upper & between

  n_pairs_between <- sum(mask, na.rm = TRUE)
  if (n_pairs_between == 0) {
    return(tibble(Site = site_level, beta_jac = NA_real_, beta_jtu = NA_real_, beta_jne = NA_real_,
                  n_pairs_between = 0L, note = "No between-treatment pairs"))
  }

  tibble(
    Site = site_level,
    beta_jac = mean(mjac[mask], na.rm = TRUE),
    beta_jtu = mean(mjtu[mask], na.rm = TRUE),
    beta_jne = mean(mjne[mask], na.rm = TRUE),
    n_pairs_between = n_pairs_between,
    note = sprintf("After filter: Control=%d, Invaded=%d", n_ctl, n_inv)
  )
})

beta_by_site






library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Presence/absence matrix
presence_pa <- as.data.frame((species_matrix_2022 > 0) * 1)
presence_pa$Plot_ID <- rownames(species_matrix_2022)

# Add metadata
pres_long <- presence_pa %>%
  pivot_longer(-Plot_ID, names_to = "Species", values_to = "Present") %>%
  left_join(metadata_2022, by = "Plot_ID")

# Function to get species overlap summary for one elevation
species_overlap_one <- function(site_level) {
  d <- pres_long %>% filter(Site == site_level)

  ctl <- d %>%
    filter(Treatment == "Control") %>%
    group_by(Species) %>%
    summarise(n_ctl = sum(Present), .groups = "drop")

  inv <- d %>%
    filter(Treatment == "Invaded") %>%
    group_by(Species) %>%
    summarise(n_inv = sum(Present), .groups = "drop")

  full <- full_join(ctl, inv, by = "Species") %>%
    mutate(across(c(n_ctl, n_inv), ~replace_na(., 0L)),
           Status = case_when(
             n_ctl > 0 & n_inv > 0 ~ "Shared",
             n_ctl > 0 & n_inv == 0 ~ "Control-only (Lost in Invaded)",
             n_ctl == 0 & n_inv > 0 ~ "Invaded-only (New/Replacement)",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(Status)) %>%
    arrange(Status, desc(pmax(n_ctl, n_inv)), Species) %>%
    mutate(Site = site_level, .before = 1)

  full
}

# Run for all sites and combine
species_overlap <- map_df(levels(metadata_2022$Site), species_overlap_one)

# Peek
species_overlap %>%
  group_by(Site, Status) %>%
  summarise(n_species = n(), .groups = "drop")

# Optional exports
#
write_csv(species_overlap, "SpeciesOverlap_byElevation_Oct2022.csv")


library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(openxlsx)

# --- Presence–absence matrix and metadata join ---
presence_pa <- as.data.frame((species_matrix_2022 > 0) * 1)
presence_pa$Plot_ID <- rownames(species_matrix_2022)

pres_long <- presence_pa %>%
  pivot_longer(-Plot_ID, names_to = "Species", values_to = "Present") %>%
  left_join(metadata_2022, by = "Plot_ID")

# --- Function: classify species status per elevation ---
species_overlap_one <- function(site_level) {
  d <- pres_long %>% filter(Site == site_level)
  ctl <- d %>% filter(Treatment == "Control") %>%
    group_by(Species) %>% summarise(n_ctl = sum(Present), .groups = "drop")
  inv <- d %>% filter(Treatment == "Invaded") %>%
    group_by(Species) %>% summarise(n_inv = sum(Present), .groups = "drop")

  full_join(ctl, inv, by = "Species") %>%
    mutate(across(c(n_ctl, n_inv), ~replace_na(., 0L)),
           Status = case_when(
             n_ctl > 0 & n_inv > 0 ~ "Shared",
             n_ctl > 0 & n_inv == 0 ~ "Control-only (Lost in Invaded)",
             n_ctl == 0 & n_inv > 0 ~ "Invaded-only (New/Replacement)"
           )) %>%
    filter(!is.na(Status)) %>%
    mutate(Site = site_level, .before = 1) %>%
    arrange(Site, Status, desc(pmax(n_ctl, n_inv)), Species)
}

# --- Combine all elevations ---
species_overlap <- map_df(levels(metadata_2022$Site), species_overlap_one)

# --- Summary counts ---
summary_counts <- species_overlap %>%
  count(Site, Status) %>%
  group_by(Site) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup()

# --- Save to Excel ---
out_path <- "C:/Users/lmalekana/Desktop/SpeciesOverlap_Oct2022.xlsx"

wb <- createWorkbook()
addWorksheet(wb, "Full_Species_List")
addWorksheet(wb, "Summary")

writeData(wb, "Full_Species_List", species_overlap)
writeData(wb, "Summary", summary_counts)

saveWorkbook(wb, out_path, overwrite = TRUE)

cat("✅ Excel file saved to:", out_path, "\n")

