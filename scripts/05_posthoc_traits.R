#########################################################################
# EXTRACT r3_pc1 DRAWS + MERGE WITH TRAIT PCs  (delta_PC1 version)
#########################################################################
# Builds a tidy, long-format (draw x species) table with the trait-
# covariance component of the delta_PC1 (temperature CHANGE) random
# slope alongside each species' TRAIT PC scores -- ready for a per-draw
# multiple regression. Same rationale as the original: trait PCs are
# orthogonal by construction, so a multiple regression on them isn't
# confounded by collinearity between raw traits.
#
# KEY CHANGE FROM ORIGINAL: r3_pc1 (the r_3_2 Stan parameter) now
# represents the trait-covariance component of sensitivity to a CHANGE
# in PC1 (delta_PC1), not the raw/absolute PC1 value. Interpretively,
# this is "how strongly does a species' trait profile predict how much
# MORE (or less) it switches per unit of environmental shift" rather
# than "...per unit of absolute environmental state."
#
# Uses shiny_app_data_delta.RDS (the delta-model analog of
# shiny_app_data.RDS) rather than re-deriving from
# stan_dirichlet_explicit_delta.RDS directly, for the same reason as the
# original script: r3_pc1 in that file should already be correctly
# extracted and species-indexed via standata() in the delta model's prep
# script.
#
# INPUT:  data/shiny_app_data_delta.RDS
# OUTPUT: model_data (in-memory tibble; not saved to file by default)
#########################################################################

library(tidyverse)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data_delta.RDS")
darling <- read.csv("~/Documents/Coral/LongitudinalZoox_local/data/Darling_coral_species_traits.csv")
darling <- darling %>% select(Species, Life.history) %>% mutate(clean.spec = gsub(" ", "_", Species))

#### ----------- r3_pc1 DRAWS (n_draws x n_species matrix) ----------- ####
# NOTE: this is the delta_PC1 trait-component slope, not raw PC1

r3_pc1 <- app_data$r3_pc3
# species identity is carried by POSITION, matching app_data$species_levels
colnames(r3_pc1) <- app_data$species_levels

#### ----------- TRAIT PCs (unchanged -- traits themselves didn't change) ----------- ####

trait_pc_cols <- grep("^PC[0-9]+$", colnames(app_data$trait_pca_lookup), value = TRUE)
message("Trait PCs available: ", paste(trait_pc_cols, collapse = ", "))

traits_sp <- app_data$trait_pca_lookup %>%
  select(clean.spec, all_of(trait_pc_cols))

#### ----------- MERGE: LONG FORMAT, ONE ROW PER (DRAW x SPECIES) ----------- ####

r3_pc1_long <- as_tibble(r3_pc1) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "clean.spec", values_to = "r3_pc1")

model_data <- r3_pc1_long %>%
  inner_join(traits_sp, by = "clean.spec") %>%
  left_join(darling, by = "clean.spec")

model_data %>% group_by(draw, Life.history) %>%
  drop_na() %>%
  mutate(mean_r3 = mean(r3_pc1)) %>%
  ggplot(aes(x = Life.history, y = mean_r3)) +
  geom_boxplot() +
  labs(y = "Trait component of \u0394PC1 (temperature change) slope")

missing_species <- setdiff(app_data$species_levels, unique(model_data$clean.spec))
if (length(missing_species) > 0) {
  message(length(missing_species), " species excluded from model_data (no trait PC data): ",
          paste(missing_species, collapse = ", "))
}

message("model_data ready: ", n_distinct(model_data$draw), " draws x ",
        n_distinct(model_data$clean.spec), " species with trait PC data (",
        nrow(model_data), " total rows).")

#### ----------- PER-DRAW REGRESSION ----------- ####

library(broom)

fits <- model_data %>%
  group_by(draw) %>%
  summarise(fit = list(lm(r3_pc1 ~ PC1 + PC2 + PC3 + PC4, data = pick(everything()))), .groups = "drop")

#### ----------- SUMMARIZE ACROSS ALL DRAWS ----------- ####

coef_draws <- fits %>%
  mutate(tidied = purrr::map(fit, tidy)) %>%
  select(draw, tidied) %>%
  unnest(tidied)

coef_summary <- coef_draws %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  summarise(
    median_estimate    = median(estimate),
    ci_low             = quantile(estimate, 0.025),
    ci_high            = quantile(estimate, 0.975),
    prop_positive      = mean(estimate > 0),
    prop_significant   = mean(p.value < 0.05),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(median_estimate)))

message("Coefficient summary across all posterior draws (delta_PC1 trait-slope model):")
print(coef_summary)

r2_summary <- fits %>%
  mutate(glanced = purrr::map(fit, glance)) %>%
  select(draw, glanced) %>%
  unnest(glanced) %>%
  summarise(median_r2 = median(r.squared), ci_low = quantile(r.squared, 0.025), ci_high = quantile(r.squared, 0.975))

message("\nModel R^2 across all posterior draws:")
print(r2_summary)

#### ----------- SPAGHETTI PLOT: PC2 RELATIONSHIP ACROSS ALL DRAWS ----------- ####

library(ggrepel)

pc2_seq <- seq(min(model_data$PC2), max(model_data$PC2), length.out = 50)
hold_means <- model_data %>% summarise(PC1 = mean(PC1), PC3 = mean(PC3), PC4 = mean(PC4))
pred_grid <- tibble(PC2 = pc2_seq, PC1 = hold_means$PC1, PC3 = hold_means$PC3, PC4 = hold_means$PC4)

full_pred <- fits %>%
  mutate(pred = purrr::map(fit, ~ pred_grid %>% mutate(r3_pc1_pred = predict(.x, newdata = pred_grid)))) %>%
  select(draw, pred) %>%
  unnest(pred)

median_line <- full_pred %>%
  group_by(PC2) %>%
  summarise(median_pred = median(r3_pc1_pred), .groups = "drop")

species_points <- model_data %>%
  group_by(clean.spec, PC2) %>%
  summarise(r3_pc1_median = median(r3_pc1), .groups = "drop")

p_spaghetti <- ggplot() +
  geom_line(data = full_pred, aes(x = PC2, y = r3_pc1_pred, group = draw),
            color = "#2166ac", alpha = 0.05, linewidth = 0.3) +
  geom_line(data = median_line, aes(x = PC2, y = median_pred), color = "black", linewidth = 1) +
  geom_point(data = species_points, aes(x = PC2, y = r3_pc1_median), color = "black", size = 1.3) +
  labs(
    x = "Trait PC2",
    y = "\u0394PC1 (temperature change) random slope\n(trait component, r3_pc1)",
    title = "Posterior uncertainty in the trait PC2 relationship",
    subtitle = paste0(n_distinct(full_pred$draw), " posterior draws (thin lines) + median fit (bold)")
  ) +
  theme_bw(base_size = 11)

if (interactive()) print(p_spaghetti)

ggsave("results/TraitPC2_spaghetti_posterior_delta.tiff", plot = p_spaghetti,
       width = 6.5, height = 5, dpi = 600, device = ragg::agg_tiff)


#########################################################################
# WHICH RAW TRAITS PREDICT THE delta_PC1 (TEMPERATURE CHANGE) RANDOM SLOPE?
#########################################################################
# Same reviewer comment as the original ("which specific traits or trait
# combinations were important"), applied to the delta model. Uses the 7
# raw traits directly (not orthogonal trait PCs) so results can be named
# in prose (e.g. "colony growth rate"), and re-checks where the two
# manuscript-named species (Porites lutea, Orbicella faveolata) sit on
# whichever trait(s) come out strongest -- now with respect to
# SENSITIVITY TO CHANGE rather than absolute switching propensity.
#
# INPUT:  data/shiny_app_data_delta.RDS, data/McWilliam_2018/data/traitbiogeography.csv
# OUTPUT: results/RawTrait_deltaPC1Slope_regression.xlsx
#########################################################################

library(tidyverse)
library(broom)
library(car)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data_delta.RDS")

#### ----------- r3_pc1 DRAWS (delta_PC1 trait component) ----------- ####

r3_pc1 <- app_data$r3_pc1
colnames(r3_pc1) <- app_data$species_levels

#### ----------- RAW TRAIT DATA (unchanged source/cleaning) ----------- ####

trait_vars <- c("cat_growthrate", "cat_skeletaldensity", "cat_corallitesize",
                "cat_colonydiameter", "cat_colonyheight", "cat_SA_vol", "cat_spacesize")
trait_labels <- c(
  cat_growthrate = "Growth rate", cat_skeletaldensity = "Skeletal density",
  cat_corallitesize = "Corallite size", cat_colonydiameter = "Colony diameter",
  cat_colonyheight = "Colony height", cat_SA_vol = "SA:vol ratio", cat_spacesize = "Space size"
)

trait_raw <- read.csv("data/McWilliam_2018/data/traitbiogeography.csv", row.names = 1)

traits_sp <- na.omit(trait_raw[, c(trait_vars, "species")]) %>%
  mutate(clean.spec = gsub(" ", "_", species)) %>%
  mutate(clean.spec = case_when(
    clean.spec == "Montastrea_cavernosa" ~ "Montastraea_cavernosa",
    clean.spec == "Pocillopora_eydouxi"  ~ "Pocillopora_grandis",
    .default = clean.spec
  )) %>%
  mutate(across(all_of(trait_vars), ~ as.numeric(scale(.)))) %>%
  select(clean.spec, all_of(trait_vars))

#### ----------- COLLINEARITY CHECK (once, on point estimates) ----------- ####

species_point_est <- tibble(clean.spec = app_data$species_levels,
                            r3_pc1_median = apply(r3_pc1, 2, median)) %>%
  inner_join(traits_sp, by = "clean.spec")

vif_check <- lm(as.formula(paste("r3_pc1_median ~", paste(trait_vars, collapse = " + "))),
                data = species_point_est)
message("Variance inflation factors (VIF > ~5 is concerning for collinearity):")
print(vif(vif_check))

#### ----------- MERGE: LONG FORMAT (draw x species) ----------- ####

r3_pc1_long <- as_tibble(r3_pc1) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "clean.spec", values_to = "r3_pc1")

model_data <- r3_pc1_long %>% inner_join(traits_sp, by = "clean.spec")

missing_species <- setdiff(app_data$species_levels, unique(model_data$clean.spec))
if (length(missing_species) > 0) {
  message(length(missing_species), " species excluded from model_data (no raw trait data): ",
          paste(missing_species, collapse = ", "))
}

#### ----------- PER-DRAW MULTIPLE REGRESSION ----------- ####

trait_formula <- as.formula(paste("r3_pc1 ~", paste(trait_vars, collapse = " + ")))

fits <- model_data %>%
  group_by(draw) %>%
  summarise(fit = list(lm(trait_formula, data = pick(everything()))), .groups = "drop")

coef_draws <- fits %>%
  mutate(tidied = purrr::map(fit, tidy)) %>%
  select(draw, tidied) %>%
  unnest(tidied)

coef_summary <- coef_draws %>%
  filter(term != "(Intercept)") %>%
  mutate(trait = trait_labels[term]) %>%
  group_by(trait, term) %>%
  summarise(
    median_estimate  = median(estimate),
    ci_low           = quantile(estimate, 0.025),
    ci_high          = quantile(estimate, 0.975),
    prop_positive    = mean(estimate > 0),
    prop_significant = mean(p.value < 0.05),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(median_estimate)))

message("\nRaw trait coefficient summary across all posterior draws (delta_PC1 slope):")
print(coef_summary)

r2_summary <- fits %>%
  mutate(glanced = purrr::map(fit, glance)) %>%
  select(draw, glanced) %>%
  unnest(glanced) %>%
  summarise(median_r2 = median(r.squared), ci_low = quantile(r.squared, 0.025), ci_high = quantile(r.squared, 0.975))

message("\nModel R^2 (all 7 raw traits) across all posterior draws:")
print(r2_summary)

#### ----------- WHERE DO P. LUTEA / O. FAVEOLATA SIT ON THE WINNING TRAIT(S)? ----------- ####

top_traits <- coef_summary %>% slice_max(abs(median_estimate), n = 2) %>% pull(term)

focal_species <- c("Porites_lutea", "Orbicella_faveolata")
missing_focal <- setdiff(focal_species, traits_sp$clean.spec)
if (length(missing_focal) > 0) {
  warning("Requested focal species not found in trait data: ", paste(missing_focal, collapse = ", "))
}

focal_check <- traits_sp %>%
  filter(clean.spec %in% focal_species) %>%
  select(clean.spec, all_of(top_traits))

all_species_trait_range <- traits_sp %>%
  summarise(across(all_of(top_traits), list(min = min, mean = mean, max = max)))

message("\nTop trait(s) by |effect| on sensitivity to environmental CHANGE: ",
        paste(trait_labels[top_traits], collapse = ", "))
message("Porites lutea / Orbicella faveolata values on these traits (standardized, z-scores):")
print(focal_check)
message("\nFor context, full-dataset range on these traits:")
print(all_species_trait_range)

#### ----------- SAVE ----------- ####

writexl::write_xlsx(
  list(coef_summary = coef_summary, r2_summary = r2_summary,
       focal_species_check = focal_check, vif = enframe(vif(vif_check), name = "trait", value = "VIF")),
  "results/RawTrait_deltaPC1Slope_regression.xlsx"
)


#########################################################################
# TRAIT PC2 vs. delta_PC1 (temperature CHANGE) RANDOM SLOPE -- TRAIT COMPONENT
#########################################################################
# Recreates the spaghetti plot showing each species' TRAIT PC2 score
# against the posterior of r3_pc1 -- the TRAIT-covariance component
# (gr(clean.spec_trait, cov = T_mat)) of the delta_PC1 random slope from
# the Dirichlet-partitioned ordinal model. This isolates trait-driven
# variation in temperature-change sensitivity, separate from the
# phylogenetic (r1_pc1) and idiosyncratic/residual (r2_pc1) components.
#
# INPUT:  data/shiny_app_data_delta.RDS
# OUTPUT: results/TraitPC2_vs_deltaPC1_slope.tiff
#########################################################################

library(tidyverse)
library(broom)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data_delta.RDS")

#### ----------- r3_pc1 DRAWS (trait component of delta_PC1 slope) ----------- ####

r3_pc1 <- app_data$r3_pc1
colnames(r3_pc1) <- app_data$species_levels  # species identity is positional -- see prep-script note

#### ----------- TRAIT PCs ----------- ####

trait_pc_cols <- grep("^PC[0-9]+$", colnames(app_data$trait_pca_lookup), value = TRUE)
traits_sp <- app_data$trait_pca_lookup %>%
  select(clean.spec, all_of(trait_pc_cols))

#### ----------- MERGE: LONG FORMAT (draw x species) ----------- ####

r3_pc1_long <- as_tibble(r3_pc1) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "clean.spec", values_to = "r3_pc1")

model_data <- r3_pc1_long %>%
  inner_join(traits_sp, by = "clean.spec")

missing_species <- setdiff(app_data$species_levels, unique(model_data$clean.spec))
if (length(missing_species) > 0) {
  message(length(missing_species), " species excluded (no trait PC data): ",
          paste(missing_species, collapse = ", "))
}

#### ----------- PER-DRAW REGRESSION (all trait PCs, to hold others constant) ----------- ####

fits <- model_data %>%
  group_by(draw) %>%
  summarise(fit = list(lm(r3_pc1 ~ PC1 + PC2 + PC3 + PC4, data = pick(everything()))),
            .groups = "drop")

#### ----------- PREDICTION GRID: VARY PC2, HOLD PC1/PC3/PC4 AT MEAN ----------- ####

pc2_seq <- seq(min(model_data$PC2), max(model_data$PC2), length.out = 50)
hold_means <- model_data %>% summarise(PC1 = mean(PC1), PC3 = mean(PC3), PC4 = mean(PC4))
pred_grid <- tibble(PC2 = pc2_seq, PC1 = hold_means$PC1, PC3 = hold_means$PC3, PC4 = hold_means$PC4)

# Predict from every draw's model -- explicitly qualified purrr::map to
# avoid the maps::map() namespace collision
full_pred <- fits %>%
  mutate(pred = purrr::map(fit, ~ pred_grid %>% mutate(r3_pc1_pred = predict(.x, newdata = pred_grid)))) %>%
  select(draw, pred) %>%
  unnest(pred)

median_line <- full_pred %>%
  group_by(PC2) %>%
  summarise(median_pred = median(r3_pc1_pred), .groups = "drop")

species_points <- model_data %>%
  group_by(clean.spec, PC2) %>%
  summarise(r3_pc1_median = median(r3_pc1), .groups = "drop")

#### ----------- PLOT ----------- ####

p_spaghetti <- ggplot() +
  geom_line(data = full_pred, aes(x = PC2, y = r3_pc1_pred, group = draw),
            color = "#2166ac", alpha = 0.05, linewidth = 0.3) +
  geom_line(data = median_line, aes(x = PC2, y = median_pred), color = "black", linewidth = 1) +
  geom_point(data = species_points, aes(x = PC2, y = r3_pc1_median), color = "black", size = 1.3) +
  labs(
    x = "Trait PC2",
    y = "\u0394PC1 (temperature change) random slope\n(trait component, r3_pc1)",
    title = "Posterior uncertainty in the trait PC2 relationship",
    subtitle = paste0(n_distinct(full_pred$draw), " posterior draws (thin lines) + median fit (bold)")
  ) +
  theme_bw(base_size = 11)

if (interactive()) print(p_spaghetti)

ggsave("results/TraitPC2_vs_deltaPC1_slope.tiff", plot = p_spaghetti,
       width = 6.5, height = 5, dpi = 600, device = ragg::agg_tiff)


#########################################################################
# TRAIT PC2 vs. delta_PC1 TRAIT-COMPONENT SLOPE -- SPECIES-LEVEL MEDIAN VERSION
#########################################################################
# Simpler alternative to the per-draw spaghetti plot: collapses each
# species' r3_pc1 posterior to its MEDIAN first, then fits a single lm()
# of that species-level point estimate against trait PC2. geom_smooth()
# gives the standard OLS 95% CI ribbon on this species-level relationship.
#
# Three species called out with triangle markers + labels, matching your
# reference plot (Acropora hyacinthus, Porites lutea, Montipora capitata)
# -- swap FOCAL_SPECIES below if you want different ones highlighted.
#
# INPUT:  data/shiny_app_data_delta.RDS
# OUTPUT: results/TraitPC2_vs_deltaPC1_slope_median.tiff
#########################################################################

library(tidyverse)
library(ggrepel)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data_delta.RDS")

FOCAL_SPECIES <- c("Acropora_hyacinthus", "Porites_lutea", "Montipora_capitata")

#### ----------- r3_pc1 DRAWS -> SPECIES-LEVEL MEDIAN ----------- ####

r3_pc1 <- app_data$r3_pc3
colnames(r3_pc1) <- app_data$species_levels

species_median <- tibble(
  clean.spec    = app_data$species_levels,
  r3_pc1_median = apply(r3_pc1, 2, median)
)

#### ----------- MERGE WITH TRAIT PC2 ----------- ####

traits_sp <- app_data$trait_pca_lookup %>% select(clean.spec, PC1:PC4)

model_data <- species_median %>%
  inner_join(traits_sp, by = "clean.spec")

missing_species <- setdiff(app_data$species_levels, model_data$clean.spec)
if (length(missing_species) > 0) {
  message(length(missing_species), " species excluded (no trait PC data): ",
          paste(missing_species, collapse = ", "))
}

#### ----------- SINGLE LM ON SPECIES-LEVEL MEDIANS ----------- ####

fit <- lm(r3_pc1_median ~ PC1 + PC2 + PC3 + PC4, data = model_data)
print(summary(fit))

#### ----------- FLAG FOCAL SPECIES FOR LABELING ----------- ####

model_data <- model_data %>%
  mutate(
    is_focal   = clean.spec %in% FOCAL_SPECIES,
    label_text = ifelse(is_focal, gsub("_", " ", clean.spec), NA_character_)
  )

missing_focal <- setdiff(FOCAL_SPECIES, model_data$clean.spec)
if (length(missing_focal) > 0) {
  warning("Focal species not found in model_data: ", paste(missing_focal, collapse = ", "))
}

#### ----------- PLOT ----------- ####

p_median <- ggplot(model_data, aes(x = PC1, y = r3_pc1_median)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = F, color = "#d7191c", fill = "#d7191c", alpha = 0.15, linewidth = 0.8) +
  geom_point(data = filter(model_data, !is_focal),
             color = "#2166ac", alpha = 0.6, size = 2) +
  geom_point(data = filter(model_data, is_focal),
             shape = 17, color = "black", size = 3) +
  geom_text_repel(aes(label = label_text), fontface = "italic", size = 3.2,
                  na.rm = TRUE, seed = 1, box.padding = 0.6) +
  labs(
    x = "Trait PC1 (Growth Rate)",
    y = "\u0394P(Switch) ~ \u0394Temp."
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

if (interactive()) print(p_median)

ggsave("results/TraitPC2_vs_deltaPC1_slope_median.tiff", plot = p_median,
       width = 4.2, height = 3.5, dpi = 600, device = ragg::agg_tiff)