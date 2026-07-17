#########################################################################
# EXTRACT r3_pc1 DRAWS + MERGE WITH TRAIT PCs
#########################################################################
# Builds a tidy, long-format (draw x species) table with the trait-
# covariance component of the PC1 (temperature) random slope alongside
# each species' TRAIT PC scores (not raw traits) -- ready for a per-draw
# multiple regression. Trait PCs are used here rather than the 7 raw
# traits specifically because they're orthogonal by construction, so a
# multiple regression on them isn't confounded by collinearity between
# traits (e.g. colony height/diameter are likely correlated with each
# other in the raw trait space).
#
# Uses shiny_app_data.RDS rather than re-deriving from
# stan_dirichlet_explicit.RDS directly: r3_pc1 in that file is already
# correctly extracted and species-indexed (verified against the model's
# real design matrix via standata() in 05_prep_shiny_data.R), so this
# avoids re-solving that same species-ordering problem here.
#
# INPUT:  data/shiny_app_data.RDS
# OUTPUT: model_data (in-memory tibble; not saved to file by default --
#         add a saveRDS()/write_csv() at the bottom if you want it on disk)
#########################################################################

library(tidyverse)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data.RDS")
darling <- read.csv("~/Documents/Coral/LongitudinalZoox/data/Darling_coral_species_traits.csv")
darling = darling %>% select(Species, Life.history) %>% mutate(clean.spec = gsub(" ", "_", Species))
#### ----------- r3_pc1 DRAWS (n_draws x n_species matrix) ----------- ####

r3_pc1 <- app_data$r3_pc1
# r3_pc1's column names as stored are raw Stan parameter names
# (e.g. "r_3_4[12]"), not species -- species identity is only carried by
# POSITION, matching app_data$species_levels. Set real names explicitly
# before reshaping, or a pivot_longer() below would silently mislabel
# every species.
colnames(r3_pc1) <- app_data$species_levels

#### ----------- TRAIT PCs (already computed, orthogonal by construction) ----------- ####

trait_pc_cols <- grep("^PC[0-9]+$", colnames(app_data$trait_pca_lookup), value = TRUE)
app_data$
message("Trait PCs available: ", paste(trait_pc_cols, collapse = ", "))

traits_sp <- app_data$trait_pca_lookup %>%
  select(clean.spec, all_of(trait_pc_cols))

#### ----------- MERGE: LONG FORMAT, ONE ROW PER (DRAW x SPECIES) ----------- ####

r3_pc1_long <- as_tibble(r3_pc1) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "clean.spec", values_to = "r3_pc1")

model_data <- r3_pc1_long %>%
  inner_join(traits_sp, by = "clean.spec") %>% left_join(darling)
model_data %>% group_by(draw, Life.history) %>%
  drop_na() %>%
  mutate(mean_r3 = mean(r3_pc1)) %>%
  ggplot(aes(x=Life.history,  y=mean_r3)) + geom_boxplot()

missing_species <- setdiff(app_data$species_levels, unique(model_data$clean.spec))
if (length(missing_species) > 0) {
  message(length(missing_species), " species excluded from model_data (no trait PC data): ",
          paste(missing_species, collapse = ", "))
}

message("model_data ready: ", n_distinct(model_data$draw), " draws x ",
        n_distinct(model_data$clean.spec), " species with trait PC data (",
        nrow(model_data), " total rows).")

# model_data columns: draw, clean.spec, r3_pc1, PC1, PC2, PC3, PC4 (however
# many trait PCs were retained)

#### ----------- PER-DRAW REGRESSION ----------- ####

library(broom)

fits <- model_data %>%
  group_by(draw) %>%
  summarise(fit = list(lm(r3_pc1 ~ PC1 + PC2 + PC3 + PC4, data = pick(everything()))), .groups = "drop")

#### ----------- SUMMARIZE ACROSS ALL DRAWS ----------- ####

# One row per (draw x coefficient) -- broom::tidy() extracts estimate,
# std.error, statistic, p.value from each fitted lm
coef_draws <- fits %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(draw, tidied) %>%
  unnest(tidied)

coef_summary <- coef_draws %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  summarise(
    median_estimate    = median(estimate),
    ci_low             = quantile(estimate, 0.025),
    ci_high            = quantile(estimate, 0.975),
    prop_positive      = mean(estimate > 0),   # analog of a Bayesian "probability of direction"
    prop_significant   = mean(p.value < 0.05), # how often this term reached p<0.05 across draws --
    # descriptive consistency, not itself a p-value
    .groups = "drop"
  ) %>%
  arrange(desc(abs(median_estimate)))

message("Coefficient summary across all posterior draws:")
print(coef_summary)

# Overall model fit (R^2) across draws -- how much of the trait-covariance
# slope's variance the 4 trait PCs collectively explain
r2_summary <- fits %>%
  mutate(glanced = map(fit, glance)) %>%
  select(draw, glanced) %>%
  unnest(glanced) %>%
  summarise(median_r2 = median(r.squared), ci_low = quantile(r.squared, 0.025), ci_high = quantile(r.squared, 0.975))

message("\nModel R^2 across all posterior draws:")
print(r2_summary)

#### ----------- SPAGHETTI PLOT: PC2 RELATIONSHIP ACROSS ALL DRAWS ----------- ####
# Predicts r3_pc1 vs. trait PC2 from EVERY draw's fitted model, holding
# PC1/PC3/PC4 at their mean -- each draw's prediction becomes one line.
# Plotting all ~4000 draws directly (heavily transparent) shows the full
# posterior spread; a random subsample would look the same but render
# faster if 4000 lines is too slow/heavy for your machine.

library(ggrepel)

# Prediction grid: vary PC2 over its observed range, hold others at mean
pc2_seq <- seq(min(model_data$PC2), max(model_data$PC2), length.out = 50)
hold_means <- model_data %>% summarise(PC1 = mean(PC1), PC3 = mean(PC3), PC4 = mean(PC4))
pred_grid <- tibble(PC2 = pc2_seq, PC1 = hold_means$PC1, PC3 = hold_means$PC3, PC4 = hold_means$PC4)

# Predict from every draw's model -- this is the expensive step (4000
# predict() calls); a few seconds is normal
full_pred <- fits %>%
  mutate(pred = map(fit, ~ pred_grid %>% mutate(r3_pc1_pred = predict(.x, newdata = pred_grid)))) %>%
  select(draw, pred) %>%
  unnest(pred)

# Median + 95% interval line, summarized across all draws at each PC2 value
median_line <- full_pred %>%
  group_by(PC2) %>%
  summarise(median_pred = median(r3_pc1_pred), .groups = "drop")

# Observed per-species point estimates (posterior median r3_pc1), for
# context against the fitted lines
species_points <- model_data %>%
  group_by(clean.spec, PC2) %>%
  summarise(r3_pc1_median = median(r3_pc1), .groups = "drop")

# If you only want a subsample of draws plotted (faster render, same
# visual impression), uncomment:
# set.seed(1)
# sample_draws <- sample(unique(full_pred$draw), 300)
# full_pred <- full_pred %>% filter(draw %in% sample_draws)

p_spaghetti <- ggplot() +
  geom_line(data = full_pred, aes(x = PC2, y = r3_pc1_pred, group = draw),
            color = "#2166ac", alpha = 0.05, linewidth = 0.3) +
  geom_line(data = median_line, aes(x = PC2, y = median_pred), color = "black", linewidth = 1) +
  geom_point(data = species_points, aes(x = PC2, y = r3_pc1_median), color = "black", size = 1.3) +
  labs(
    x = "Trait PC2",
    y = "PC1 (temperature) random slope\n(trait component, r3_pc1)",
    title = "Posterior uncertainty in the trait PC2 relationship",
    subtitle = paste0(n_distinct(full_pred$draw), " posterior draws (thin lines) + median fit (bold)")
  ) +
  theme_bw(base_size = 11)

if (interactive()) print(p_spaghetti)
hist(species_points$PC2)

ggsave("results/TraitPC2_spaghetti_posterior.tiff", plot = p_spaghetti, width = 6.5, height = 5, dpi = 600)




#########################################################################
# WHICH RAW TRAITS PREDICT THE PC1 (TEMPERATURE) RANDOM SLOPE?
#########################################################################
# Addresses reviewer comment: "It would seem worthwhile to say more about
# which specific traits or trait combinations were important." Replaces
# the earlier trait-PC-based regression (12_extract_r3pc1_merge_traits.R)
# with the 7 RAW traits directly -- PCs are orthogonal and were the right
# choice for a clean regression, but "PC2" isn't something you can name
# in a discussion section the way "colony growth rate" is.
#
# Also specifically checks where Porites lutea and Orbicella faveolata
# (the two species already named in the manuscript text) sit on whichever
# raw trait(s) come out strongest -- since confirming they're both
# extreme on the same trait axis is what would let the discussion section
# say something concrete rather than "host traits... may constrain
# symbiotic flexibility."
#
# INPUT:  data/shiny_app_data.RDS, data/McWilliam_2018/data/traitbiogeography.csv
# OUTPUT: results/RawTrait_PC1Slope_regression.xlsx
#########################################################################

library(tidyverse)
library(broom)
library(car)   # for VIF collinearity check -- install.packages("car") if needed

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data.RDS")

#### ----------- r3_pc1 DRAWS ----------- ####

r3_pc1 <- app_data$r3_pc1
colnames(r3_pc1) <- app_data$species_levels  # see note in script 12 -- stored colnames are Stan
# parameter names, not species; species identity is
# positional, matching species_levels

#### ----------- RAW TRAIT DATA (same source/cleaning as elsewhere in this project) ----------- ####

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
# Raw traits, unlike the orthogonal trait PCs, can be correlated with
# each other (e.g. colony height/diameter) -- check VIF before trusting
# individual coefficients. VIF > ~5 is generally considered concerning.

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
  mutate(tidied = map(fit, tidy)) %>%
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

message("\nRaw trait coefficient summary across all posterior draws:")
print(coef_summary)

r2_summary <- fits %>%
  mutate(glanced = map(fit, glance)) %>%
  select(draw, glanced) %>%
  unnest(glanced) %>%
  summarise(median_r2 = median(r.squared), ci_low = quantile(r.squared, 0.025), ci_high = quantile(r.squared, 0.975))

message("\nModel R^2 (all 7 raw traits) across all posterior draws:")
print(r2_summary)

#### ----------- WHERE DO P. LUTEA / O. FAVEOLATA SIT ON THE WINNING TRAIT(S)? ----------- ####
# Checks the top 2 traits by |median_estimate| against the two species
# already named in the manuscript text.

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

message("\nTop trait(s) by |effect|: ", paste(trait_labels[top_traits], collapse = ", "))
message("Porites lutea / Orbicella faveolata values on these traits (standardized, z-scores):")
print(focal_check)
message("\nFor context, full-dataset range on these traits:")
print(all_species_trait_range)

#### ----------- SAVE ----------- ####

writexl::write_xlsx(
  list(coef_summary = coef_summary, r2_summary = r2_summary,
       focal_species_check = focal_check, vif = enframe(vif(vif_check), name = "trait", value = "VIF")),
  "results/RawTrait_PC1Slope_regression.xlsx"
)



#########################################################################
# STATISTICAL TEST OF PHYLOGENETIC/TRAIT CLUSTERING (Figure 3b)
#########################################################################
# Addresses reviewer comment: "It would be stronger to report a
# statistical test of this clustering" -- for species-level predicted
# switching intensity at median environmental conditions, mapped onto
# both the phylogenetic tree and the trait dendrogram.
#
# Uses Blomberg's K and Pagel's lambda (phytools::phylosig), which work
# identically on a true phylogeny (tree_final) and on the trait
# dendrogram (trait_tree -- itself a proper `phylo` object built via
# as.phylo(hclust(...)) in 05_prep_shiny_data.R) -- same test, same
# interpretation, applied to both axes shown in Fig 3b.
#
# "Median environmental conditions" here = PC1 = PC2 = PC3 = 0 (the
# environmental PCA's dataset average, since these are mean-centered by
# construction), applied as ONE SHARED scenario across all species, each
# at its own default/most-common genus. If Figure 3b was actually built
# some other way (e.g. species held at their own historical mean rather
# than the global average), adjust the pc1/pc2/pc3 values below to match.
#
# INPUT:  data/shiny_app_data.RDS
# OUTPUT: results/PhyloTraitSignal_SwitchingIntensity.xlsx
#########################################################################

library(tidyverse)
library(phytools)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

app_data <- readRDS("data/shiny_app_data.RDS")

#### ----------- SELF-CONTAINED COPY OF predict_switching() FROM app.R ----------- ####
# Duplicated rather than sourced, since app.R launches the Shiny app on
# source() (shinyApp(ui, server) at the bottom) -- keep this in sync with
# app.R if that function's logic changes there.

predict_switching <- function(pc1, pc2, pc3, species, app_data, delta_days = 365, genus_choice = "default") {
  sp_idx <- match(species, app_data$species_levels)
  valid  <- !is.na(sp_idx)
  species <- species[valid]; sp_idx <- sp_idx[valid]
  n_draws <- nrow(app_data$b_mat)
  if (length(pc1) == 1) pc1 <- rep(pc1, length(species))
  if (length(pc2) == 1) pc2 <- rep(pc2, length(species))
  if (length(pc3) == 1) pc3 <- rep(pc3, length(species))
  
  log_dt_offset <- log(delta_days) - app_data$mean_log_delta_t
  b_log_dt <- if ("log_delta_t" %in% colnames(app_data$b_mat)) app_data$b_mat[, "log_delta_t"] else 0
  
  genus_per_species <- if (identical(genus_choice, "default")) {
    lookup <- app_data$species_default_genus
    lookup$default_genus[match(species, lookup$clean.spec)]
  } else {
    rep(genus_choice, length(species))
  }
  
  genus_dummy_names <- names(app_data$genus_col_means)
  genus_contribution <- function(genus) {
    raw_dummy <- setNames(rep(0, length(genus_dummy_names)), genus_dummy_names)
    target_col <- paste0("last_genus", genus)
    if (target_col %in% genus_dummy_names) raw_dummy[target_col] <- 1
    centered <- raw_dummy - app_data$genus_col_means[genus_dummy_names]
    as.numeric(app_data$b_mat[, genus_dummy_names, drop = FALSE] %*% centered)
  }
  
  purrr::map_dfr(seq_along(species), function(i) {
    sp <- species[i]; si <- sp_idx[i]
    b_genus <- genus_contribution(genus_per_species[i])
    mu <- app_data$r1_int[, si] + app_data$r1_pc1[, si] * pc1[i] +
      app_data$r1_pc2[, si] * pc2[i] + app_data$r1_pc3[, si] * pc3[i] +
      app_data$r2_int[, si] + app_data$r2_pc1[, si] * pc1[i] +
      app_data$r2_pc2[, si] * pc2[i] + app_data$r2_pc3[, si] * pc3[i] +
      app_data$r3_int[, si] + app_data$r3_pc1[, si] * pc1[i] +
      app_data$r3_pc2[, si] * pc2[i] + app_data$r3_pc3[, si] * pc3[i] +
      b_log_dt * log_dt_offset + b_genus
    thres <- app_data$int_mat
    cum_p <- plogis(thres - matrix(mu, nrow = n_draws, ncol = ncol(thres)))
    cat_p <- cbind(cum_p, 1) - cbind(0, cum_p)
    e_switch <- as.numeric(cat_p %*% app_data$category_values)
    tibble(clean.spec = sp, draw = seq_len(n_draws), e_switch = e_switch)
  })
}

#### ----------- PREDICT SWITCHING INTENSITY AT MEDIAN CONDITIONS ----------- ####

species <- app_data$species_levels
pred <- predict_switching(pc1 = 0, pc2 = 0, pc3 = 0, species = species,
                          app_data = app_data, delta_days = 365, genus_choice = "default")

species_switching <- pred %>%
  group_by(clean.spec) %>%
  summarise(e_switch_median = median(e_switch), .groups = "drop")

#### ----------- PHYLOGENETIC SIGNAL ----------- ####

phylo_vec <- setNames(species_switching$e_switch_median, species_switching$clean.spec)
tree_phylo <- keep.tip(app_data$tree_final, intersect(names(phylo_vec), app_data$tree_final$tip.label))
phylo_vec  <- phylo_vec[tree_phylo$tip.label]

missing_phylo <- setdiff(species_switching$clean.spec, tree_phylo$tip.label)
if (length(missing_phylo) > 0) {
  message(length(missing_phylo), " species excluded from phylogenetic signal test (not in tree_final): ",
          paste(missing_phylo, collapse = ", "))
}

phylo_K      <- phylosig(tree_phylo, phylo_vec, method = "K", test = TRUE, nsim = 999)
phylo_lambda <- phylosig(tree_phylo, phylo_vec, method = "lambda", test = TRUE)

message("\n=== PHYLOGENETIC SIGNAL (switching intensity @ median conditions) ===")
message("Blomberg's K = ", round(phylo_K$K, 3), ", p = ", signif(phylo_K$P, 3), " (n = ", length(phylo_vec), " species)")
message("Pagel's lambda = ", round(phylo_lambda$lambda, 3), ", p = ", signif(phylo_lambda$P, 3))

#### ----------- TRAIT SIGNAL (same test, trait dendrogram instead of phylogeny) ----------- ####

trait_vec  <- setNames(species_switching$e_switch_median, species_switching$clean.spec)
tree_trait <- keep.tip(app_data$trait_tree, intersect(names(trait_vec), app_data$trait_tree$tip.label))
trait_vec  <- trait_vec[tree_trait$tip.label]

missing_trait <- setdiff(species_switching$clean.spec, tree_trait$tip.label)
if (length(missing_trait) > 0) {
  message(length(missing_trait), " species excluded from trait signal test (not in trait_tree): ",
          paste(missing_trait, collapse = ", "))
}

trait_K      <- phylosig(tree_trait, trait_vec, method = "K", test = TRUE, nsim = 999)
trait_lambda <- phylosig(tree_trait, trait_vec, method = "lambda", test = TRUE)

message("\n=== TRAIT SIGNAL (switching intensity @ median conditions) ===")
message("Blomberg's K = ", round(trait_K$K, 3), ", p = ", signif(trait_K$P, 3), " (n = ", length(trait_vec), " species)")
message("Pagel's lambda = ", round(trait_lambda$lambda, 3), ", p = ", signif(trait_lambda$P, 3))

#### ----------- SAVE RESULTS ----------- ####

signal_results <- tibble(
  axis      = c("Phylogeny", "Phylogeny", "Trait dendrogram", "Trait dendrogram"),
  statistic = c("Blomberg's K", "Pagel's lambda", "Blomberg's K", "Pagel's lambda"),
  estimate  = c(phylo_K$K, phylo_lambda$lambda, trait_K$K, trait_lambda$lambda),
  p_value   = c(phylo_K$P, phylo_lambda$P, trait_K$P, trait_lambda$P),
  n_species = c(length(phylo_vec), length(phylo_vec), length(trait_vec), length(trait_vec))
)

signal_results
writexl::write_xlsx(list(signal_test = signal_results, species_switching = species_switching),
                    "results/PhyloTraitSignal_SwitchingIntensity.xlsx")